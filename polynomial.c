#include "group.c"
#include <stdlib.h>

typedef struct {
    fe *coes;
    size_t deg;
} poly;

// requires r is unitialized
int poly_init(poly *r, size_t deg) {
    r->coes = calloc(deg + 1, sizeof(fe));

    if (!r->coes) {
        return 0;
    }

    r->deg = deg;
    return 1;
}

void poly_free(poly *r) {
    free(r->coes);
    
    r->coes = 0;
    r->deg = 0;
}

void poly_ensure_deg(poly *r) {
    while (r->coes[r->deg].n == 0 && r->deg > 0) {
        --r->deg;
    }

    fe *temp = realloc(r->coes, (r->deg + 1) * sizeof(fe));

    if (temp) {
        r->coes = temp;
    }
}

int poly_is_zero(const poly *a) {
    for (uint64_t i = 0; i < a->deg + 1; ++i) {
        if (a->coes[i].n != 0) {
            return 0;
        }
    }
    return 1;
}

int poly_add(poly *r, const poly *a, const poly *b) {
    size_t deg = a->deg > b->deg ? a->deg : b->deg;

    fe *coes = calloc(deg + 1, sizeof(fe));

    if (!coes) {
        return 0;
    }

    for (uint64_t i = 0; i < deg + 1; ++i) {
        if (i < a->deg + 1 && i < b->deg + 1) {
            fe_add(coes + i, a->coes + i, b->coes + i);
        } else if (i < a->deg + 1) {
            coes[i] = a->coes[i];
        } else {
            coes[i] = b->coes[i];
        }
    }

    poly_free(r);
    
    r->coes = coes;
    r->deg = deg;

    poly_ensure_deg(r);

    return 1;
}

int poly_mul(poly *r, const poly *a, const poly *b) {
    if (poly_is_zero(a) || poly_is_zero(b)) {
        fe *coes = calloc(1, sizeof(fe));
        
        if (!coes) {
            return 0;
        }
        
        poly_free(r);
        
        r->coes = coes;
        r->deg = 0;
        
        return 1;
    }
    
    size_t deg = a->deg + b->deg; // can overflow

    fe *coes = calloc(deg + 1, sizeof(fe));

    if (!coes) {
        return 0;
    }

    for (uint64_t i = 0; i < a->deg + 1; ++i) {
        for (uint64_t j = 0; j < b->deg + 1; ++j) {
            fe c;
            fe_mul(&c, a->coes + i, b->coes + j);
            fe_add(coes + i + j, coes + i + j, &c);
        }
    }

    poly_free(r);
    
    r->coes = coes;
    r->deg = deg;

    poly_ensure_deg(r);

    return 1;
}

int poly_mul_fe(poly *r, const poly *a, const fe *b) {
    fe *coes = calloc(a->deg + 1, sizeof(fe));
    size_t deg = a->deg;

    if (!coes) {
        return 0;
    }

    for (uint64_t i = 0; i < a->deg + 1; ++i) {
        fe_mul(coes + i, a->coes + i, b);
    }

    poly_free(r);
    
    r->coes = coes;
    r->deg = deg;

    poly_ensure_deg(r);

    return 1;
}

int poly_mul_int(poly *r, const poly *a, const uint64_t b) {
    fe *coes = calloc(a->deg + 1, sizeof(fe));
    size_t deg = a->deg;

    if (!coes) {
        return 0;
    }

    for (uint64_t i = 0; i < a->deg + 1; ++i) {
        fe_mul_int(coes + i, a->coes + i, b);
    }

    poly_free(r);
    
    r->coes = coes;
    r->deg = deg;

    poly_ensure_deg(r);

    return 1;
}

int poly_sub(poly *r, const poly *a, const poly *b) {
    fe c = {p - 1};
    
    poly temp;
    
    if (!poly_init(&temp, 0)) {
        return 0;
    }
    
    if (!poly_mul_fe(&temp, b, &c)) {
        poly_free(&temp);
        return 0;
    }
    
    int success = poly_add(r, a, &temp);
    
    poly_free(&temp);
    
    return success;
}

int poly_divrem(poly *q, poly *r, const poly *a, const poly *b) {
    assert(q != r);
    
    if (poly_is_zero(b)) {
        return 0;
    }
    
    poly q_temp, r_temp;
    
    size_t q_deg = a->deg > b->deg ? a->deg - b->deg : 0;
    
    if (!poly_init(&q_temp, q_deg)) {
        return 0;
    }
    
    if (!poly_init(&r_temp, a->deg)) {
        poly_free(&q_temp);
        return 0;
    }
    
    for (uint64_t i = 0; i < a->deg + 1; ++i) {
        r_temp.coes[i] = a->coes[i];
    }
    
    while (r_temp.deg >= b->deg && !(poly_is_zero(&r_temp))) {
        size_t expo = r_temp.deg - b->deg;
        
        fe c;
        fe_inv(&c, b->coes + b->deg);
        fe_mul(&c, &r_temp.coes[r_temp.deg], &c);
        
        q_temp.coes[expo] = c;
        
        for (uint64_t i = 0; i < b->deg + 1; ++i) {
            fe d;
            fe_mul(&d, &c, b->coes + i);
            fe_sub(&r_temp.coes[i + expo], &r_temp.coes[i + expo], &d);
        }
        
        poly_ensure_deg(&r_temp);
    }
    
    poly_ensure_deg(&q_temp);
    
    poly_free(q);
    poly_free(r);
    
    *q = q_temp;
    *r = r_temp;
    
    return 1;
}

int poly_rem(poly *r, const poly *a, const poly *b) {
    poly q;
    
    if (!poly_init(&q, 0)) {
        return 0;
    }
    
    if (!poly_divrem(&q, r, a, b)) {
        poly_free(&q);
        return 0;
    }
    
    poly_free(&q);
    
    return 1;
}

int poly_powmod(poly *r, const poly *a, uint64_t n, const poly *m) {
    poly r_temp, b;
    
    if (!poly_init(&r_temp, 0)) {
        return 0;
    }
    
    if (!poly_init(&b, 0)) {
        poly_free(&r_temp);
        return 0;
    }
    
    r_temp.coes[0].n = 1;
    
    if (!poly_rem(&b, a, m)) {
        poly_free(&r_temp);
        poly_free(&b);
        return 0;
    }
    
    while (n > 0) {
        
        if (n & 1) {
            
            if (!poly_mul(&r_temp, &r_temp, &b) || !poly_rem(&r_temp, &r_temp, m)) {
                poly_free(&r_temp);
                poly_free(&b);
            }

        }
        
        if (!poly_mul(&b, &b, &b) || !poly_rem(&b, &b, m)) {
            poly_free(&r_temp);
            poly_free(&b);
        }
        
        n >>= 1;
        
    }
    
    poly_free(r);
    *r = r_temp;
    
    poly_free(&b);
    
    return 1;
}

int poly_div_polys(poly **r, uint64_t n) {
    poly *r_temp;
    
    if ((__uint128_t)n * n >= (__uint128_t)2 * (UINT64_MAX - 1)) {
        
        return 0;
        
    }
    
    r_temp = calloc(n + 1, sizeof(poly));
    
    if (!r_temp) {
        
        return 0;
        
    }
    
    poly f = {0};
    
    if (!poly_init(&f, 3)) {
        
        goto memory;
        
    }
    
    fe a = {A};
    fe b = {B};
    
    f.coes[0] = b;
    f.coes[1] = a;
    f.coes[3].n = 1;
    
    if (!poly_mul(&f, &f, &f) || !poly_mul_int(&f, &f, 16)) {
        
        goto memory;
        
    }
    
    poly p1 = {0};
    poly p2 = {0};
    
    fe a2, a3, b2, ab;
    
    for (uint64_t i = 0; i < n + 1; ++i) {
        
        if (!poly_init(&r_temp[i], i * i / 2)) {
            
            goto memory;
            
        }
        
        if (i == 1) {
            
            r_temp[i].coes[0].n = 1;
            
        }
        
        if (i == 2) {
            
            r_temp[i].coes[0].n = 1;
            
        }
        
        if (i == 3) {
            
            fe_sqr(&a2, &a);
            
            fe_neg(&r_temp[i].coes[0], &a2);
            fe_mul_int(&r_temp[i].coes[1], &b, 12);
            fe_mul_int(&r_temp[i].coes[2], &a, 6);
            
            r_temp[i].coes[4].n = 3;
            
        }
        
        if (i == 4) {
            
            fe_sqr(&b2, &b);
            
            fe_mul(&a3, &a2, &a);
            
            fe_mul(&ab, &a, &b);
            
            fe_mul_int(&r_temp[i].coes[0], &b2, 8);
            fe_add(&r_temp[i].coes[0], &r_temp[i].coes[0], &a3);
            fe_mul_int(&r_temp[i].coes[0], &r_temp[i].coes[0], 2);
            fe_neg(&r_temp[i].coes[0], &r_temp[i].coes[0]);
            
            fe_mul_int(&r_temp[i].coes[1], &ab, 8);
            fe_neg(&r_temp[i].coes[1], &r_temp[i].coes[1]);
            
            fe_mul_int(&r_temp[i].coes[2], &a2, 10);
            fe_neg(&r_temp[i].coes[2], &r_temp[i].coes[2]);
            
            fe_mul_int(&r_temp[i].coes[3], &b, 40);
            
            fe_mul_int(&r_temp[i].coes[4], &a, 10);
            
            r_temp[i].coes[6].n = 2;
            
        }
        
        if (i > 4) {
            
            if (!poly_init(&p1, 0) || !poly_init(&p2, 0)) {
                
                goto memory;
                
            }
            
            if (i & 1) {
                
                uint64_t k = (i - 1) / 2;
                
                if (!poly_mul(&p1, &r_temp[k], &r_temp[k]) || !poly_mul(&p1, &p1, &r_temp[k]) || !poly_mul(&p1, &p1, &r_temp[k + 2])) {
                    
                    goto memory;
                    
                }
                
                if (!poly_mul(&p2, &r_temp[k + 1], &r_temp[k + 1]) || !poly_mul(&p2, &p2, &r_temp[k + 1]) || !poly_mul(&p2, &p2, &r_temp[k - 1])) {
                    
                    goto memory;
                    
                }
                
                
                if (i & 2) {
                    
                    if (!poly_mul(&p2, &p2, &f) || !poly_sub(&r_temp[i], &p1, &p2)) {
                        
                        goto memory;
                        
                    }
                    
                } else {
                    
                    if (!poly_mul(&p1, &p1, &f) || !poly_sub(&r_temp[i], &p1, &p2)) {
                        
                        goto memory;
                        
                    }
                    
                }
                
            } else {
                
                uint64_t k = i / 2;
                
                if (!poly_mul(&p1, &r_temp[k - 1], &r_temp[k - 1]) || !poly_mul(&p1, &p1, &r_temp[k + 2])) {
                    
                    goto memory;
                    
                }
                
                if (!poly_mul(&p2, &r_temp[k + 1], &r_temp[k + 1]) || !poly_mul(&p2, &p2, &r_temp[k - 2])) {
                    
                    goto memory;
                    
                }
                
                if (!poly_sub(&r_temp[i], &p1, &p2) || !poly_mul(&r_temp[i], &r_temp[i], &r_temp[k])) {
                    
                    goto memory;
                    
                }
                
            }
            
        }
        
        poly_ensure_deg(&r_temp[i]);
        
        poly_free(&p1);
        poly_free(&p2);
        
    }
    
    poly_free(&f);
    
    *r = r_temp;
    
    return 1;
    
memory:
    
    for (uint64_t i = 0; i < n + 1; ++i) {
        
        poly_free(&r_temp[i]);
        
    }
    
    free(r_temp);
    
    poly_free(&f);
    poly_free(&p1);
    poly_free(&p2);
    
    return 0;
    
}
