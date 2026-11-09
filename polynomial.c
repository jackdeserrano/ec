#include "field.c"
#include <stdlib.h>

typedef struct {
    fe *coes;
    size_t deg;
} poly;

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
    for (int i = 0; i < a->deg + 1; ++i) {
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

    free(r->coes);
    
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
        
        free(r->coes);
        
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

    free(r->coes);
    
    r->coes = coes;
    r->deg = deg;

    poly_ensure_deg(r);

    return 1;
}

int poly_mul_fe(poly *r, const poly *a, const fe *b) {
    fe *coes = calloc(a->deg + 1, sizeof(fe));

    if (!coes) {
        return 0;
    }

    for (uint64_t i = 0; i < a->deg + 1; ++i) {
        fe_mul(coes + i, a->coes + i, b);
    }

    free(r->coes);
    
    r->coes = coes;
    r->deg = a->deg;

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
