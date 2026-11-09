#include <stdlib.h>
#include "field.c"

uint64_t A = 0;
uint64_t B = 7;

/* the elliptic curve y^2 = x^3 + 7 has bad reduction at
 * the primes p = 2, 3, 7. when it has good reduction at
 * p and p is 2 mod 3, the curve has p + 1 points. */

typedef struct {
    fe x;
    fe y;
    int infinity;
} ge;

const ge inf = {{0}, {0}, 1};

int ge_valid(const ge *a) {
    if (a->infinity) {
        return 1;
    }

    fe y2,x3;

    fe_sqr(&y2, &(a->y));

    fe_sqr(&x3, &(a->x));
    fe_mul(&x3, &x3, &(a->x));
    fe_add_int(&x3, &x3, B);

    return fe_equal(&y2, &x3);
}

void ge_neg(ge *r, const ge *a) {
    assert(ge_valid(a));
    
    *r = *a;
    fe_neg(&(r->y), &(r->y));
}

void ge_double(ge *r, const ge *a) {
    assert(ge_valid(a));
    
    ge c = *a;
    fe num, den, lambda;
    
    if ((c.y).n == 0) {
        r->infinity = 1;
    } else {
        
        r->infinity = 0;
        
        fe_sqr(&num, &(c.x));
        fe_mul_int(&num, &num, 3);
        fe_mul_int(&den, &(c.y), 2);
        
        fe_inv(&lambda, &den);
        fe_mul(&lambda, &lambda, &num);
        
        fe_sqr(&(r->x), &lambda);
        fe_sub(&(r->x), &(r->x), &(c.x));
        fe_sub(&(r->x), &(r->x), &(c.x));
        
        fe_sub(&(r->y), &(c.x), &(r->x));
        fe_mul(&(r->y), &(r->y), &lambda);
        fe_sub(&(r->y), &(r->y), &(c.y));
        
    }
}

void ge_add(ge *r, const ge *a, const ge *b) {
    assert(ge_valid(a));
    assert(ge_valid(b));
    
    fe num, den, lambda;
    
    if (a->infinity) {
        *r = *b;
    } else if (b->infinity) {
        *r = *a;
    } else if ((a->x).n == (b->x).n) {
        
        if ((a->y).n != (b->y).n) {
            r->infinity = 1;
        } else {
            ge_double(r, a);
        }
        
    } else {
        
        ge c = *a;
        ge d = *b;
        
        r->infinity = 0;
        fe_sub(&num, &(d.y), &(c.y));
        fe_sub(&den, &(d.x), &(c.x));
        
        fe_inv(&lambda, &den);
        fe_mul(&lambda, &num, &lambda);
        
        fe_sqr(&(r->x), &lambda);
        fe_sub(&(r->x), &(r->x), &(c.x));
        fe_sub(&(r->x), &(r->x), &(d.x));
        
        fe_sub(&(r->y), &(c.x), &(r->x));
        fe_mul(&(r->y), &(r->y), &lambda);
        fe_sub(&(r->y), &(r->y), &(c.y));
        
    }
}

void ge_mul_int(ge *r, const ge *a, uint64_t n) {
    assert(ge_valid(a));
    
    ge b = *a;
    r->infinity = 1;

    while (n > 0) {
        if (n & 1) {
            ge_add(r, r, &b);
        }
        
        ge_double(&b, &b);
        n >>= 1;
    }
}

ge *ge_points(size_t *size) {
    size_t maxlen = 1;
    size_t len = 1;
    
    ge *points = malloc(maxlen * sizeof(ge));
    
    points[0] = inf;
    
    for (uint64_t m = 0; m < p; m++) {
        for (uint64_t n = 0; n < p; n++) {
            
            fe x = {m};
            fe y = {n};
            ge g = {x,y};
            
            if (ge_valid(&g)) {
                
                if (len == maxlen) {
                    maxlen *= 2;
                    points = realloc(points, maxlen * sizeof(ge));
                }
                
                ++len;
                points[len - 1] = g;
                
            }
            
            points = realloc(points, len * sizeof(ge));
            
        }
    }
    
    *size = len;
    return points;
}
