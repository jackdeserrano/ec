#include <stdint.h>
#include <stdio.h>
#include <assert.h>

uint64_t p;

typedef struct {
    uint64_t n;
} fe;

void fe_add(fe *r, const fe *a, const fe *b) {
    r->n = ((__uint128_t)a->n + b->n) % p;
}

void fe_add_int(fe *r, const fe *a, const uint64_t b) {
    r->n = ((__uint128_t)a->n + b) % p;
}

void fe_sub(fe *r, const fe *a, const fe *b) {
    r->n = ((__uint128_t)p + a->n - b->n) % p;
}

void fe_neg(fe *r, const fe *a) {
    r->n = (p - a->n) % p;
}

void fe_mul(fe *r, const fe *a, const fe *b) {
    r->n = ((__uint128_t)a->n * b->n) % p;
}

void fe_mul_int(fe *r, const fe *a, const uint64_t b) {
    r->n = ((__uint128_t)a->n * b) % p;
}

void fe_sqr(fe *r, const fe *a) {
    fe_mul(r, a, a);
}

void fe_exp(fe *r, const fe *a, uint64_t e) {
    fe x = *a;
    r->n = 1;

    while (e > 0) {
        if (e & 1) {
	    fe_mul(r,r,&x);
	}

	fe_sqr(&x,&x);
	e >>= 1;
    }
}

void fe_inv(fe *r, const fe *a) {
    assert(a->n);
    fe_exp(r,a,p-2);
}

int fe_equal(const fe *a, const fe *b) {
    return a->n == b->n;
}
