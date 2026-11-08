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
