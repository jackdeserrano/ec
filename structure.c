#include "group.c"

int has_good_reduction(void) {
    fe discriminant;

    fe a = {A};
    fe b = {B};

    fe_exp(&a, &a, 3);
    fe_mul_int(&a, &a, 4);

    fe_sqr(&b, &b);
    fe_mul_int(&b, &b, 27);

    fe_add(&discriminant, &a, &b);
    fe_mul_int(&discriminant, &discriminant, 16);

    return discriminant.n != 0;
}

uint64_t abs_trace(void) {
    size_t size;
    ge_points(&size);
    return p + 1 > size ? p + 1 - size : size - p - 1;
}

int is_supersingular(void) {
    return abs_trace() % p == 0;
}

uint64_t ge_order(const ge *a) {
    uint64_t order = 1;
    ge b = *a;
    while (!(b.infinity)) {
        ge_add(&b, &b, a);
        ++order;
    }
    return order;
}

// E(F_p) is isomorphic to Z/n_1Z \times Z/n_2Z
void group_structure(uint64_t *n_1, uint64_t *n_2) {
    size_t size;
    ge *points = ge_points(&size);
    uint64_t max_order = 1;
    for (int i = 0; i < size; ++i) {
        ge point = points[i];
        uint64_t order = ge_order(&point);
        max_order = order > max_order ? order : max_order;
    }
    *n_1 = size / max_order;
    *n_2 = max_order;
    free(points);
}
