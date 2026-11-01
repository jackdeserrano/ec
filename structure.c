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
