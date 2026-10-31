#include <math.h>
#include <primesieve.h>
#include "group.c"

int main(void) {
    uint64_t *pp = primesieve_generate_n_primes(100, 0, UINT64_PRIMES);
    
    for (int i = 0; i < 100; i++) {
        p = pp[i];
        
        if (p == 2 || p == 3 || p == 7) {
            continue; // bad reduction
        }
        
        size_t len;
        
        ge *l = ge_points(&len);
        
        printf("%zu,",len); // card(E(F_{p}))
        free(l);
    }
    
    printf("\n");
    
    for (int i = 0; i < 100; i++) {
        p = pp[i];
        
        if (p == 2 || p == 3 || p == 7) {
            continue;
        }
        
        if ((p + 1) % 3 == 0) {
            printf("%llu,", p + 1);
        } else {
            
            size_t len;
            
            ge *l = ge_points(&len);
            
            printf("%zu,",len);
            free(l);
            
        }
    }
    
    printf("\n");
    primesieve_free(pp);
    
    return 0;
}
