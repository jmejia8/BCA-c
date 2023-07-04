#include "xoshiro256plusplus.h"

void set_seed(uint64_t*, uint64_t);
double randm(uint64_t*);
uint64_t randuint(uint64_t*, uint64_t, uint64_t);
double* randmat(uint64_t*, int, int);
double* randvec(uint64_t*, int);
void randperm(uint64_t *, int*, int);
void random_in_bounds(uint64_t*, double *, double* , double*, int, int);
