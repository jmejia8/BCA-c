#include <stdint.h>
#include "objective.h"

void optimize(Problem, double*, double*, int, int, double, int, int);

void reproduction(
        uint64_t *s,
        double* X_new,
        double* population,
        double* fitness,
        double* a,
        double* b,
        int N,
        int K,
        int eta_max,
        int D,
        double p
    );

int check_ftol(double* X, double* fitness, int N, int D, double tol);

