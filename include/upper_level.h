#include "objective.h"

void upper_level(
        BLProblem problem,
        double* X,
        double* Y,
        double* fitness,
        double* fitness_ll,
        int N,
        int Nll,
        int K,
        double eta_max,
        int max_iter,
        int max_iter_ll,
        int seed,
        int verbose,
        double tol
    );
