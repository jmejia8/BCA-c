#include "utils.h"
#include "upper_level.h"

void bca(
        double (*F)(double*, double*, int, int),
        double (*f)(double*, double*, int, int),
        // problem
        double* ul_low_bound,
        double* ul_up_bound,
        double* ll_low_bound,
        double* ll_up_bound,
        int nx,
        int ny,
        // output
        double* X,
        double* Y,
        double* fitness,
        double* fitness_ll,
        // BCA parameters
        int N,
        int Nll,
        int K,
        double eta_max,
        // budget and options
        int max_iter,
        int max_iter_ll,
        int seed,
        int verbose,
        double tol
    )
{

    // UPPER LEVEL 
    Problem ul;
    ul.f = F;
    ul.lb = ul_low_bound;
    ul.ub = ul_up_bound;
    ul.n_vars = nx;
    // LOWER LEVEL 
    Problem ll;
    ll.f = f;
    ll.lb = ll_low_bound;
    ll.ub = ll_up_bound;
    ll.n_vars = ny;
    ll.nx = nx;
    // bilevel problem
    BLProblem problem;
    problem.ul = ul;
    problem.ll = ll;

    upper_level(problem, X, Y, fitness, fitness_ll, N, Nll, K, eta_max, max_iter, max_iter_ll, seed, verbose, tol);
}
