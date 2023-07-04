#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "utils.h"
#include "eca.h"
#include "rng.h"


void _copy(double* a, double* b, int n)
{
    int i;
    for (i = 0; i < n; ++i) { a[i] = b[i]; }
}

void lower_level(
        Problem problem,
        double* X,
        double* Y,
        double* fitness,
        int N,
        int K,
        double eta_max,
        int max_iter,
        uint64_t *s
    )
{
    int i, j, ny = problem.n_vars;
    int nx = problem.nx;

    double* Y_new = zeros(N, ny);
    double* fitness_new = ones(N, 1);

    // TODO use threads
    for (i = 0; i < N; ++i) {
        // parametrize lower-level
        problem.x = &X[i*nx];
        // lower-level optimization
        optimize(problem, Y_new, fitness_new, N, K, eta_max, max_iter, 1);

        j = minind(fitness_new, N);
        _copy(&Y[i*ny], &Y_new[j*ny], ny);
        fitness[i] = fitness_new[j];
    }

    free(Y_new);
    free(fitness_new);
}

void evaluate_ul(BLProblem problem, double* X, double* Y, double* fitness, int N)
{
    int i, nx = problem.ul.n_vars, ny = problem.ll.n_vars;
    for (i = 0; i < N; ++i) {
        fitness[i] = problem.ul.f(&X[i*nx], &Y[i*ny], nx, ny);
    }

}

void replace_ul(
        double* X, double* Y,
        double* X_new, double* Y_new,
        double* fitness, double* fitness_new,
        double* fitness_ll, double* fitness_new_ll,
        int N, int nx, int ny
        )
{
    int i;
    for (i = 0; i < N; ++i) {
        if (fitness[i] > fitness_new[i]) {
            int w = maxind(fitness, N);
            _copy(&X[w*nx], &X_new[i*nx], nx);
            _copy(&Y[w*ny], &Y_new[i*ny], ny);
            fitness[w] = fitness_new[i];
            fitness_ll[w] = fitness_new_ll[i];
        }
    }

}

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
    )
{

    uint64_t s[4];
    int nx = problem.ul.n_vars;
    int ny = problem.ll.n_vars;
    double* a = problem.ul.lb, *b = problem.ul.ub;

    double* X_new = zeros(N, nx);
    double* Y_new = zeros(Nll, ny);
    double* fitness_new = zeros(1, N);
    double* fitness_new_ll = ones(1, N);

    // initialization
    random_in_bounds(s, X, a, b, N, nx);
    lower_level(problem.ll, X, Y, fitness_ll, N, K, eta_max, max_iter_ll, s);
    evaluate_ul(problem, X, Y, fitness, N);

    int iter;
    for (iter = 0; iter < max_iter; ++iter) {
        double p = (double) iter / (double) max_iter;
        reproduction(s, X_new, X, fitness, a, b, N, K, eta_max, nx, p);
        lower_level(problem.ll, X_new, Y_new, fitness_new_ll, Nll, K, eta_max, max_iter_ll, s);
        evaluate_ul(problem, X_new, Y_new, fitness_new, N);
        replace_ul(X, Y, X_new, Y_new, fitness, fitness_new, fitness_ll, fitness_new_ll, N, nx, ny);
        // stop criterion
        if (check_ftol(X, fitness, N, nx, tol)) break;

        if (verbose) {
            int jj = minind(fitness, N);
            printf("iter: %d  F = %f  f = %f\n", iter, fitness[jj], fitness_ll[jj]);
        }
    }

    free(X_new);
    free(Y_new);
    free(fitness_new);
    free(fitness_new_ll);
}


