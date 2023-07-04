#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "rng.h"
#include "objective.h"

void get_subset(uint64_t *s, int *U, int N, int K){
    int i;
    for (i = 0; i < N*K; ++i)
        U[i] = i % N;

    randperm(s, U, N*K);
}


void mutation(double *h, double *x, double *c, double* u_best, double* u_worst,
              double eta, double* a, double* b, double p_bin, double p_exploit,
              double p, int D, uint64_t *s)
{
    int i, exploit = p < p_exploit;

    for (i = 0; i < D; ++i) {
        if (randm(s) < p_bin){
            h[i] = u_best[i];
            continue;
        }

        if (exploit)
            h[i] = x[i] + eta*(c[i] - u_worst[i]);
        else
            h[i] = x[i] + eta*(u_best[i] - c[i]);

        if ( h[i] < a[i] || b[i] < h[i]){
            double l = b[i] - a[i];
            h[i] = a[i] + l*randm(s);
        }

    }
}

void center(double *c, double *population, double *m, int *U, int* best,
            int* worst, int K, int D){
    int i, j;
    double M = 0.0;
    double m_max = fabs(m[U[0]]);

    for (i = 0; i < D; ++i) c[i] = 0.0;
    for (i = 1; i < K; ++i){
        double mm = fabs(m[U[i]]);
        if (m_max < mm) m_max = mm;
    }
    m_max *= 2.0;


    best[0]  = 0;
    worst[0] = 0;
    
    double m_best  = m_max - m[U[best[0]]];
    double m_worst = m_max - m[U[worst[0]]];

    for (i = 0; i < K; ++i) {
        double *x = &population[ U[i]*D ];
        double mass = m_max - m[U[i]];

        for (j = 0; j < D; ++j)
            c[j] += x[j] * mass;

        M += mass;

        if (mass > m_best){
            m_best = mass;
            best[0] = U[i];
        }

        if (mass < m_worst){
            m_worst = mass;
            worst[0] = U[i];
        }

    }

    for (i = 0; i < D; ++i) c[i] /= M;

}

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
    )
{
    const double p_bin = 0.02;
    const double p_exploit = 0.95;
    int best, worst, i;
    int *U = (int *) malloc(sizeof(int) * N*K);
    double *c  = zeros(1, D);

    get_subset(s, U, N, K);

    // use threads in this part
    for (i = 0; i < N; ++i) {
        // generate a center of mass
        center(c, population, fitness, &U[i*K], &best, &worst, K, D);
        // stepsize
        double eta = eta_max*randm(s);
        double *u_worst = &population[worst*D];
        double *u_best = &population[best*D];
        double *h = &X_new[i*D];
        // current solution
        double *x = &population[i*D];
        // variation operator
        mutation(h, x, c, u_best, u_worst, eta, a, b, p_bin, p_exploit, p, D, s);
    }

    free(U); free(c);
}


void replace(double* X, double* fitness, double* X_new, double* fitness_new, int N, int D){
    int i,j;
    for (i = 0; i < N; ++i) {
        if (fitness[i] > fitness_new[i]) {
            int w = maxind(fitness, N);
            for (j = 0; j < D; ++j)
                X[w*D + j] = X_new[i*D + j];
            fitness[w] = fitness_new[i];
        }
    }

}

int check_ftol(double* X, double* fitness, int N, int D, double tol){
    int imax = maxind(fitness, N);
    int imin = minind(fitness, N);
    return fitness[imax] - fitness[imin] < tol;
}

void initialize(
        Problem problem,
        double* X,
        double* fitness,
        int N,
        int K,
        double eta_max,
        uint64_t *s,
        int seed
    )
{

    set_seed(s, (uint64_t) seed);
    int D = problem.n_vars;

    double* a = problem.lb, *b = problem.ub;

    random_in_bounds(s, X, a, b, N, D);
    evaluate_population(problem, X, fitness, N);
}

void update_population(
        Problem problem,
        double* X_new,
        double* fitness_new,
        double* X,
        double* fitness,
        int N,
        int K,
        double eta_max,
        int iter,
        int max_iter,
        uint64_t *s
    )
{
    double p = (double) iter / (double) max_iter;
    int D = problem.n_vars;
    double* a = problem.lb, *b = problem.ub;
    // optimization steps
    reproduction(s, X_new, X, fitness, a, b, N, K, eta_max, D, p);
    evaluate_population(problem, X_new, fitness_new, N);
    replace(X, fitness, X_new, fitness_new, N, D);
}

void optimize(
        Problem problem,
        double* X,
        double* fitness,
        int N,
        int K,
        double eta_max,
        int max_iter,
        int seed
    )
{

    uint64_t s[4];
    const double tol = 1e-8;
    int D = problem.n_vars;

    double* X_new = zeros(N, D);
    double* fitness_new = zeros(N, 1);

    initialize(problem, X, fitness, N, K, eta_max, s, seed);

    int iter;
    for (iter = 0; iter < max_iter; ++iter) {
        update_population(problem, X_new, fitness_new, X, fitness, N, K, eta_max, iter, max_iter, s);
        // stop criterion
        if (check_ftol(X, fitness, N, D, tol)) break;
    }

    free(X_new);
    free(fitness_new);
}

