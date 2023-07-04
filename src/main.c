#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils.h"
#include "upper_level.h"
#include "bca.h"

double objective_function(double* x, int D)
{
    int i;
    double s = 0.0;
    for (i = 0; i < D; ++i)
        s += pow(x[i] - 1.0 / ((double) i+1), 2);
    return 100.0 + s;
}

double leader(double *x, double *y, int nx, int ny)
{
    int i;
    double s = 0.0;
    for (i = 0; i < nx; ++i)
        s += pow(x[i] - 1.0, 2);

    for (i = 0; i < ny; ++i)
        s += pow(y[i] - 1.0, 2);
    return s;
}

double follower(double *x, double *y, int nx, int ny)
{
    int i, n = nx < ny ? nx : ny;
    double s = 0.0;
    for (i = 0; i < n; ++i)
        s += pow(x[i] - y[i], 2.0);
    for (i = n; i < ny; ++i)
        s += pow(y[i] - 1.0, 2.0);
    return x[0] - 1.0 + s;
}

void test_eca()
{
    
    int D = 10;
    int K = 3;
    int N = 2*D*K;
    int max_iter = D*10000 / N;
    double eta_max = 2.0;
    int seed = time(NULL);
    double* low_bound = zeros(1,D);
    double* up_bound = zeros(1,D);
    int i;
    for (i = 0; i < D; ++i) {
        low_bound[i] = -100.0;
        up_bound[i] = 100.0;
    }

    double* X = zeros(N, D);
    double* fitness = (double *) malloc(sizeof(double)*N);

    // define problem
    Problem problem;
    problem.f = follower;
    problem.lb = low_bound;
    problem.ub = up_bound;
    problem.n_vars = D;
    problem.x = zeros(1, D);
    problem.nx = D;

    optimize(problem, X, fitness, N, K, eta_max, max_iter, seed);

    // get best solutions
    int i_best = minind(fitness, N);
    printf("x = ");
    print_vec(&X[i_best*D], D);
    printf("fx = %g\n", fitness[i_best]);
    printf("-----------\n");

    free(X);
    free(fitness);
    free(low_bound);
    free(up_bound);
}

void test_ul()
{
    
    int D = 5;
    int K = 7;
    int N = D*K;
    int Nll = D*K;
    double eta_max = 2.0;
    int seed = time(NULL);
    int max_iter = 100;
    int max_iter_ll = 5000;
    int verbose = 1;

    double* low_bound = zeros(1,D);
    double* up_bound = zeros(1,D);

    double* X = zeros(N, D);
    double* Y = zeros(N, D);
    double* fitness = zeros(N, 1);
    double* fitness_ll = zeros(N, 1);

    int i;
    for (i = 0; i < D; ++i) { low_bound[i] = -100.0; up_bound[i] = 100.0; }

    // UPPER LEVEL 
    Problem ul;
    ul.f = leader;
    ul.lb = low_bound;
    ul.ub = up_bound;
    ul.n_vars = D;
    // LOWER LEVEL 
    Problem ll;
    ll.f = follower;
    ll.lb = low_bound;
    ll.ub = up_bound;
    ll.n_vars = D;
    ll.nx = D;
    // bilevel problem
    BLProblem problem;
    problem.ul = ul;
    problem.ll = ll;


    upper_level(problem, X, Y, fitness, fitness_ll, N, Nll, K, eta_max, max_iter, max_iter_ll, seed, verbose, 1e-8);



    free(X);
    free(Y);
    free(fitness);
    free(fitness_ll);
    free(low_bound);
    free(up_bound);
}


void test_bca()
{
    
    int D = 2;
    int Dll = 3;
    int K = 7;
    int N = D*K;
    int Nll = D*K;
    double eta_max = 2.0;
    int seed = time(NULL);
    int max_iter = 100;
    int max_iter_ll = 5000;
    int verbose = 1;
    double tol = 1e-8;

    double* low_bound = zeros(1,D);
    double* up_bound = zeros(1,D);
    double* ll_low_bound = zeros(1, Dll);
    double* ll_up_bound =  zeros(1, Dll);

    double* X = zeros(N, D);
    double* Y = zeros(N, Dll);
    double* fitness = zeros(N, 1);
    double* fitness_ll = zeros(N, 1);

    int i;
    for (i = 0; i < D; ++i) { low_bound[i] = -100.0; up_bound[i] = 100.0; }
    for (i = 0; i < Dll; ++i) { ll_low_bound[i] = -100.0; ll_up_bound[i] = 100.0; }


    bca(leader, follower, 
                low_bound, up_bound,
                ll_low_bound, ll_up_bound,
                D, Dll, 
                X, Y, 
                fitness, fitness_ll, 
                N, Nll, K, eta_max, 
                max_iter, max_iter_ll, 
                seed, verbose, tol);

    printf("=====================\n");
    printf("    RESULT\n");
    printf("=====================\n");
    i = minind(fitness, N);
    printf("x = ");
    print_vec(&X[i*D], D);
    printf("y = ");
    print_vec(&Y[i*Dll], Dll);
    printf("F = %f\nf = %f\n", fitness[i], fitness_ll[i]);
    printf("=====================\n");


    free(X);
    free(Y);
    free(fitness);
    free(fitness_ll);
    free(low_bound);
    free(up_bound);
}

int main(int argc, char *argv[])
{
    test_bca();
    return 0;
}
