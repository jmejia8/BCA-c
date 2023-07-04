#include <math.h>
#include "objective.h"


/*
double evaluate(Problem problem, double *x)
{
    return problem.f(x, problem.n_vars);
}
*/

double evaluate(Problem problem, double *y)
{
    return problem.f(problem.x, y, problem.nx, problem.n_vars);
}

void evaluate_population(Problem problem, double* X, double* fitness, int N){
    int i, D = problem.n_vars;
    for (i = 0; i < N; ++i) {
        fitness[i] = evaluate(problem, &X[i*D]);
    }
}

