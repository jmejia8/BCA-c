/*
typedef double (*ObjectiveFunction)(double *x, int nx);

typedef struct _str {
    ObjectiveFunction f;
    double* lb; // lower bounds
    double* ub; // upper bounds
    int n_vars;
} Problem;
*/

/*
typedef struct BLProblem_t {
    Problem ul;
    Problem ll;
} BLProblem;
*/


typedef double (*ObjectiveFunction)(double *x, double *y, int nx, int ny);
typedef struct _str {
    ObjectiveFunction f;
    double* lb; // lower bounds
    double* ub; // upper bounds
    int n_vars; // LL dim

    // for lower-level problem
    double* x; // current UL decision
    int nx; // size of the x
} Problem;


typedef struct BLProblem_t {
    Problem ul;
    Problem ll;
} BLProblem;

void evaluate_population(Problem, double*, double*, int);
