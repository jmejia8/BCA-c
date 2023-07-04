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
    );
