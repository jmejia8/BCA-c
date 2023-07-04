#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "rng.h"


void set_seed(uint64_t *s, uint64_t _s) {
	s[0] = _s + 111;
	s[1] = _s + 121;
	s[2] = _s + 113;
	s[3] = _s + 411;
	jump(s);
}

double randm(uint64_t *s) {
	return (double)next(s) / (double)(UINT64_MAX);
}

uint64_t randuint(uint64_t *s, uint64_t a, uint64_t b) {
	return a + next(s) % (b - a + 1);
}

double* randvec(uint64_t *s, int len) {
	double* v = (double*) malloc(sizeof(double)*len);
	int i;
	for (i = 0; i < len; ++i)
		v[i] = randm(s);
	return v;
}

double* randmat(uint64_t *s, int rows, int cols) {
	int len = rows * cols;
	double* v = (double*) malloc(sizeof(double)*len);
	int i;
	for (i = 0; i < len; ++i)
		v[i] = randm(s);
	return v;
}


void randperm(uint64_t *s, int *r, int n){
    int i;
    for (i = n-1; i >= 0; --i){
        //generate a random number [0, n-1]
        int j = randuint(s, 0, i) ;//% (i+1);

        //swap the last element with element at random index
        int temp = r[i];
        r[i] = r[j];
        r[j] = temp;
    }
}

void random_in_bounds(uint64_t* s, double *X, double* a, double* b, int N, int D){
    int i,j;

	for (i = 0; i < N; ++i){
		for (j = 0; j < D; ++j) {
			double l = b[j] - a[j];
			X[i*D + j] = a[j] + l*randm(s);
		}
	}

}
