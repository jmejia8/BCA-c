#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

void print_vec(double* v, int len) {
	int i;
	for (i = 0; i < len; ++i) {
		printf("%f ", v[i]);
	}
	printf("\n");
}

void print_mat(double* v, int rows, int cols) {
	int i;
	for (i = 0; i < rows; ++i) {
		print_vec(&v[i*cols], cols);
	}
}

double* zeros(int rows, int cols) {
	int len = rows * cols;
	double* v = (double*) malloc(sizeof(double)*len);
	int i;
	for (i = 0; i < len; ++i) v[i] = 0.0;
	return v;
}

double* ones(int rows, int cols) {
	int len = rows * cols;
	double* v = (double*) malloc(sizeof(double)*len);
	int i;
	for (i = 0; i < len; ++i) v[i] = 1.0;
	return v;
}

int maxind(double *f, int N){
    int i;
    int max = 0;
    for (i = 1; i < N; ++i) {
        if (f[i] > f[max]) {
            max = i;
        }
    }

    return max;
}

int minind(double *f, int N){
    int i;
    int min = 0;
    for (i = 1; i < N; ++i) {
        if (f[i] < f[min]) {
            min = i;
        }
    }

    return min;
}
