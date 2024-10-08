#include "collision.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    size_t Nm = 2;   // Number of margins
    size_t Np = 100;   // Number of samples per margin
    size_t dim = 2;  // Dimensions
    size_t MinIter = 100;
    size_t MaxIter = 1000;
    double tol = 1e-3;
    size_t avg_window = 20;
    int ISvectorized = 1; // 1 for true, 0 for false

    // Allocate memory for a 3D array (Nm x Np x dim)
    double ***x = (double***)malloc(Nm * sizeof(double**));
    for (size_t i = 0; i < Nm; ++i) {
        x[i] = (double**)malloc(Np * sizeof(double*));
        for (size_t j = 0; j < Np; ++j) {
            x[i][j] = (double*)malloc(dim * sizeof(double));
        }
    }

    // Initialize x with some values for testing
    for (size_t i = 0; i < Nm; ++i) {
        for (size_t j = 0; j < Np; ++j) {
            for (size_t k = 0; k < dim; ++k) {
                x[i][j][k] = rand() / ( 1.*RAND_MAX ); // uniform samples in [0,1]
                //printf("x[%d][%d][%d]:%lf\n",i,j,k,x[i][j][k]);
            }
        }
    }

    // Call the function to test
    find_OT_collision_nd_nmargins(x, Nm, Np, dim, MinIter, MaxIter, tol, avg_window, ISvectorized);


    // Free allocated memory
    for (size_t i = 0; i < Nm; ++i) {
        for (size_t j = 0; j < Np; ++j) {
            free(x[i][j]); // Free each 1D array
        }
        free(x[i]); // Free each 2D array
    }
    free(x); // Free the top-level pointer

    return 0; // Successful execution
}
