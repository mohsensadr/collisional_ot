#include "collision.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <string.h>


void get_max_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    // ru_maxrss is in kilobytes on most systems
    printf("mem: %ld\n", usage.ru_maxrss);
    //return usage.ru_maxrss;
}


long get_peak_memory_usage() {
    FILE *file = fopen("/proc/self/status", "r");
    if (file == NULL) {
        fprintf(stderr, "Could not open /proc/self/status\n");
        return -1;
    }

    char line[256];
    long peak_memory_kb = -1;

    while (fgets(line, sizeof(line), file)) {
        if (strncmp(line, "VmPeak:", 7) == 0) {
            sscanf(line, "VmPeak: %ld kB", &peak_memory_kb);
            break;
        }
    }

    fclose(file);
    return peak_memory_kb;
}

int main(int argc, char *argv[]) {
    size_t Nm = atoi(argv[1]); // = 2;   // Number of margins
    size_t Np = atoi(argv[2]); // = 1000;   // Number of samples per margin
    size_t dim = 2;  // Dimensions
    size_t MinIter = 1000;
    size_t MaxIter = 1000;
    double tol = 1e-6;
    size_t avg_window = 20;
    int ISvectorized = 1; // 1 for true, 0 for false
    // Allocate memory for a 3D array (Nm x Np x dim)
    double *dists = (double*)malloc(MaxIter * sizeof(double));

    printf("Nm: %d Np: %d\n", Nm, Np);

    long peak_memory_kb1 = get_peak_memory_usage();

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
            }
        }
    }

    // Call the function to test
    find_OT_collision_nd_nmargins(x, Nm, Np, dim, MinIter, MaxIter, tol, avg_window, dists, ISvectorized);

    //long max_memory_kb = get_max_memory_usage();
    //printf("Maximum memory usage: %ld B\n", max_memory_kb*1024);

    long peak_memory_kb2 = get_peak_memory_usage();
    printf("Peak memory usage: %ld KB\n", peak_memory_kb2 - peak_memory_kb1);

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
