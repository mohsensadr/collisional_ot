// collision.c

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

double total_cost(double*** x, int** ids, int Nm, int Np, int dim) {
    double cost = 0;
    for (int i = 0; i < Nm; i++) {
        for (int j = i + 1; j < Nm; j++) {
            for (int n = 0; n < Np; n++) {
                for (int d = 0; d < dim; d++) {
                    double diff = x[i][n][d] - x[j][n][d];
                    cost += diff * diff;
                }
            }
        }
    }
    return cost / Np;
}

double change_cost(double*** x, int i, int* i1s, int* i2s, int k, int Nm, int Np, int dim) {
    double s_m = 0.0;
    for (int j = 0; j < Nm; j++) {
        if(j!=k){
                    for (int d = 0; d < dim; d++) {
                        double diff1 = x[j][i1s[i]][d] - x[k][i1s[i]][d];
                        double diff2 = x[j][i2s[i]][d] - x[k][i2s[i]][d];
                        s_m -= ( diff1 * diff1 + diff2 * diff2 );

                        double new_diff1 = x[j][i1s[i]][d] - x[k][i2s[i]][d];
                        double new_diff2 = x[j][i2s[i]][d] - x[k][i1s[i]][d];
                        s_m += new_diff1 * new_diff1 + new_diff2 * new_diff2;
                    }
        }
    }
    return s_m;
}

/* Arrange the N elements of ARRAY in random order.
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. */
void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

// The main function for OT collision
int find_OT_collision_nd_nmargins(double*** x, int Nm, int Np, int dim, int MinIter, int MaxIter, double tol, int avg_window, double *dists_coll) {
    int tries = Np;
    double sum_ = 1e6, s0, s1;
    //double dists_coll[MaxIter + 1];
    int* iss = (int*) malloc(tries * sizeof(int));
    int** ids = (int**) malloc(Nm * sizeof(int*));
    int* i1s;
    int* i2s;
    int nt;

    for (int i = 0; i < Nm; i++) {
        ids[i] = (int*) malloc(Np * sizeof(int));
        for (int j = 0; j < Np; j++) {
            ids[i][j] = j;
        }
    }

    for (int i = 0; i < tries; i++) iss[i] = i;

    s0 = total_cost(x, ids, Nm, Np, dim)*Np;
    dists_coll[0] = s0/Np;
    for (nt = 1; nt <= MaxIter; nt++) {
        for (int k = 0; k < Nm; k++) {
            shuffle(iss, tries);

            i1s = iss;
            i2s = iss + tries / 2;

            for (int i = 0; i < tries / 2; i++) {
                    s1 = change_cost(x, i, i1s, i2s, k, Nm, Np, dim);
                    //printf("i1s:%d i2s:%d, s1: %lf\n",i1s[i],i2s[i], s1);
                    if (s1 < 0) {
                        for (int d = 0; d < dim; d++) {
                            double temp = x[k][i1s[i]][d];
                            x[k][i1s[i]][d] = x[k][i2s[i]][d];
                            x[k][i2s[i]][d] = temp;
                        }
                        s0 += s1;
                    }
            }
        }

        dists_coll[nt] = s0/Np;
        //printf("\n\n");
        //printf("%e\n", dists_coll[nt]);
        if (nt > avg_window && nt > MinIter) {
            double sum_0 = 0;
            for (int i = nt - avg_window; i < nt; i++) {
                sum_0 += dists_coll[i];
            }
            if (fabs(sum_ - sum_0) / sum_0 < tol) {
                break;
                MaxIter = nt;
            }
            sum_ = sum_0;
        }
    }

    for (int i = 0; i < Nm; i++) free(ids[i]);
    free(ids);
    return nt;
}
