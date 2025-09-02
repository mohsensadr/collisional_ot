// collision.c

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

double total_cost(double*** x, int Nm, int Np, int dim) {
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

int min(int a, int b){
  if(a<b){
     return a;
  }
  else{
     return b;
  }
}

int max(int a, int b){
  if(a<b){
     return b;
  }
  else{
     return a;
  }
}

double change_cost(double*** x, int i1, int i2, int k, int Nm, int Np, int dim, int bi) {
    double s_m = 0.0;
    double diff1, diff2, new_diff1, new_diff2;
    if(bi == 1){
        for (int j = max(0,k-1); j < min(Nm,k+2); j++) {
            if(j!=k){
                    for (int d = 0; d < dim; d++) {
                        diff1 = x[j][i1][d] - x[k][i1][d];
                        diff2 = x[j][i2][d] - x[k][i2][d];
                        s_m -= ( diff1 * diff1 + diff2 * diff2 );

                        new_diff1 = x[j][i1][d] - x[k][i2][d];
                        new_diff2 = x[j][i2][d] - x[k][i1][d];
                        s_m += new_diff1 * new_diff1 + new_diff2 * new_diff2;
                    }
            }
        }
    }
    else{
        for (int j = 0; j < Nm; j++) {
            if(j!=k){
                    for (int d = 0; d < dim; d++) {
                        diff1 = x[j][i1][d] - x[k][i1][d];
                        diff2 = x[j][i2][d] - x[k][i2][d];
                        s_m -= ( diff1 * diff1 + diff2 * diff2 );

                        new_diff1 = x[j][i1][d] - x[k][i2][d];
                        new_diff2 = x[j][i2][d] - x[k][i1][d];
                        s_m += new_diff1 * new_diff1 + new_diff2 * new_diff2;
                    }
            }
        }
    }
    return s_m;
}

void shuffle(int *array, size_t n) {
    size_t i, j, t;
    if (n > 1) {
        for (i = 0; i < n - 1; i++) {
            j = i + 1 + rand() / (RAND_MAX / (n - i - 1) + 1 );
            t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}


int find_OT_collision_nd_nmargins(double*** x, int Nm, int Np, int dim, int MinIter, int MaxIter, double tol, int avg_window, double *dists_coll, int Track, int bi, int resh_freq) {
    int Ncoll = floor(Np / 2);
    double sum_ = 1e6, s0, s1, temp, sum_0;
    int* iss = (int*) malloc(2*Ncoll * sizeof(int));
    int k, i, d, nt;
    int i1, i2;
    int r, q, R;

    for (i = 0; i < 2*Ncoll; i++) iss[i] = i;

    if( Track == 1){
        s0 = total_cost(x, Nm, Np, dim)*Np;
    }
    else{
        s0 = 0.0;
    }
    dists_coll[0] = s0/Np;
    R = 1;
    for (nt = 1; nt <= MaxIter; nt++) {
        for (k = 0; k < Nm; k++) {
            for(q=0; q<floor(Np/(2*R)); q++){
                 for(r=0; r<R; r++){
                     i1 = iss[2*R*q+r];
                     i2 = iss[(2*R*q+R+r)%(Np+1)];
                     s1 = change_cost(x, i1, i2, k, Nm, Np, dim, bi);
                     if (s1 < 0) {
                        for (d = 0; d < dim; d++) {
                            temp = x[k][i1][d];
                            x[k][i1][d] = x[k][i2][d];
                            x[k][i2][d] = temp;
                        }
                        s0 += s1;
                     }
                 }
            }
        }
        R++;
        if( (R==floor(Np/2)) || nt%resh_freq==0 ){
            R = 1;
	    shuffle(iss, 2*Ncoll);
        }

        dists_coll[nt] = s0/Np;
        if (nt > avg_window && nt > MinIter) {
            sum_0 = 0;
            for (i = nt - avg_window; i < nt; i++) {
                sum_0 += dists_coll[i];
            }
            if (fabs(sum_ - sum_0) / sum_0 < tol) {
                break;
                MaxIter = nt;
            }
            sum_ = sum_0;
        }
    }
    free(iss);
    return nt;
}


// The main function for OT collision
int find_OT_ISA_c(double*** x, int Nm, int Np, int dim, int MinIter, int MaxIter, double tol, int avg_window, double *dists_coll, int Track, int bi) {
    double sum_ = 1e6, s0, s1;
    int nt;

    if( Track == 1){
        s0 = total_cost(x, Nm, Np, dim)*Np;
    }
    else{
        s0 = 0.0;
    }
    dists_coll[0] = s0/Np;
    for (nt = 1; nt <= MaxIter; nt++) {
        for (int k = 1; k < Nm; k++) {
	    for (int i = 0; i < Np; i++) {
		for (int j = i; j < Np; j++) {
                        s1 = change_cost(x, i, j, k, Nm, Np, dim, bi);
                        if (s1 < 0) {
                            for (int d = 0; d < dim; d++) {
                                double temp = x[k][i][d];
                            	x[k][i][d] = x[k][j][d];
                            	x[k][j][d] = temp;
                            }
                            s0 += s1;
                        }
		}
            }
        }

        dists_coll[nt] = s0/Np;
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
    return nt;
}
