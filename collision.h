#ifndef COLLISION_H
#define COLLISION_H

#include <stdint.h>
#include <stddef.h>  // for size_t
#include <stdio.h>

// Function prototypes
int find_OT_collision_nd_nmargins(double*** x, size_t Nm, size_t Np, size_t dim,
                                  size_t MinIter, size_t MaxIter, double tol,
                                  size_t avg_window, double *dists_coll);

double Total_Cost(double*** x, size_t Nm, size_t Np, int* ids);

double Change_Cost(double*** x, int* ids, size_t Nm, size_t Np, 
                   size_t* i1s, size_t* i2s, size_t k, 
                   int before, double* result);

#endif // COLLISION_H
