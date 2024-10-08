# cython: language_level=3

import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from cython cimport boundscheck, wraparound

cdef extern from "collision.h":
    int find_OT_collision_nd_nmargins(double*** x, int Nm, int Np, int dim, int MinIter, int MaxIter, double tol, int avg_window, double *dists_coll)

@boundscheck(False)
@wraparound(False)
def find_OT_collision(cnp.ndarray[cnp.float64_t, ndim=3, mode="c"] x, int MinIter=20, int MaxIter=100000, double tol=1e-3, int avg_window=20, int ISvectorized=1):
    cdef int Nm = x.shape[0]
    cdef int Np = x.shape[1]
    cdef int dim = x.shape[2]
    cdef int nt = 0

    # Allocate memory for the pointers to the arrays
    cdef double*** x_c = <double***>malloc(Nm * sizeof(double**))
    for i in range(Nm):
        x_c[i] = <double**>malloc(Np * sizeof(double*))
        for j in range(Np):
            x_c[i][j] = &x[i, j, 0]

    cdef double* dists_coll = <double*>malloc((MaxIter+1)*sizeof(double))

    nt = find_OT_collision_nd_nmargins(x_c, Nm, Np, dim, MinIter, MaxIter, tol, avg_window, dists_coll)
    MaxIter = nt-1

    dists_coll_np = np.empty(MaxIter+1, dtype=np.float64)
    for i in range(MaxIter+1):
         dists_coll_np[i] = dists_coll[i]

    # Free allocated memory
    for i in range(Nm):
        free(x_c[i])
    free(x_c)
    free(dists_coll)
    return x, dists_coll_np
