import numpy as np
from numpy.random import default_rng

def total_cost(x_m, ids_m):
  Nm_m = x_m.shape[0]
  Np_m = x_m.shape[1]
  y = 0
  for i in range(Nm_m):
    for j in range(i+1, Nm_m):
        y += np.sum( np.sum( (x_m[i, :, :]-x_m[j, :, :])**2 , axis=1) )
  return y/Np_m

def change_cost(x_m, ids_m, i1s_m, i2s_m, k_m, before=True):
  Nm_m = x_m.shape[0]
  s_m = np.ones(len(i1s_m))
  for j_m in range(Nm_m):
    if j_m != k_m:
        s_m += np.sum( np.sum( (x_m[j_m,:,:]-x_m[k_m,:,:])**2, axis=1) )

        if before is False:
          s_m -= np.sum( (x_m[j_m,i1s_m,:]-x_m[k_m,i1s_m,:])**2, axis=1)
          s_m -= np.sum( (x_m[j_m,i2s_m,:]-x_m[k_m,i2s_m,:])**2, axis=1)
          s_m += np.sum( (x_m[j_m,i1s_m,:]-x_m[k_m,i2s_m,:])**2, axis=1)
          s_m += np.sum( (x_m[j_m,i2s_m,:]-x_m[k_m,i1s_m,:])**2, axis=1)

  return s_m

def collOT_numpy(x, Total_Cost=total_cost, Change_Cost=change_cost, MinIter=20, MaxIter=100000, tol = 1e-3, avg_window=20):
    # x with size (Nmarginals, Np, dim)
    Nm = x.shape[0] # number of marginals
    Np = x.shape[1] # number of samples per marginal
    rng = default_rng()
    tries = Np
    sum_ = 1000000
    ids = np.array([[i for i in range(Np)] for j in range(Nm)]) # Nmarginals, Np

    dists_coll = [Total_Cost(x, ids)]

    for nt in range(1, MaxIter+1):
        for k in range(Nm):
            iss = rng.choice(Np, size=tries, replace=False)
            i1s = iss[:int(tries/2)]
            i2s = iss[int(tries/2):]

            s0 = Change_Cost(x, ids, i1s, i2s, k, before=True)
            s1 = Change_Cost(x, ids, i1s, i2s, k, before=False)

            # Determine which swaps to accept
            mask = s1 < s0

            # Perform the swaps for accepted cases
            accepted_i1s = i1s[mask]
            accepted_i2s = i2s[mask]

            x[k,accepted_i1s,:], x[k,accepted_i2s,:] = x[k,accepted_i2s,:].copy(), x[k,accepted_i1s,:].copy()

        dists_coll.append( Total_Cost(x, ids) )
            
        if nt>avg_window and nt > MinIter:
            sum_0 = np.sum(np.array(dists_coll)[-avg_window:])
            if abs(sum_ - sum_0)/sum_0 < tol:
              break
            sum_ = sum_0

    return x, dists_coll, nt



