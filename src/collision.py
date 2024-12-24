import torch

def collOT_pytorch(x, y, avg_window=10, tol=1e-6, MinIter=10, MaxIter=100):
      tries = int(len(x))
      sum_ = 1000000
      dists_coll = torch.zeros(MaxIter+1)
      dists_coll[0] = torch.mean((x - y)**2)
      for nt in range(MaxIter):
        iss = torch.randperm(len(x))
        i1s = iss[:int(tries/2)]
        i2s = iss[int(tries/2):]

        # Calculate the initial distances
        s0 = torch.sum( (x[i1s] - y[i1s])**2, axis=-1) + torch.sum( (x[i2s] - y[i2s])**2, axis=-1)

        # Calculate the distances after the swap
        s1 = torch.sum( (x[i1s] - y[i2s])**2, axis=-1) + torch.sum( (x[i2s] - y[i1s])**2, axis=-1)

        # Determine which swaps to accept
        mask = s1 < s0

        # Perform the swaps for accepted cases
        accepted_i1s = i1s[mask]
        accepted_i2s = i2s[mask]

        y[accepted_i1s], y[accepted_i2s] = y[accepted_i2s], y[accepted_i1s]

        dists_coll[nt+1] = torch.mean((x - y)**2)
        if nt>avg_window and nt > MinIter:
          sum_0 = torch.sum(dists_coll[-avg_window:])
          if abs(sum_ - sum_0)/sum_0 < tol:
            break
          sum_ = sum_0
      return x, y, dists_coll
