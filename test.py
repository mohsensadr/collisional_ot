import numpy as np
from matplotlib import pyplot as plt
from collision_wrapper import find_OT_collision

x = np.random.rand(2, 1000, 2)  # Example data

print( np.mean( np.sum( (x[0,:,:] - x[1,:,:])**2 , axis=1 ) ) )

x, dists_coll = find_OT_collision(x, 100, 10000)

print( np.mean( np.sum( (x[0,:,:] - x[1,:,:])**2 , axis=1 ) ) )

plt.plot(abs(dists_coll))
plt.yscale("log")
plt.show()
