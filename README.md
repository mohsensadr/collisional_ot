# Collision-based dynamics for the optimal transport problem

In this repository, we present an implementation of collision-based dynamics for the optimal transport problem. This git repository has been used to produce results in the paper:

## C implementation
For faster runs on the CPU, we have prepared an implementation in C with a Python wrapper. To use it, first, compile the code by

```
cd src/
python3 setup.py build_ext --inplace
```

Then, in the Python code, simply import the collOT_c via
```
from collision_wrapper import collOT_c
```
and call the function via something like
```
X, dists_coll_xy, nt = collOT_c(X, MinIter=100, MaxIter=1000, tol = 1e-6, avg_window=20, Track=1)
```

For examples of using this implementation, see the Jupyter Notebooks in ```examples/``` directory.
