![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

# Collisional Optimal Transport

In this repository, we present an implementation of collision-based dynamics for the optimal transport problem. This git repository has been used to produce results in this paper:

## PyTorch implementation
We provide a PyTorch implementation of the collision-based dynamics to solve optimal transport, with the application as a loss function. First import it
```
from collision import collOT_pytorch
```
and then call it simply by
```
x, y, log_loss = collOT_pytorch(x, y, MinIter=100, MaxIter=1000, tol = 1e-6, avg_window=20, Track=1)
```

## C implementation
For faster runs on the CPU, we have prepared an implementation in C with a Python wrapper. To use it, first, compile the code by

```
cd src/
python3 setup.py build_ext --inplace
```

Then, in the Python code, simply import the ```collOT_c``` via
```
from collision_wrapper import collOT_c
```
and call the function via something like
```
X, log_loss, nsteps = collOT_c(X, MinIter=100, MaxIter=1000, tol = 1e-6, avg_window=20, Track=1)
```

For examples of using this implementation, see the Jupyter Notebooks in ```examples/``` directory.
