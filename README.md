![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

# Collisional Optimal Transport

In this repository, we present an implementation of collision-based dynamics for the optimal transport problem. This git repository has been used to produce results in the following paper:

Sadr, Mohsen, and Hossein Gorji. "Collision-based Dynamics for Multi-Marginal Optimal Transport." arXiv preprint at [arXiv:2412.16385 (2024)](https://doi.org/10.48550/arXiv.2412.16385)

##  Implementation in PyTorch
We provide a PyTorch implementation of the collision-based dynamics to solve the optimal transport problem to be used as a loss function in training statistical models. First, import the function via
```
from collision import collOT_pytorch
```
and then call it simply by
```
x, y, log_loss = collOT_pytorch(x, y, MinIter=100, MaxIter=1000, tol = 1e-6, avg_window=20, Track=1)
```

## Implementation in C
For faster runs on the CPU, we have prepared an implementation in C with a Python wrapper. To use it, first, compile the code by

```
cd src/
python3 setup.py build_ext --inplace
```

Then, in the Python code, import the ```collOT_c``` via
```
from collision_wrapper import collOT_c
```
and call the function via something like
```
X, log_loss, nsteps = collOT_c(X, MinIter=100, MaxIter=1000, tol = 1e-6, avg_window=20, Track=1)
```

For examples of how this implementation can be used, see the Jupyter Notebooks in ```examples/``` directory.
