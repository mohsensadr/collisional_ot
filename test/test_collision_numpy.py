import numpy as np
import pytest
from src import collision_numpy as cn

def test_total_cost_known_case():
    # Simple deterministic case: 2 marginals, 2 points, 1D
    x_m = np.array([
        [[0.0], [1.0]],  # marginal 0
        [[1.0], [2.0]],  # marginal 1
    ])
    ids_m = np.array([[0, 1], [0, 1]])
    # Manually compute: pairwise distance squared for each sample index
    # For marginal 0 and 1:
    # Point 0: (0-1)^2 = 1
    # Point 1: (1-2)^2 = 1
    # Sum = 2, divide by Np=2 => 1.0
    expected_cost = 1.0
    assert np.isclose(cn.total_cost(x_m, ids_m), expected_cost)


def test_change_cost_behavior():
    # Create 3 marginals, 4 points each, 2D
    rng = np.random.default_rng(0)
    x_m = rng.normal(size=(3, 4, 2))
    ids_m = np.array([np.arange(4)] * 3)
    i1s = np.array([0, 1])
    i2s = np.array([2, 3])
    k_m = 0

    before_cost = cn.change_cost(x_m, ids_m, i1s, i2s, k_m, before=True)
    after_cost = cn.change_cost(x_m, ids_m, i1s, i2s, k_m, before=False)

    # These arrays should have same shape and differ numerically
    assert before_cost.shape == after_cost.shape
    assert not np.allclose(before_cost, after_cost)


def test_collOT_numpy_reduces_cost():
    # Random but seeded for reproducibility
    rng = np.random.default_rng(1)
    x = rng.normal(size=(3, 6, 2))  # 3 marginals, 6 samples each, dim=2

    # Get initial cost
    ids = np.array([np.arange(6)] * 3)
    init_cost = cn.total_cost(x.copy(), ids)

    # Run optimization
    x_out, dists_coll, nt = cn.collOT_numpy(x.copy(), MinIter=5, MaxIter=50, tol=1e-6)

    # The cost sequence should be non-empty and start near init_cost
    assert len(dists_coll) >= 2
    assert np.isclose(dists_coll[0], init_cost)

    # Final cost should be <= initial cost
    assert dists_coll[-1] <= init_cost + 1e-8

    # Should have iterated at least once
    assert nt >= 1

