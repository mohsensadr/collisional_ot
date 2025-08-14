import torch
import pytest
from src import collision_pytorch as cp

def test_collOT_pytorch_shapes_and_types():
    # Fixed seed for deterministic behavior
    torch.manual_seed(0)

    # Create two small tensors (simulate "marginals")
    x = torch.randn(6, 2)  # 6 samples, 2D
    y = torch.randn(6, 2)

    # Run optimization
    x_out, y_out, dists_coll = cp.collOT_pytorch(
        x.clone(),
        y.clone(),
        avg_window=3,
        tol=1e-6,
        MinIter=2,
        MaxIter=10
    )

    # Shapes must match
    assert x_out.shape == x.shape
    assert y_out.shape == y.shape
    assert dists_coll.shape[0] >= 2  # At least initial + one iteration

    # All outputs must be torch tensors
    assert isinstance(x_out, torch.Tensor)
    assert isinstance(y_out, torch.Tensor)
    assert isinstance(dists_coll, torch.Tensor)

    # First distance matches initial mean squared error
    init_cost = torch.mean((x - y) ** 2)
    assert torch.isclose(dists_coll[0], init_cost, rtol=1e-6)

    # Final cost should be finite and non-negative
    assert torch.isfinite(dists_coll[-1])
    assert dists_coll[-1] >= 0

    # Should have at least one change in the distance array
    assert not torch.allclose(dists_coll[0], dists_coll[-1])

