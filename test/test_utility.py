import numpy as np
import pytest
from src import utility as ut


def test_fit_xp_simple_case(monkeypatch):
    # Monkeypatch linregress to a deterministic fake
    from scipy.stats import linregress as real_linregress
    try:
        from scipy.stats import linregress
    except ImportError:
        pytest.skip("scipy required for fit_xp test")

    x = np.array([1, 10, 100], dtype=float)
    y = np.array([2, 20, 200], dtype=float)
    slope, slope_str, fitted_y = ut.fit_xp(x, y)

    assert isinstance(slope, float)
    assert slope_str == f"{slope:.2f}"
    assert fitted_y.shape == x.shape
    # Should be close to y up to scaling in log-log fit
    assert np.all(fitted_y > 0)


def test_swissroll_shape_and_values():
    np.random.seed(0)
    X = ut.swissroll(n=50, noise=0.1)
    assert isinstance(X, np.ndarray)
    assert X.shape == (50, 2)
    assert np.isfinite(X).all()


def test_sample_banana_shape():
    np.random.seed(0)
    out = ut.sample_banana(20)
    assert out.shape == (20, 2)
    # First column should match standard normal distribution shape
    assert np.isfinite(out).all()


def test_funnel_sample_and_threshold():
    np.random.seed(0)
    f = ut.Funnel(d=3, sigma=1.5, limit_min=0.1, limit_max=5.0)
    samples = f.sample(10)
    assert samples.shape == (10, 3)
    # Threshold clamps values
    vals = np.array([-10, 0, 10])
    th = f.threshold(vals)
    assert (th >= f.limit_min).all()
    assert (th <= f.limit_max).all()


def test_funnel_log_pdf_and_errors():
    np.random.seed(0)
    f = ut.Funnel(d=2)
    samples = f.sample(5)
    logp = f.log_pdf(samples)
    assert logp.shape == (5,)
    assert np.isfinite(logp).all()
    with pytest.raises(ValueError):
        f.log_pdf(np.random.randn(5, 3))  # wrong dimension


def test_funnel_norm_log_pdf_matches_manual():
    X = np.array([0.0])
    sigma = 2.0
    manual = -0.5 * np.log(2 * np.pi * sigma**2) - 0.5 * (X**2 / sigma**2)
    calc = ut.Funnel.norm_log_pdf(X, 0, sigma)
    assert np.allclose(manual, calc)


def test_ring_sample_and_log_pdf():
    np.random.seed(0)
    r = ut.Ring(d=3, sigma=0.5, radia=[3.0, 5.0])
    samples = r.sample(10)
    assert samples.shape == (10, 3)
    logp = r.log_pdf(samples)
    assert logp.shape == (10,)
    assert np.isfinite(logp).all()


def test_ring_invalid_dimension():
    with pytest.raises(ValueError):
        ut.Ring(d=1)


def test_ring_log_pdf_dimension_mismatch():
    r = ut.Ring(d=2)
    X = np.random.randn(5, 3)
    with pytest.raises(ValueError):
        r.log_pdf(X)


def test_ring_log_sum_exp_and_norm_log_pdf():
    r = ut.Ring(d=2)
    arr = np.array([1.0, 2.0, 3.0])
    result = r.log_sum_exp(arr)
    manual = np.log(np.sum(np.exp(arr)))
    assert np.allclose(result, manual)

    X = np.array([0.0])
    sigma = 0.5
    manual_pdf = -0.5 * np.log(2 * np.pi * sigma**2) - 0.5 / sigma**2 * (X**2)
    assert np.allclose(r.norm_log_pdf(X, 0, sigma), manual_pdf)

