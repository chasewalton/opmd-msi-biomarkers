import numpy as np

from opmd_msi.contrasts.effect_size import (
    cliffs_delta,
    hedges_g,
    paired_log2_fold_change,
)


def test_hedges_g_zero_for_identical_distributions():
    rng = np.random.default_rng(0)
    a = rng.normal(size=200)
    b = rng.normal(size=200)
    est = hedges_g(a, b)
    assert abs(est.estimate) < 0.2
    assert est.p_value > 0.05


def test_hedges_g_positive_for_shifted():
    rng = np.random.default_rng(0)
    a = rng.normal(loc=1.0, size=80)
    b = rng.normal(loc=0.0, size=80)
    est = hedges_g(a, b)
    assert est.estimate > 0.7
    assert est.p_value < 1e-5


def test_paired_log2fc():
    rng = np.random.default_rng(0)
    base = rng.normal(loc=5.0, size=40)
    case = base + 1.5 + rng.normal(scale=0.2, size=40)
    est = paired_log2_fold_change(case, base)
    assert abs(est.estimate - 1.5) < 0.1
    assert est.p_value < 1e-10


def test_cliffs_delta_sign():
    rng = np.random.default_rng(0)
    a = rng.normal(loc=1.0, size=80)
    b = rng.normal(loc=0.0, size=80)
    est = cliffs_delta(a, b)
    assert est.estimate > 0
