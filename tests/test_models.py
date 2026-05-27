import numpy as np
import pandas as pd

from opmd_msi.models.logistic import fit_l1_logistic


def test_l1_logistic_recovers_strong_signal():
    rng = np.random.default_rng(0)
    n = 200
    X = rng.normal(size=(n, 5))
    # Only feature 0 is informative
    logits = 2.0 * X[:, 0]
    p = 1 / (1 + np.exp(-logits))
    y = (rng.uniform(size=n) < p).astype(int)
    expr = pd.DataFrame(X, columns=[f"g{i}" for i in range(5)],
                        index=[f"s{i}" for i in range(n)])
    meta = pd.DataFrame({"risk_class": ["high_risk_OPMD" if v else "stable_OPMD"
                                        for v in y]},
                        index=expr.index)
    fit = fit_l1_logistic(expr, meta, list(expr.columns))
    assert fit.coefs["g0"] > 0
    assert abs(fit.coefs["g0"]) > max(abs(fit.coefs[g]) for g in ("g1", "g2", "g3", "g4"))
