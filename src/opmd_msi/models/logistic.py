"""L1-regularized logistic regression for the binary 'transformed vs stable' outcome."""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from sklearn.linear_model import LogisticRegressionCV


@dataclass(frozen=True)
class LogisticFit:
    coefs: pd.Series
    selected_C: float
    feature_columns: tuple[str, ...]
    sample_index: tuple[str, ...]
    train_X: pd.DataFrame
    train_y: pd.Series


def fit_l1_logistic(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    Cs: tuple[float, ...] = (0.01, 0.1, 1.0, 10.0),
    max_iter: int = 5000,
    target_col: str = "risk_class",
    case_value: str = "high_risk_OPMD",
    control_value: str = "stable_OPMD",
) -> LogisticFit:
    """Binary L1 logistic: high-risk OPMD vs stable OPMD."""
    relevant = metadata[metadata[target_col].isin([case_value, control_value])]
    if relevant.empty or relevant[target_col].nunique() < 2:
        return LogisticFit(
            coefs=pd.Series(0.0, index=genes), selected_C=float("nan"),
            feature_columns=tuple(genes), sample_index=(),
            train_X=pd.DataFrame(columns=genes), train_y=pd.Series(dtype=int),
        )
    X = expression.loc[relevant.index, genes].astype(float)
    y = (relevant[target_col] == case_value).astype(int)
    model = LogisticRegressionCV(
        Cs=list(Cs), penalty="l1", solver="liblinear",
        class_weight="balanced", max_iter=max_iter,
        cv=min(5, max(2, len(y) // 10)),
        random_state=0,
    )
    model.fit(X.values, y.values)
    coefs = pd.Series(model.coef_.ravel(), index=genes)
    return LogisticFit(
        coefs=coefs, selected_C=float(model.C_[0]),
        feature_columns=tuple(genes), sample_index=tuple(relevant.index),
        train_X=X, train_y=y,
    )
