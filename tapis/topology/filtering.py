from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


@dataclass
class ReadQualityFilter:
    model: Pipeline
    feature_names: tuple[str, ...]

    @classmethod
    def fit(
        cls,
        features: NDArray[np.float64],
        labels: NDArray[np.int_],
        feature_names: Sequence[str],
    ) -> "ReadQualityFilter":
        pipeline = Pipeline(
            [
                ("scaler", StandardScaler()),
                ("classifier", LogisticRegression(max_iter=500, random_state=0)),
            ]
        )
        pipeline.fit(features, labels)
        return cls(model=pipeline, feature_names=tuple(feature_names))

    def score(self, features: NDArray[np.float64]) -> NDArray[np.float64]:
        probabilities = self.model.predict_proba(features)
        return np.asarray(probabilities[:, 1], dtype=float)

    def predict_keep(self, features: NDArray[np.float64]) -> NDArray[np.bool_]:
        probabilities = self.score(features)
        return probabilities >= 0.5
