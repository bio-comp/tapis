from __future__ import annotations

import numpy as np

from tapis.topology.filtering import ReadQualityFilter


def test_read_quality_filter_fit_and_score() -> None:
    features = np.array(
        [
            [0.1, 0.2],
            [0.2, 0.1],
            [1.0, 1.1],
            [1.2, 1.0],
        ]
    )
    labels = np.array([0, 0, 1, 1])
    model = ReadQualityFilter.fit(features, labels, feature_names=["edr", "mapq"])

    scores = model.score(features)
    keep_flags = model.predict_keep(features)

    assert scores.shape == (4,)
    assert keep_flags.dtype == bool
    assert keep_flags.tolist().count(True) >= 2
