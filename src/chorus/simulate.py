"""Simulated multi-study omics data with study effects and unstable (study specific) signals.

The generator mimics the failure modes documented in the literature: study specific
location and scale shifts per feature (platform and pipeline), differing case to control
ratios, a small set of genes with an invariant causal effect, and a larger set of
"shortcut" genes whose association with the label exists only in some studies (for
example because they track hospitalisation, which is confounded with LC status in one
study but not another). A classifier trained on pooled unharmonised data latches on to the
shortcuts and the study offsets; leave-one-study-out exposes it.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .io import Cohort, CohortSet


def simulate_cohorts(
    n_studies: int = 6,
    n_per_study: Optional[List[int]] = None,
    n_features: int = 2000,
    n_causal: int = 30,
    n_shortcut: int = 60,
    causal_effect: float = 0.6,
    shortcut_effect: float = 1.2,
    study_shift_sd: float = 1.0,
    study_scale_sd: float = 0.3,
    case_fraction: Optional[List[float]] = None,
    noise_sd: float = 1.0,
    n_factors: int = 0,
    factor_sd: float = 0.0,
    seed: int = 0,
    layer: str = "bulk_rna",
) -> CohortSet:
    rng = np.random.default_rng(seed)
    n_per_study = n_per_study or [rng.integers(40, 200) for _ in range(n_studies)]
    case_fraction = case_fraction or list(np.clip(rng.normal(0.5, 0.15, n_studies), 0.2, 0.8))
    features = [f"G{i:05d}" for i in range(n_features)]
    causal = features[:n_causal]
    shortcut = features[n_causal : n_causal + n_shortcut]
    # causal genes have detectable effects of either sign: |beta| uniform on [0.4, 1.0] x causal_effect
    beta_causal = rng.uniform(0.4, 1.0, n_causal) * causal_effect * rng.choice([-1, 1], n_causal)
    cohorts = []
    for s in range(n_studies):
        n = int(n_per_study[s])
        y = (rng.random(n) < case_fraction[s]).astype(float)
        X = rng.normal(0, noise_sd, (n, n_features))
        # invariant causal signal, same direction and size in every study
        X[:, :n_causal] += np.outer(y, beta_causal)
        # shortcut signal: present with a study specific sign and size (zero in roughly a third of studies)
        active = rng.random(n_shortcut) > 0.33
        beta_short = rng.normal(0, shortcut_effect, n_shortcut) * active
        X[:, n_causal : n_causal + n_shortcut] += np.outer(y, beta_short)
        # sample level technical factors shared across genes (composition, library quality), with
        # study specific means so that they also act as study effects; off by default
        if n_factors and factor_sd > 0:
            loadings = rng.normal(0, factor_sd, (n_factors, n_features))
            scores = rng.normal(0, 1, (n, n_factors)) + rng.normal(0, 1, n_factors)
            X += scores @ loadings
        # study location and scale effects per feature (platform, pipeline, tissue mix)
        shift = rng.normal(0, study_shift_sd, n_features)
        scale = np.exp(rng.normal(0, study_scale_sd, n_features))
        X = X * scale + shift + 8.0
        sid = f"S{s+1}"
        idx = [f"{sid}_{i:04d}" for i in range(n)]
        meta = pd.DataFrame(
            {
                "sample_id": idx,
                "subject_id": idx,
                "study_id": sid,
                "lc_status": np.where(y == 1, "LC", "COVID_recovered"),
                "case_definition": rng.choice(["WHO_2021", "NICE_NG188", "CDC_4wk"]),
                "sample_matrix": "whole_blood",
                "omics_layer": layer,
                "age": rng.normal(50, 12, n).round(),
                "sex": rng.choice(["female", "male"], n),
                "days_since_infection": rng.integers(90, 400, n),
            },
            index=idx,
        )
        cohorts.append(Cohort(sid, layer, pd.DataFrame(X, index=idx, columns=features), meta, "simulated", False, {"accession": "simulated", "causal": ",".join(causal)}))
    return CohortSet(cohorts)


def truth(cs: CohortSet) -> Dict[str, List[str]]:
    c = cs.cohorts[0]
    return {"causal": c.source["causal"].split(",")}
