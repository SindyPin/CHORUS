"""Release v1: every exported matrix reproduces the values the benchmark loads for the same cohort, samples and features.

For each matrix row of release/manifest.csv the source config is loaded with chorus.cli.load_config (or the paired
Su2022 layers with scripts/paired_layers.load_layers), the samples and features of the file are taken from the pooled
matrix, the harmoniser of the file (none, center_scale or fsqn) is fitted on those samples with the cohort as the group
and without labels, and the values are compared (7 significant digits are written). Skipped when the local data or the
release is absent.
"""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "release" / "manifest.csv"

pytestmark = pytest.mark.skipif(not MANIFEST.exists() or not (ROOT / "data" / "meta").exists(), reason="release or local data absent")


def _rows():
    if not MANIFEST.exists():
        return []
    m = pd.read_csv(MANIFEST)
    return m[m["kind"] == "matrix"].to_dict("records")


_CACHE = {}


def _pooled(config):
    from chorus.cli import load_config

    if config not in _CACHE:
        cs, _, _ = load_config(str(ROOT / config))
        _CACHE[config] = cs
    return _CACHE[config]


def _paired(layer):
    if "paired" not in _CACHE:
        spec = importlib.util.spec_from_file_location("paired_layers", ROOT / "scripts" / "paired_layers.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _CACHE["paired"] = mod.load_layers("Su2022")[0]
    return _CACHE["paired"][layer]


@pytest.mark.parametrize("row", _rows(), ids=lambda r: r["path"])
def test_release_matrix_matches_benchmark_values(row):
    from chorus.harmonise import make_harmoniser

    F = pd.read_csv(ROOT / "release" / row["path"], index_col=0)
    src = json.loads(row["source"])
    if src["loader"] == "paired":
        ref = _paired(src["layer"]).loc[F.index.astype(str), list(F.columns)].apply(pd.to_numeric, errors="coerce")
    else:
        X, M = _pooled(src["config"]).pooled(list(F.columns))
        ref = X.loc[F.index]
        h = src.get("harmoniser", "none")
        if h != "none":
            s = M.loc[F.index, "study_id"].astype(str)
            ref = make_harmoniser(h).fit(ref, s).transform(ref, s)
    assert ref.shape == F.shape and list(ref.columns) == list(F.columns)
    np.testing.assert_allclose(F.to_numpy(float), ref.to_numpy(float), rtol=2e-6, atol=1e-9, equal_nan=True)
