import numpy as np
import pandas as pd
import pytest

from chorus.eval import LosoConfig, run_loso, run_within_study_cv, summarise
from chorus.features import InvariantSelector, VarianceSelector, meta_analyse, per_study_effects
from chorus.harmonise import make_harmoniser
from chorus.schema import ClinicalSchema, binary_label, validate_metadata
from chorus.simulate import simulate_cohorts, truth


@pytest.fixture(scope="module")
def sim():
    return simulate_cohorts(n_studies=5, n_per_study=[80, 100, 60, 120, 90], n_features=400, n_causal=20, n_shortcut=40, seed=1)


def test_schema_validation(sim):
    _, M = sim.pooled()
    assert validate_metadata(M, ClinicalSchema(), strict=False) == []
    y = binary_label(M)
    assert set(y.dropna().unique()) == {0.0, 1.0}


def test_harmonisers_addon_mode(sim):
    X, M = sim.pooled()
    s = M["study_id"]
    train = s != "S5"
    for name in ["identity", "center_scale", "quantile_ref", "fsqn", "rank", "combat_addon"]:
        H = make_harmoniser(name)
        H.fit(X[train], s[train], binary_label(M)[train])
        Xt = H.transform(X[~train], s[~train])  # study S5 never seen at fit time
        assert Xt.shape[0] == (~train).sum()
        assert np.isfinite(Xt.values).all()
    # per study centring removes study means
    H = make_harmoniser("center_scale").fit(X[train], s[train])
    Z = H.transform(X, s)
    assert np.allclose(Z.groupby(s).mean().values, 0, atol=1e-6)


def test_gene_pair_requires_labels(sim):
    X, M = sim.pooled()
    s = M["study_id"]
    with pytest.raises(ValueError):
        make_harmoniser("gene_pair").fit(X, s)
    H = make_harmoniser("gene_pair", n_pairs=20, n_candidate_features=40).fit(X, s, binary_label(M))
    assert H.transform(X, s).shape[1] == 20


def test_invariant_selector_recovers_causal(sim):
    X, M = sim.pooled()
    y, s = binary_label(M), M["study_id"]
    H = make_harmoniser("center_scale").fit(X, s)
    Xh = H.transform(X, s)
    t = meta_analyse(per_study_effects(Xh, y, s))
    assert {"effect", "I2", "Q_p", "sign_consistency"} <= set(t.columns)
    inv = InvariantSelector(k=20).fit(Xh, y, s)
    var = VarianceSelector(k=20).fit(Xh, y, s)
    causal = set(truth(sim)["causal"])
    assert len(set(inv.selected_) & causal) > len(set(var.selected_) & causal)
    assert len(set(inv.selected_) & causal) >= 12


def test_loso_harmonisation_helps(sim):
    cfg = LosoConfig(harmonisers=["identity", "center_scale"], selectors=[{"name": "variance", "k": 100}, {"name": "invariant", "k": 20}], classifiers=["logreg_en"], n_boot=50)
    res, preds, selected = run_loso(sim, cfg, verbose=False)
    assert "error" not in res.columns or res["error"].isna().all()
    summ = summarise(res)
    best_identity = summ.xs("identity", level="harmoniser")["mean"].max()
    best_harm = summ.xs("center_scale", level="harmoniser")["mean"].max()
    assert best_harm > best_identity
    inv = summ.loc[("center_scale", "invariant", "logreg_en"), "mean"]
    var = summ.loc[("center_scale", "variance", "logreg_en"), "mean"]
    assert inv > var
    assert preds["held_out"].nunique() == 5


def test_within_study_cv_runs(sim):
    cfg = LosoConfig(selectors=[{"name": "variance", "k": 50}], classifiers=["logreg_en"], n_boot=20)
    w = run_within_study_cv(sim, cfg, n_splits=3)
    assert set(w["study"]) == set(sim.study_ids)


def test_cross_study_ensemble(sim):
    cfg = LosoConfig(harmonisers=["center_scale"], selectors=[{"name": "invariant", "k": 20}], classifiers=["cross_study"], n_boot=20)
    res, _, _ = run_loso(sim, cfg, verbose=False)
    assert res["auroc"].mean() > 0.6


def test_nested_selection(sim):
    cfg = LosoConfig(harmonisers=["identity", "center_scale"], selectors=[{"name": "invariant", "k": 20}], classifiers=["logreg_en"], nested_select_harmoniser=True, n_boot=20)
    res, _, _ = run_loso(sim, cfg, verbose=False)
    nested = res[res["harmoniser"] == "nested"]
    assert len(nested) == 5
    assert set(nested["harmoniser_chosen"]) <= {"identity", "center_scale"}


def test_loader_log_transforms_linear_scale_author_matrix(tmp_path):
    from chorus.io import Cohort

    X = pd.DataFrame({"s1": [0.0, 120.5, 3.0], "s2": [10.0, 60.0, 7.25]}, index=["ENSG01.1", "ENSG02.3", "ENSG03.1"])
    X.to_csv(tmp_path / "m.tsv.gz", sep="\t")
    meta = pd.DataFrame({"sample_id": ["s1", "s2"], "subject_id": ["a", "b"], "lc_status": ["LC", "COVID_recovered"], "case_definition": "study_specific", "sample_matrix": "PBMC"})
    meta.to_csv(tmp_path / "meta.csv", index=False)
    c = Cohort.from_files("T", "bulk_rna", tmp_path / "m.tsv.gz", tmp_path / "meta.csv", is_counts=False, provenance="author_matrix")
    assert c.source["provenance"] == "author_normalised"
    assert np.isclose(c.X.loc["s1", "ENSG02"], np.log2(121.5))
    small = Cohort.from_files("T", "bulk_rna", tmp_path / "m.tsv.gz", tmp_path / "meta.csv", is_counts=True)
    assert small.X.loc["s1", "ENSG02"] == 120.5  # counts are left for log CPM
