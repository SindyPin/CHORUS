import dataclasses
from types import SimpleNamespace

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


def test_loader_rounds_estimated_counts_and_logs_normalised_values(tmp_path):
    from chorus.io import Cohort

    X = pd.DataFrame({"s1": [0.0, 120.4, 3.0], "s2": [10.0, 60.0, 7.6]}, index=["ENSG01.1", "ENSG02.3", "ENSG03.1"])
    X.to_csv(tmp_path / "m.tsv.gz", sep="\t")
    meta = pd.DataFrame({"sample_id": ["s1", "s2"], "subject_id": ["a", "b"], "lc_status": ["LC", "COVID_recovered"], "case_definition": "study_specific", "sample_matrix": "PBMC"})
    meta.to_csv(tmp_path / "meta.csv", index=False)
    est = Cohort.from_files("T", "bulk_rna", tmp_path / "m.tsv.gz", tmp_path / "meta.csv", is_counts=True, provenance="author_matrix")
    assert est.source["provenance"] == "author_estimated_counts"
    assert est.X.loc["s1", "ENSG02"] == 120 and est.X.loc["s2", "ENSG03"] == 8
    pd.DataFrame({"s1": [0, 120, 3], "s2": [10, 60, 7]}, index=X.index).to_csv(tmp_path / "int.tsv", sep="\t")
    raw = Cohort.from_files("T", "bulk_rna", tmp_path / "int.tsv", tmp_path / "meta.csv", is_counts=True, provenance="author_matrix")
    assert raw.source["provenance"] == "author_matrix"
    ours = Cohort.from_files("T", "bulk_rna", tmp_path / "m.tsv.gz", tmp_path / "meta.csv", is_counts=True, provenance="requantify")
    assert ours.source["provenance"] == "requantify" and ours.X.loc["s1", "ENSG02"] == 120
    norm = Cohort.from_files("T", "bulk_rna", tmp_path / "m.tsv.gz", tmp_path / "meta.csv", is_counts=False)
    assert norm.source["provenance"] == "author_normalised"
    assert np.isclose(norm.X.loc["s1", "ENSG02"], np.log2(121.4))


def test_train_only_cohort_is_never_held_out():
    cs = simulate_cohorts(n_studies=4, n_per_study=[60, 60, 60, 60], n_features=200, n_causal=10, n_shortcut=20, seed=2)
    cs.cohorts[0].source["train_only"] = "true"
    cfg = LosoConfig(harmonisers=["center_scale"], selectors=[{"name": "variance", "k": 50}], classifiers=["logreg"], n_boot=10)
    res, _, _ = run_loso(cs, cfg, verbose=False)
    assert cs.study_ids[0] not in set(res["held_out"])
    assert res["n_train_studies"].eq(3).all()  # still used for training


def test_one_sample_per_subject_keeps_earliest_post_acute():
    from chorus.io import Cohort, CohortSet

    meta = pd.DataFrame(
        {
            "sample_id": ["a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1"],
            "subject_id": ["a", "a", "a", "b", "b", "c", "c", "d"],
            "timepoint": ["Acute", "3_6M", ">12M", "3", "1", "7/29/2021", "4/2/2021", "acute"],
            "days_since_infection": [27, 140, 400, None, None, None, None, 10],
            "lc_status": "LC",
            "case_definition": "study_specific",
            "sample_matrix": "whole_blood",
        }
    ).set_index("sample_id", drop=False)
    X = pd.DataFrame(np.arange(16.0).reshape(8, 2), index=meta.index, columns=["g1", "g2"])
    cs = CohortSet([Cohort("T", "bulk_rna", X, meta)]).one_sample_per_subject()
    kept = set(cs["T"].X.index)
    assert kept == {"a2", "b2", "c2"}  # acute draws dropped, d has no post acute sample
    assert cs["T"].meta.index.equals(cs["T"].X.index)


@pytest.mark.parametrize("dtype", [object, "string", "str", "category"])
def test_covariate_matrix_one_hot_encodes_text_columns(dtype):
    from chorus.eval.loso import _covariate_matrix

    meta = pd.DataFrame({"sample_matrix": pd.Series(["PBMC", "whole_blood", "PBMC"], dtype=dtype), "age": [30, 40, 50]})
    cov = _covariate_matrix(meta, ["age", "sample_matrix"])
    assert list(cov.columns) == ["cov_age", "cov_sample_matrix_PBMC", "cov_sample_matrix_whole_blood"]
    assert cov["cov_sample_matrix_PBMC"].tolist() == [1.0, 0.0, 1.0]


def test_gene_list_ignores_inline_comments(tmp_path):
    from chorus.features.priors import load_gene_list

    (tmp_path / "g.txt").write_text("# header\nENSG00000141510  # TP53\nENSG00000005339.12 # CREBBP\n\n")
    assert load_gene_list(tmp_path / "g.txt") == ["ENSG00000141510", "ENSG00000005339"]


def test_chorus_operator_partition_and_addon():
    from chorus.harmonise import make_harmoniser

    cs = simulate_cohorts(n_studies=5, n_per_study=[100, 90, 110, 80, 100], n_features=800, n_causal=20, n_shortcut=40, case_fraction=[0.2, 0.8, 0.5, 0.3, 0.7], n_factors=2, factor_sd=0.5, seed=7)
    X, M = cs.pooled()
    y, s = binary_label(M), M["study_id"]
    train = s != "S5"
    H = make_harmoniser("chorus").fit(X[train], s[train], y[train])
    causal, short = set(truth(cs)["causal"]), set(X.columns[20:60])
    assert len(set(H.controls_) & causal) <= 1 and len(set(H.controls_) & short) == 0  # weak causal genes can pass the null screen
    assert len(set(H.invariant_) & short) <= 2 and len(set(H.invariant_) & causal) >= 10  # a shortcut active with one sign in every training study is indistinguishable
    assert len(set(H.unstable_) & short) >= 30
    Z = H.transform(X[~train], s[~train])  # unseen study, label free
    assert Z.shape == (int((~train).sum()), X.shape[1]) and np.isfinite(Z.values).all()
    with pytest.raises(ValueError):
        make_harmoniser("chorus").fit(X, s)


def test_chorus_operator_competitive_under_imbalance():
    cs = simulate_cohorts(n_studies=5, n_per_study=[100, 90, 110, 80, 100], n_features=800, n_causal=20, n_shortcut=40, case_fraction=[0.2, 0.8, 0.5, 0.3, 0.7], n_factors=2, factor_sd=0.5, seed=7)
    cfg = LosoConfig(harmonisers=["identity", "combat_addon", "chorus"], selectors=[{"name": "invariant", "k": 30}], classifiers=["logreg"], n_boot=20)
    res, _, _ = run_loso(cs, cfg, verbose=False)
    m = res.groupby("harmoniser")["auroc"].mean()
    assert m["chorus"] >= m["identity"] - 0.02 and m["chorus"] >= m["combat_addon"] - 0.02
    assert res["sign_agreement"].notna().any()


def test_taco_without_tabpfn_when_disabled(monkeypatch):
    from chorus.models import make_classifier

    monkeypatch.setenv("CHORUS_DISABLE_TABPFN", "1")
    clf = make_classifier("taco", 0)
    assert [n for n, _ in clf.estimators] == ["rf", "extratrees", "svc", "mlp"]


def test_prior_selectors_keep_distinct_labels(sim, tmp_path):
    cs = sim
    X = cs.pooled()[0]
    a, b = tmp_path / "list_a.txt", tmp_path / "list_b.txt"
    a.write_text("\n".join(X.columns[:15]) + "\n")
    b.write_text("\n".join(X.columns[15:25]) + "\n")
    sels = [{"name": "prior", "genes": str(a), "mode": "prior_only"}, {"name": "prior", "genes": str(a), "mode": "union", "k": 20},
            {"name": "prior", "genes": str(b), "mode": "union", "k": 20}]
    from chorus.eval.loso import selector_label
    from chorus.features.priors import load_gene_list

    sels = [{**c, "label": selector_label(c), "genes": load_gene_list(c["genes"])} for c in sels]  # as chorus.cli does
    cfg = LosoConfig(harmonisers=["identity"], selectors=sels, classifiers=["logreg"], n_boot=10)
    res, _, selected = run_loso(cs, cfg, verbose=False)
    labels = {"prior:list_a:prior_only", "prior:list_a:union:k20", "prior:list_b:union:k20"}
    assert set(res["selector"]) == labels
    assert {k[2] for k in selected} == labels
    assert (summarise(res)["n_studies"] == len(cs.study_ids)).all()


def _single_study(n_pos, n_neg, shift, seed=0):
    from chorus.io import Cohort, CohortSet

    rng = np.random.default_rng(seed)
    n = n_pos + n_neg
    ids = [f"s{i}" for i in range(n)]
    meta = pd.DataFrame({"sample_id": ids, "subject_id": ids, "study_id": "P", "lc_status": ["LC"] * n_pos + ["COVID_recovered"] * n_neg,
                         "case_definition": "study_specific", "sample_matrix": "PBMC"}).set_index("sample_id", drop=False)
    X = pd.DataFrame(rng.normal(size=(n, 200)), index=meta.index, columns=[f"g{j}" for j in range(200)])
    X.iloc[:n_pos, :10] += shift  # planted signal in the first 10 genes
    return CohortSet([Cohort("P", "bulk_rna", X, meta)])


def test_within_study_cv_per_fold_mean_recovers_planted_signal():
    cfg = LosoConfig(selectors=[{"name": "differential", "k": 20}], classifiers=["logreg"], n_boot=10)
    w, folds = run_within_study_cv(_single_study(20, 12, shift=1.5), cfg, n_splits=5, n_repeats=4, return_folds=True)
    row = w.iloc[0]
    assert row["n_folds"] == 20 and row["n_failed_folds"] == 0 and len(folds) == 20
    assert row["auroc"] > 0.9 and np.isclose(row["auroc"], folds["auroc"].mean())
    assert {"auroc_sd", "auroc_pooled", "auroc_pooled_repeat0"} <= set(w.columns)
    null = run_within_study_cv(_single_study(20, 12, shift=0.0, seed=1), cfg, n_splits=5, n_repeats=4).iloc[0]
    assert 0.3 < null["auroc"] < 0.7



def _sex_logcpm(n_f=20, n_m=12, seed=0):
    from chorus.sex import RPS4Y1, XIST

    rng = np.random.default_rng(seed)
    xist = np.r_[rng.normal(6, 0.5, n_f), rng.normal(-3, 0.5, n_m)]
    rps = np.r_[rng.normal(-4, 0.5, n_f), rng.normal(7, 0.5, n_m)]
    idx = [f"F{i}" for i in range(n_f)] + [f"M{i}" for i in range(n_m)]
    return pd.DataFrame({XIST: xist, RPS4Y1: rps, "ENSG_OTHER": rng.normal(size=n_f + n_m)}, index=idx)


def test_infer_sex_splits_two_components_and_flags_outliers():
    from chorus.sex import RPS4Y1, XIST, infer_sex

    L = _sex_logcpm()
    L.loc["odd"] = {XIST: -0.2, RPS4Y1: -0.4, "ENSG_OTHER": 0.0}  # neither XIST nor RPS4Y1, as GSM7069947
    t = infer_sex(L)
    assert (t.loc[[i for i in L.index if i.startswith("F")], "inferred_sex"] == "female").all()
    assert (t.loc[[i for i in L.index if i.startswith("M")], "inferred_sex"] == "male").all()
    assert t.loc["odd", "inferred_sex"] == "ambiguous" and t.loc["odd", "reason"] == "outlier"
    assert (t["inferred_sex_method"] == "gmm2:XIST+RPS4Y1").all()
    # without an XIST row the split uses RPS4Y1 alone and says so
    t1 = infer_sex(_sex_logcpm().drop(columns=XIST))
    assert (t1["inferred_sex"] == np.where(t1.index.str.startswith("F"), "female", "male")).all()
    assert t1["inferred_sex_method"].str.startswith("gmm2:RPS4Y1").all()


def test_sex_chromosome_genes_include_pseudoautosomal(tmp_path):
    from chorus.sex import gene_chromosomes, sex_chromosome_genes

    p = tmp_path / "tx2gene.tsv"
    pd.DataFrame(
        {
            "transcript_id": ["T1", "T2", "T3", "T4", "T5"],
            "gene_id": ["ENSG01.3", "ENSG02.1", "ENSG03.2", "ENSG03.2_PAR_Y", "ENSG04.1"],
            "chrom": ["chr1", "chrX", "chrX", "chrY", "chrY"],
        }
    ).to_csv(p, sep="\t", index=False)
    ch = gene_chromosomes(p)
    assert ch["ENSG03"] == "chrX,chrY"
    assert sex_chromosome_genes(ch, ["ENSG01", "ENSG02", "ENSG03", "ENSG04", "ENSG99"]) == ["ENSG02", "ENSG03", "ENSG04"]


def test_harmonisation_metrics_on_shifted_studies():
    from chorus.eval.harmonisation import MARKERS, harmonisation_table
    from chorus.sex import RPS4Y1, XIST

    cs = simulate_cohorts(n_studies=3, n_per_study=[40, 50, 60], n_features=120, n_causal=10, n_shortcut=10, seed=3)
    X, M = cs.pooled()
    y, s = binary_label(M), M["study_id"]
    rng = np.random.default_rng(0)
    female = pd.Series(rng.random(len(X)) < 0.5, index=X.index)
    X[XIST] = np.where(female, 6.0, -3.0) + rng.normal(0, 0.3, len(X)) + s.map({"S1": 0.0, "S2": 2.0, "S3": -1.0}).values
    X[RPS4Y1] = np.where(female, -4.0, 7.0) + rng.normal(0, 0.3, len(X))
    neut = list(MARKERS["neutrophil"].values())
    base = rng.normal(size=len(X))
    for g in neut:
        X[g] = base + rng.normal(0, 0.2, len(X))
    feats = [c for c in X.columns if c.startswith("g") or c.startswith("f") or c not in [XIST, RPS4Y1] + neut]
    sex = pd.Series(np.where(female, "female", "male"), index=X.index)
    Xa = make_harmoniser("center_scale").fit(X, s).transform(X, s)
    t = harmonisation_table(X, Xa, s, y, feats, sex=sex, n_pcs=10, k=10).set_index(["metric", "scope"])
    assert t.loc[("silhouette_study", "pool"), "after"] < t.loc[("silhouette_study", "pool"), "before"]
    assert t.loc[("knn_other_study", "pool"), "after"] > t.loc[("knn_other_study", "pool"), "before"]
    # a per study affine map per gene leaves the within study effect unchanged
    assert np.isclose(t.loc[("lc_effect_spearman", "mean_over_studies"), "before_vs_after"], 1.0)
    assert t.loc[("sex_auroc_XIST", "pool"), "after"] > 0.99 and t.loc[("sex_auroc_RPS4Y1", "pool"), "before"] > 0.99
    assert t.loc[("marker_score_spearman_neutrophil", "pool"), "n_genes"] == len(neut)
    ident = harmonisation_table(X, X, s, y, feats, sex=sex, n_pcs=10, k=10).set_index(["metric", "scope"])
    assert np.isclose(ident.loc[("marker_score_spearman_neutrophil", "pool"), "before_vs_after"], 1.0)
    assert np.isclose(ident.loc[("silhouette_study", "pool"), "before"], ident.loc[("silhouette_study", "pool"), "after"])


def test_harmonisation_diagnostics_write_one_table_per_fold_and_harmoniser(sim, tmp_path):
    from chorus.eval.harmonisation import run_harmonisation_diagnostics

    cfg = LosoConfig(harmonisers=["identity", "center_scale"])
    feats = list(sim.cohorts[0].features[:150])
    res = run_harmonisation_diagnostics(sim, cfg, feats, out=tmp_path, n_pcs=10, k=10, verbose=False)
    folds = sorted(res["held_out"].unique())
    assert len(folds) == len(sim.study_ids)
    for f in folds:
        for h in cfg.harmonisers:
            assert (tmp_path / f / f"{h}.csv").exists()
    assert (tmp_path / "all_metrics.csv").exists()
    assert "error" not in set(res["family"])

def _load_script(name):
    return _load_script_path(name)


def _load_script_path(name):
    """Load scripts/<name>.py (``name`` may hold a subfolder, ``metabolomics/build_cohorts``)."""
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name.replace("/", "_"), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_concordance_effects_and_spearman():
    conc = _load_script("concordance")
    rng = np.random.default_rng(0)
    genes = [f"g{j}" for j in range(300)]
    effect = np.r_[np.linspace(3.0, 1.0, 30), np.zeros(270)]  # graded, so ranks within the top genes carry signal
    Xs, Ms = [], []
    for sid, sign in [("A", 1), ("B", 1), ("C", -1)]:
        n = 40
        y = np.r_[np.ones(n // 2), np.zeros(n // 2)]
        X = rng.normal(size=(n, 300)) + sign * np.outer(y, effect)
        idx = [f"{sid}{i}" for i in range(n)]
        Xs.append(pd.DataFrame(X, index=idx, columns=genes))
        Ms.append(pd.DataFrame({"study_id": sid, "lc_status": np.where(y == 1, "LC", "COVID_recovered"), "case_definition": "WHO_2021",
                                "sample_matrix": "PBMC", "time_anchor": "positive_test", "days_since_infection": 120}, index=idx))
    X, M = pd.concat(Xs), pd.concat(Ms)
    eff = conc.cohort_effects(X, M)
    assert list(eff.columns) == ["A", "B", "C"]
    r, p = conc.spearman_matrix(eff)
    assert r.loc["A", "B"] > 0.2 and r.loc["A", "C"] < -0.2 and np.allclose(np.diag(r), 1.0)
    top = conc.strongest_genes(eff, 30)
    assert len(set(top.index) & set(genes[:30])) >= 27
    rt, _ = conc.spearman_matrix(eff.loc[top.index])
    assert rt.loc["A", "B"] > r.loc["A", "B"]
    ann = conc.cohort_annotation(M).set_index("study_id")
    assert ann.loc["A", "case_definition"] == "WHO_2021" and ann.loc["A", "negative_class_lc_status"] == "COVID_recovered"
    assert ann.loc["A", "control_type"] == "missing"  # no schema control_type column in M
    M2 = M.assign(control_type=np.where(M["lc_status"] == "LC", "not_applicable", "convalescent_unverified"))
    assert conc.cohort_annotation(M2).set_index("study_id").loc["A", "control_type"] == "convalescent_unverified"
    assert ann.loc["A", "time_anchor"] == "positive_test" and ann.loc["A", "n_LC"] == 20


def test_concordance_min_shared_and_pair_flags():
    conc = _load_script("concordance")
    rng = np.random.default_rng(1)
    E = pd.DataFrame(rng.normal(size=(40, 3)), columns=["A", "B", "C"])
    E.loc[4:, "C"] = np.nan  # C shares 4 features with A and B
    r, p = conc.spearman_matrix(E, min_shared=10)
    assert np.isnan(r.loc["A", "C"]) and np.isnan(p.loc["B", "C"]) and not np.isnan(r.loc["A", "B"])
    r0, _ = conc.spearman_matrix(E)
    assert not np.isnan(r0.loc["A", "C"])  # without the minimum the 4 shared features give a value
    n = conc.shared_counts(E)
    t = conc.pair_table(r0, p, n, "all", 10).set_index(["cohort_a", "cohort_b"])
    assert t.loc[("A", "C"), "n_shared"] == 4 and t.loc[("A", "C"), "fewer_than_10_shared"]
    assert not t.loc[("A", "B"), "fewer_than_10_shared"] and len(t) == 3


def test_paired_layer_pairs_flag():
    from chorus.eval.paired import layer_pairs

    E = pd.DataFrame({"a": np.arange(20.0), "b": np.r_[np.arange(5.0), [np.nan] * 15]})
    t = layer_pairs(E)
    assert t.loc[0, "n_shared"] == 5 and bool(t.loc[0, "fewer_than_10_shared"])


def test_enrichment_hypergeometric_and_mapping():
    from scipy import stats

    enr = _load_script("enrichment")
    universe = [f"G{i}" for i in range(100)]
    selected = universe[:10] + ["NOT_IN_UNIVERSE"]
    gene_list = universe[5:15] + ["ALSO_OUT"]
    r = enr.hypergeom_enrichment(selected, gene_list, universe)
    assert (r["universe"], r["list_in_universe"], r["n_selected"], r["overlap"]) == (100, 10, 10, 5)
    assert np.isclose(r["expected"], 1.0) and np.isclose(r["fold_enrichment"], 5.0)
    assert np.isclose(r["p_hypergeom"], stats.hypergeom.sf(4, 100, 10, 10))
    q = enr.bh(pd.Series([0.01, 0.04, 0.03, np.nan]))
    assert np.allclose(q.iloc[:3], [0.03, 0.04, 0.04]) and np.isnan(q.iloc[3])
    hgnc = pd.DataFrame({"symbol": ["C1QA", "C4A", "X"], "ensembl_gene_id": ["ENSG1", "ENSG2", "ENSG3"], "uniprot_ids": ["P02745", "P0C0L4|P0C0L5", np.nan]})
    assert enr.uniprot_to_ensembl(["P0C0L4", "P02745", "Q99999"], hgnc) == ["ENSG1", "ENSG2"]


def test_combat_addon_protected_covariate_constant_within_study():
    rng = np.random.default_rng(0)
    n, p = 40, 50
    idx = [f"a{i}" for i in range(n)] + [f"b{i}" for i in range(n)]
    s = pd.Series(["A"] * n + ["B"] * n, index=idx)
    X = pd.DataFrame(rng.normal(size=(2 * n, p)) + np.r_[np.zeros(n), np.full(n, 5.0)][:, None], index=idx, columns=[f"g{j}" for j in range(p)])
    y = pd.Series(rng.integers(0, 2, 2 * n).astype(float), index=idx)
    protect = pd.DataFrame({"m_A": (s == "A").astype(float), "m_B": (s == "B").astype(float)})  # one hot by study, as sample_matrix
    means = {}
    for name, prot in [("protect", protect), ("none", None)]:
        H = make_harmoniser("combat_addon", protect=prot) if prot is not None else make_harmoniser("combat_addon")
        Z = H.fit(X, s, y).transform(X, s)
        means[name] = Z.groupby(s).mean().mean(axis=1)
    for name in means:
        assert abs(means[name]["A"] - means[name]["B"]) < 0.3, means
        assert np.allclose(means[name], 2.5, atol=0.3), means
    assert np.allclose(means["protect"], means["none"])


def test_sex_control_forced_genes_and_incomplete_feature_drop(sim):
    from chorus.io import Cohort, CohortSet
    from chorus.sex import RPS4Y1, XIST

    ctl = _load_script("control_runs")
    rng = np.random.default_rng(0)
    cohorts = []
    for i, c in enumerate(sim.cohorts):
        X = c.X.iloc[:, :60].copy()
        meta = c.meta.copy()
        female = rng.random(len(X)) < 0.5
        meta["inferred_sex"] = np.where(female, "female", "male")
        if i != 0:  # the first cohort has no XIST row, as GSE270045
            X[XIST] = np.where(female, 6.0, -3.0) + rng.normal(0, 0.3, len(X))
        X[RPS4Y1] = np.where(female, -4.0, 7.0) + rng.normal(0, 0.3, len(X))
        X["YEXPR"] = np.where(female, -5.0, 4.0)  # chrY gene expressed in males everywhere
        X["YLOW"] = np.where(female, -5.0, -2.0 if i == 1 else 4.0)  # not expressed in the males of cohort 1
        cohorts.append(Cohort(c.study_id, c.layer, X, meta, c.feature_space, False, dict(c.source)))
    chroms = pd.Series({"YEXPR": "chrY", "YLOW": "chrY", RPS4Y1: "chrY", XIST: "chrX"})
    cs, cand, presence = ctl.forced_sex_genes(CohortSet(cohorts), chroms)
    assert set(cand) == {XIST, RPS4Y1, "YEXPR", "YLOW"}
    p = presence.set_index(["study_id", "gene_id"])["used"]
    assert not p[(sim.study_ids[0], XIST)] and p[(sim.study_ids[1], XIST)]
    assert not p[(sim.study_ids[1], "YLOW")] and p[(sim.study_ids[0], "YLOW")] and p.xs("YEXPR", level="gene_id").all()
    assert cs[sim.study_ids[1]].X["YLOW"].isna().all()
    cfg = LosoConfig(harmonisers=["identity"], selectors=[{"name": "variance", "k": 10}], classifiers=["logreg"], n_boot=10, drop_incomplete_features=True)
    feats = list(sim.cohorts[0].features[:60]) + [XIST, RPS4Y1, "YEXPR", "YLOW"]
    res, _, _ = run_loso(cs, cfg, features=feats, verbose=False)
    assert "error" not in res.columns or res["error"].isna().all()
    assert (res["n_features_dropped_incomplete"] == 2).all()  # XIST and YLOW are missing in some cohort of every fold


def test_parallel_analysis_counts_factors():
    from chorus.harmonise.anchored import parallel_analysis

    rng = np.random.default_rng(0)
    noise = rng.normal(size=(150, 300))
    assert parallel_analysis(noise) == 0
    signal = noise + rng.normal(size=(150, 3)) @ rng.normal(0, 1.0, size=(3, 300))
    assert parallel_analysis(signal) == 3


def test_chorus_zero_factors_and_protected_marker_scores():
    from chorus.harmonise.anchored import BLOOD_MARKERS

    rng = np.random.default_rng(1)
    n_s, n, p = 3, 60, 400
    genes = [f"g{j}" for j in range(p)]
    s = pd.Series(np.repeat([f"S{i}" for i in range(n_s)], n), index=[f"x{i}" for i in range(n_s * n)])
    y = pd.Series(rng.integers(0, 2, len(s)).astype(float), index=s.index)
    # no study shift and no factor: parallel analysis on the null block finds none
    X0 = pd.DataFrame(rng.normal(size=(len(s), p)), index=s.index, columns=genes)
    H0 = make_harmoniser("chorus").fit(X0, s, y)
    assert H0.k_ == 0 and np.isfinite(H0.transform(X0, s).values).all()
    assert make_harmoniser("chorus", k_method="studies_plus").fit(X0, s, y).k_ == n_s + 4
    # a composition score that varies within studies and moves most genes, plus one technical factor
    neut = list(BLOOD_MARKERS["neutrophil"].values())
    comp = rng.normal(size=len(s))
    tech = rng.normal(size=len(s))
    L = rng.normal(0, 0.8, size=p)
    T = rng.normal(0, 0.8, size=p)
    X = pd.DataFrame(np.outer(comp, L) + np.outer(tech, T) + rng.normal(size=(len(s), p)), index=s.index, columns=genes)
    for g in neut:
        X[g] = comp + rng.normal(0, 0.3, len(s))

    def score(Z):
        z = (Z[neut] - Z[neut].mean()) / Z[neut].std()
        return z.mean(axis=1)

    rho = {}
    for name, kw in [("plain", {}), ("protected", {"protect_scores": ["neutrophil"]})]:
        H = make_harmoniser("chorus", **kw).fit(X, s, y)
        Z = H.transform(X, s)
        rho[name] = np.mean([pd.Series(score(X)[s == k]).corr(score(Z)[s == k], method="spearman") for k in s.unique()])
        if name == "protected":
            assert not set(neut) & set(H.controls_) and H.score_coef_.shape == (1, X.shape[1])
    assert rho["protected"] > 0.9 and rho["protected"] > rho["plain"] + 0.2, rho
    H = make_harmoniser("chorus", protect_scores={"custom": neut[:2]}).fit(X, s, y)
    assert list(H.score_genes_) == ["custom"]


def test_rerun_of_some_harmonisers_replaces_only_their_rows(sim, tmp_path):
    import json

    from chorus.cli import _merge_rerun

    cfg = LosoConfig(harmonisers=["identity", "center_scale"], selectors=[{"name": "variance", "k": 20}], classifiers=["logreg"], n_boot=10, nested_select_harmoniser=True)
    feats = list(sim.cohorts[0].features[:80])
    res, preds, sel = run_loso(sim, cfg, features=feats, verbose=False)
    res.to_csv(tmp_path / "loso_results.csv", index=False)
    preds.to_csv(tmp_path / "loso_predictions.csv", index=False)
    (tmp_path / "selected_features.json").write_text(json.dumps({"|".join(k): v for k, v in sel.items()}))
    cfg.outer_harmonisers = ["center_scale", "nested"]
    res2, preds2, sel2 = run_loso(sim, cfg, features=feats, verbose=False)
    assert set(res2["harmoniser"]) == {"center_scale", "nested"}
    res2["auroc"] = -1.0  # mark the rerun rows
    r, p, s = _merge_rerun(tmp_path, cfg.outer_harmonisers, res2, preds2, {"|".join(k): v for k, v in sel2.items()})
    assert len(r) == len(res) and len(p) == len(preds) and set(s) == {"|".join(k) for k in sel}
    assert (r.loc[r["harmoniser"].isin(["center_scale", "nested"]), "auroc"] == -1).all()
    assert (r.loc[r["harmoniser"] == "identity", "auroc"].values == res.loc[res["harmoniser"] == "identity", "auroc"].values).all()
    assert list(r[["held_out", "harmoniser"]].itertuples(index=False)) == list(res[["held_out", "harmoniser"]].itertuples(index=False))


def test_load_config_extra_features_deduplicated(sim, tmp_path):
    import yaml

    from chorus.cli import load_config

    cohorts = []
    g = sim.cohorts[0].features[0]
    for i, c in enumerate(sim.cohorts):
        X = c.X.drop(columns=g) if i == 1 else c.X  # g fails the shared filter, so it can only enter as an extra
        X.T.to_csv(tmp_path / f"{c.study_id}.tsv", sep="\t")
        c.meta.to_csv(tmp_path / f"{c.study_id}_meta.csv", index=False)
        cohorts.append({"study_id": c.study_id, "matrix": str(tmp_path / f"{c.study_id}.tsv"), "meta": str(tmp_path / f"{c.study_id}_meta.csv"), "is_counts": False})
    (tmp_path / "c.yaml").write_text(yaml.safe_dump({"layer": "bulk_rna", "cohorts": cohorts, "preprocessing": {"deduplicate": False}, "benchmark": {}}))
    cs, _, cfg = load_config(str(tmp_path / "c.yaml"), extra_features=["NOT_THERE", g, g])
    assert g not in cfg["filtered_features"] and list(cs.cohorts[0].X.columns).count(g) == 1
    assert cs.cohorts[1].X[g].isna().all()


def test_load_config_exclude_samples(sim, tmp_path):
    import yaml

    from chorus.cli import load_config

    cohorts = []
    for c in sim.cohorts:
        c.X.T.to_csv(tmp_path / f"{c.study_id}.tsv", sep="\t")
        c.meta.to_csv(tmp_path / f"{c.study_id}_meta.csv", index=False)
        cohorts.append({"study_id": c.study_id, "matrix": str(tmp_path / f"{c.study_id}.tsv"), "meta": str(tmp_path / f"{c.study_id}_meta.csv"), "is_counts": False})
    drop = list(map(str, sim.cohorts[0].X.index[:3]))
    cohorts[0]["exclude_samples"] = drop
    (tmp_path / "c.yaml").write_text(yaml.safe_dump({"layer": "bulk_rna", "cohorts": cohorts, "preprocessing": {"deduplicate": False}, "benchmark": {}}))
    cs, _, _ = load_config(str(tmp_path / "c.yaml"))
    assert cs.cohorts[0].n == sim.cohorts[0].n - 3 and not set(drop) & set(map(str, cs.cohorts[0].X.index))
    assert cs.cohorts[1].n == sim.cohorts[1].n
    cohorts[0]["exclude_samples"] = ["NOT_A_SAMPLE"]
    (tmp_path / "c.yaml").write_text(yaml.safe_dump({"layer": "bulk_rna", "cohorts": cohorts, "preprocessing": {"deduplicate": False}, "benchmark": {}}))
    with pytest.raises(ValueError, match="exclude_samples"):
        load_config(str(tmp_path / "c.yaml"))


def test_combat_addon_training_study_means_equal_after_transform():
    # study offsets per gene, different case fractions, label effect on half of the genes and a
    # protected covariate that varies within studies: after fit and transform every training
    # study has the same per gene mean (the pool mean) and the same per gene variance
    rng = np.random.default_rng(1)
    p, frac, sizes = 60, [0.2, 0.5, 0.8], [30, 50, 40]
    blocks, studies, ys, cov = [], [], [], []
    for k, (n, f) in enumerate(zip(sizes, frac)):
        y = (rng.random(n) < f).astype(float)
        c = rng.normal(size=n)
        B = rng.normal(size=(n, p)) * np.exp(rng.normal(0, 0.4, p)) + rng.normal(0, 2, p)
        B[:, : p // 2] += 1.5 * y[:, None]
        B += 0.8 * c[:, None]
        blocks.append(B)
        studies += [f"S{k}"] * n
        ys.append(y)
        cov.append(c)
    idx = [f"x{i}" for i in range(sum(sizes))]
    X = pd.DataFrame(np.vstack(blocks), index=idx, columns=[f"g{j}" for j in range(p)])
    s = pd.Series(studies, index=idx)
    y = pd.Series(np.concatenate(ys), index=idx)
    prot = pd.DataFrame({"c": np.concatenate(cov)}, index=idx)
    for H in [make_harmoniser("combat_addon"), make_harmoniser("combat_addon", protect=prot)]:
        Z = H.fit(X, s, y).transform(X, s)
        mu = Z.groupby(s).mean()
        assert np.allclose(mu.values, Z.mean().values[None, :], atol=1e-8), (mu - Z.mean()).abs().max().max()
    # shrinkage, when asked for, keeps part of the offsets
    Z = make_harmoniser("combat_addon", shrink=0.5).fit(X, s, y).transform(X, s)
    assert (Z.groupby(s).mean() - Z.mean()).abs().max().max() > 0.1


def test_chorus_protected_variant_and_nested_candidates(sim):
    from chorus.harmonise import make_harmoniser as mk

    H = mk("chorus_protected")
    assert H.protect_scores == ["neutrophil", "lymphocyte"] and mk("chorus").protect_scores in (None, [], ())
    assert mk("chorus_protected", seed=3).seed == 3
    feats = list(sim.cohorts[0].features[:80])
    cfg = LosoConfig(harmonisers=["identity", "center_scale", "rank"], nested_candidates=["identity", "center_scale"], selectors=[{"name": "variance", "k": 20}], classifiers=["logreg"], n_boot=10, nested_select_harmoniser=True)
    res, _, _ = run_loso(sim, cfg, features=feats, verbose=False)
    assert set(res["harmoniser"]) == {"identity", "center_scale", "rank", "nested"}
    chosen = res.loc[res["harmoniser"] == "nested", "harmoniser_chosen"]
    assert len(chosen) and set(chosen) <= {"identity", "center_scale"}


def test_rerun_merge_places_a_new_harmoniser_in_loop_order(sim, tmp_path):
    import json

    from chorus.cli import _merge_rerun

    cfg = LosoConfig(harmonisers=["identity", "rank"], selectors=[{"name": "variance", "k": 20}], classifiers=["logreg"], n_boot=10, nested_select_harmoniser=True)
    feats = list(sim.cohorts[0].features[:80])
    res, preds, sel = run_loso(sim, cfg, features=feats, verbose=False)
    res.to_csv(tmp_path / "loso_results.csv", index=False)
    preds.to_csv(tmp_path / "loso_predictions.csv", index=False)
    (tmp_path / "selected_features.json").write_text(json.dumps({"|".join(k): v for k, v in sel.items()}))
    cfg2 = LosoConfig(harmonisers=["identity", "center_scale", "rank"], nested_candidates=["identity", "rank"], selectors=cfg.selectors, classifiers=["logreg"], n_boot=10, nested_select_harmoniser=True, outer_harmonisers=["center_scale"])
    res2, preds2, sel2 = run_loso(sim, cfg2, features=feats, verbose=False)
    r, p, s = _merge_rerun(tmp_path, ["center_scale"], res2, preds2, {"|".join(k): v for k, v in sel2.items()}, order=cfg2.harmonisers + ["nested"])
    assert len(r) == len(res) + len(res2)
    for _, g in r.groupby("held_out", sort=False):
        assert list(g["harmoniser"]) == ["identity", "center_scale", "rank", "nested"]


def test_forced_features_reach_the_classifier_for_every_selector(sim):
    # a forced feature that the selector would not pick (pure noise) is appended to the selected
    # features in every fold and for every selector, and is dropped when absent from the fold
    rng = np.random.default_rng(0)
    from chorus.io import Cohort, CohortSet

    cohorts = []
    for c in sim.cohorts:
        X = c.X.iloc[:, :60].copy()
        X["FORCED"] = rng.normal(0, 1e-3, len(X))  # near zero variance: never in the variance top k
        cohorts.append(Cohort(c.study_id, c.layer, X, c.meta, c.feature_space, False, dict(c.source)))
    cs = CohortSet(cohorts)
    feats = list(sim.cohorts[0].features[:60]) + ["FORCED"]
    sels = [{"name": "variance", "k": 10}, {"name": "invariant", "k": 10}]
    base = LosoConfig(harmonisers=["identity", "center_scale"], selectors=sels, classifiers=["logreg"], n_boot=10, nested_select_harmoniser=True)
    r0, _, s0 = run_loso(cs, base, features=feats, verbose=False)
    assert not any("FORCED" in v for v in s0.values())
    cfg = dataclasses.replace(base, forced_features=["FORCED", "NOT_IN_DATA"])
    r1, _, s1 = run_loso(cs, cfg, features=feats, verbose=False)
    assert all(v[-1] == "FORCED" and "NOT_IN_DATA" not in v for v in s1.values())
    assert (r1["n_forced_added"] == 1).all()
    assert (r1["n_features"] == r0["n_features"] + 1).all()


def test_layer_adapters_map_each_layer_to_its_shared_space():
    from chorus.adapters import ADAPTERS, canonical_uniprot, to_feature_space, uniprot_key

    assert {"bulk_rna", "olink", "somascan", "ms_proteomics"} <= set(ADAPTERS)
    assert ADAPTERS["olink"].feature_space == ADAPTERS["somascan"].feature_space == ADAPTERS["ms_proteomics"].feature_space == "uniprot"
    assert canonical_uniprot("P05231-2") == "P05231" and canonical_uniprot("CON__P02768") is None
    assert uniprot_key(["P29460", "P29459", "P29459"]) == "P29459_P29460"
    # bulk RNA: version suffixes removed, coinciding genes summed
    raw = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["s1", "s2"], columns=["ENSG00000000003.14", "ENSG00000000003.15", "ENSG00000000005.6"])
    a = to_feature_space("bulk_rna", raw)
    assert list(a.X.columns) == ["ENSG00000000003", "ENSG00000000005"] and a.X.loc["s2", "ENSG00000000003"] == 9
    # Olink long export: IL6 on two panels averaged, QC failures masked, a complex keyed by both accessions
    long = pd.DataFrame(
        {
            "SampleID": ["s1", "s1", "s1", "s2", "s2", "s2", "s1", "s2"],
            "OlinkID": ["OID1", "OID2", "OID3", "OID1", "OID2", "OID3", "OID4", "OID4"],
            "UniProt": ["P05231", "P05231", "P29460;P29459", "P05231", "P05231", "P29460;P29459", "", ""],
            "NPX": [1.0, 3.0, 5.0, 2.0, 4.0, 6.0, 9.0, 9.0],
            "QC": ["PASS", "PASS", "PASS", "PASS", "FAIL", "PASS", "PASS", "PASS"],
        }
    )
    o = to_feature_space("olink", long, qc_col="QC")
    assert sorted(o.X.columns) == ["P05231", "P29459_P29460"]
    assert o.X.loc["s1", "P05231"] == 2.0 and o.X.loc["s2", "P05231"] == 2.0  # s2 OID2 failed QC
    assert o.features.set_index("feature").loc["P05231", "n_sources"] == 2 and o.notes["assays_without_uniprot"] == "1"
    # SomaScan: linear RFU to log2, non-human aptamer dropped, two aptamers of one protein averaged
    rfu = pd.DataFrame({"10000-1": [4.0, 16.0], "10001-2": [16.0, 64.0], "10002-3": [256.0, 256.0]}, index=["s1", "s2"])
    ann = pd.DataFrame({"SeqId": ["10000-1", "10001-2", "10002-3"], "UniProt": ["P01024", "P01024", "P00000"], "Organism": ["Human", "Human", "Mouse"]})
    sm = to_feature_space("somascan", rfu, annotation=ann, keep=ann.set_index("SeqId")["Organism"] == "Human")
    assert list(sm.X.columns) == ["P01024"] and np.allclose(sm.X["P01024"], [3.0, 5.0])
    # MS protein groups: contaminants and decoys dropped, isoform suffix removed, zeros missing, log2
    pg = pd.DataFrame(
        {
            "Protein IDs": ["P02768-2;P02768", "CON__P02769", "REV__P12345", "sp|P01024|CO3_HUMAN", "P02768"],
            "Intensity a": [1024.0, 5.0, 5.0, 0.0, 4096.0],
            "Intensity b": [256.0, 5.0, 5.0, 64.0, 1024.0],
        }
    )
    ms = to_feature_space("ms_proteomics", pg, id_col="Protein IDs", sample_cols=["Intensity a", "Intensity b"])
    assert sorted(ms.X.columns) == ["P01024", "P02768"] and list(ms.X.index) == ["Intensity a", "Intensity b"]
    assert np.isnan(ms.X.loc["Intensity a", "P01024"]) and ms.X.loc["Intensity b", "P01024"] == 6.0
    assert ms.X.loc["Intensity a", "P02768"] == 11.0  # mean of log2 1024 and 4096
    assert ms.notes["groups_contaminant_or_decoy"] == "2"
    # DIA-NN sorts a group's accessions alphabetically; a preferred (reviewed) accession wins
    dia = pd.DataFrame({"Protein.Group": ["A0A087WWT3;P02768", "A0A0B4J2F0"], "r1": [8.0, 4.0]})
    d = to_feature_space("ms_proteomics", dia, id_col="Protein.Group", sample_cols=["r1"], preferred_accessions={"P02768"})
    assert sorted(d.X.columns) == ["A0A0B4J2F0", "P02768"] and d.notes["groups_keyed_outside_preferred"] == "1"
    assert sorted(to_feature_space("ms_proteomics", dia, id_col="Protein.Group", sample_cols=["r1"]).X.columns) == ["A0A087WWT3", "A0A0B4J2F0"]
    with pytest.raises(KeyError):
        to_feature_space("unknown_layer", raw)


def test_proteomics_write_cohort_contract_and_shared_protein_table(tmp_path, monkeypatch):
    import sys
    from pathlib import Path

    from chorus.adapters import to_feature_space

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "proteomics"))
    common = _load_script("proteomics/common")
    build_all = _load_script("proteomics/build_all")
    for mod in (common, build_all):
        monkeypatch.setattr(mod, "DATA", tmp_path / "proteomics")
        monkeypatch.setattr(mod, "META", tmp_path / "meta")

    def cohort(study, proteins, phases, status):
        long = pd.DataFrame(
            [{"SampleID": f"{study}_{i}", "OlinkID": f"O{p}", "UniProt": p, "NPX": float(i)} for i in range(len(phases)) for p in proteins]
        )
        a = to_feature_space("olink", long)
        meta = pd.DataFrame(
            {"subject_id": [f"{study}_s{i // 2}" for i in range(len(phases))], "lc_status": status, "case_definition": "study_specific", "sample_matrix": "plasma", "sample_phase": phases},
            index=a.X.index,
        )
        return a, meta

    a, m = cohort("A", ["P05231", "P01024", "P02768"], ["acute", "post_acute", "acute", "post_acute"], ["LC", "LC", "COVID_recovered", "COVID_recovered"])
    row = common.write_cohort("A", a, m)
    assert (row["n_LC"], row["n_COVID_recovered"], row["n_acute"], row["n_proteins"]) == (1, 1, 2, 3)
    written = pd.read_csv(tmp_path / "proteomics" / "A" / "A_matrix.csv", index_col=0)
    assert list(written.columns) == list(a.X.index) and written.index.name == "uniprot"
    meta = pd.read_csv(tmp_path / "meta" / "A_meta.csv")
    assert list(meta.columns[:7]) == ["sample_id", "subject_id", "study_id", "lc_status", "case_definition", "sample_matrix", "omics_layer"]
    assert (meta["omics_layer"] == "olink").all()
    b, mb = cohort("B", ["P05231", "P02768", "Q99999"], ["post_acute", "not_applicable"], ["LC", "never_infected"])
    common.write_cohort("B", b, mb)
    text = build_all.tables({"A": "no"})
    assert "| A | olink | 4 (2) | 1 (1) | 1 (1) | 0 | 2 | 0 |" in text
    assert "| A | 3 | 2 |" in text and "| B | 2 | 3 |" in text  # shared proteins, diagonal = cohort size
    # a vocabulary violation and a matrix / metadata mismatch are refused
    with pytest.raises(ValueError):
        common.write_cohort("C", b, mb.assign(lc_status="long_covid"))
    with pytest.raises(ValueError):
        common.write_cohort("C", b, mb.iloc[:1])


def _prot_cohort_files(tmp_path, study, X, meta):
    """Write a UniProt x sample matrix and its metadata as the proteomics builders do."""
    X.T.rename_axis("uniprot").to_csv(tmp_path / f"{study}_matrix.csv")
    meta.rename_axis("sample_id").assign(study_id=study, omics_layer="proteomics").to_csv(tmp_path / f"{study}_meta.csv")
    return {"study_id": study, "matrix": str(tmp_path / f"{study}_matrix.csv"), "meta": str(tmp_path / f"{study}_meta.csv"), "feature_space": "uniprot"}


def test_proteomics_preprocessing_phase_filter_timepoint_order_and_imputation(tmp_path):
    import yaml

    from chorus.cli import load_config

    rng = np.random.default_rng(0)
    base = {"case_definition": "study_specific", "sample_matrix": "plasma"}
    # A: subject a0 has an acute and a post acute LC sample; a1 has pre and post intervention samples (post first in file)
    ma = pd.DataFrame(
        [
            {"subject_id": "a0", "lc_status": "LC", "sample_phase": "acute", "timepoint": "T1"},
            {"subject_id": "a0", "lc_status": "LC", "sample_phase": "post_acute", "timepoint": "T3"},
            {"subject_id": "a1", "lc_status": "LC", "sample_phase": "post_acute", "timepoint": "post_intervention"},
            {"subject_id": "a1", "lc_status": "LC", "sample_phase": "post_acute", "timepoint": "pre_intervention"},
            {"subject_id": "a2", "lc_status": "COVID_recovered", "sample_phase": "post_acute", "timepoint": "single_visit"},
            {"subject_id": "a3", "lc_status": "never_infected", "sample_phase": "not_applicable", "timepoint": "single_visit"},
        ],
        index=[f"A{i}" for i in range(6)],
    ).assign(**base)
    Xa = pd.DataFrame(rng.normal(size=(6, 3)), index=ma.index, columns=["P1", "P2", "P3"])
    Xa.loc["A4", "P1"] = np.nan  # one missing value in a kept sample
    Xa.loc[:, "P3"] = np.nan  # P3 has no value in A: not shared
    mb = pd.DataFrame({"subject_id": ["b0", "b1"], "lc_status": ["LC", "COVID_recovered"], "sample_phase": "post_acute", "timepoint": "single_visit"}, index=["B0", "B1"]).assign(**base)
    Xb = pd.DataFrame([[1.0, 2.0, 3.0], [5.0, np.nan, 4.0]], index=mb.index, columns=["P1", "P2", "P4"])
    cfg = {
        "layer": "proteomics",
        "cohorts": [_prot_cohort_files(tmp_path, "A", Xa, ma), _prot_cohort_files(tmp_path, "B", Xb, mb)],
        "preprocessing": {
            "matrices": ["plasma"],
            "drop_sample_phases": ["acute"],
            "keep_lc_status": ["LC", "COVID_recovered"],
            "one_sample_per_subject": True,
            "acute_timepoints": [],
            "timepoint_order": ["T1", "T2", "T3", "pre_intervention", "post_intervention"],
            "impute": "study_min",
        },
        "benchmark": {},
    }
    cfg["preprocessing"].pop("matrices")  # select_tissue keeps cohorts with at least 6 samples of the matrix
    path = tmp_path / "c.yaml"
    path.write_text(yaml.safe_dump(cfg))
    cs, _, out = load_config(str(path))
    # a0: acute sample dropped; a1: pre_intervention kept although later in file; a3 (never_infected) dropped
    assert set(cs["A"].meta.index) == {"A1", "A3", "A4"}
    assert out["filtered_features"] == ["P1", "P2"]  # P3 has no value in A, P4 is absent from A
    a = cs["A"].X
    assert a.loc["A4", "P1"] == pytest.approx(Xa.loc[["A1", "A3"], "P1"].min())  # study minimum of the kept samples
    assert cs["B"].X.loc["B1", "P2"] == 2.0 and cs["B"].source["n_values_imputed"] == 1
    assert out["observed_fraction"].loc["P1", "A"] == pytest.approx(2 / 3)
    # union keeps every feature with a value somewhere; a feature absent from a cohort stays missing there
    cfg["preprocessing"].update(features="union", impute="study_median")
    path.write_text(yaml.safe_dump(cfg))
    cs, _, out = load_config(str(path))
    assert out["filtered_features"] == ["P1", "P2", "P4"]
    assert cs["A"].X["P4"].isna().all() and cs["B"].X["P4"].notna().all()


def test_relabel_acute_keeps_the_earliest_acute_draw_and_the_original_status():
    from chorus.io import Cohort, CohortSet

    meta = pd.DataFrame(
        {
            "subject_id": ["s0", "s0", "s0", "s1", "s2"],
            "lc_status": ["LC", "LC", "LC", "never_infected", "unknown"],
            "sample_phase": ["acute", "acute", "post_acute", "not_applicable", "acute"],
            "timepoint": ["T2", "T1", "T3", "single_visit", "T1"],
            "case_definition": "study_specific",
            "sample_matrix": "plasma",
        },
        index=["x0", "x1", "x2", "x3", "x4"],
    )
    cs = CohortSet([Cohort("S", "proteomics", pd.DataFrame({"P1": range(5)}, index=meta.index, dtype=float), meta)])
    cs = cs.relabel_acute().filter_samples(keep_lc_status=["acute", "never_infected"]).one_sample_per_subject([], ["T1", "T2", "T3"])
    m = cs["S"].meta
    assert set(m.index) == {"x1", "x3", "x4"}
    assert m.loc["x1", "lc_status"] == "acute" and m.loc["x1", "lc_status_original"] == "LC"
    assert m.loc["x4", "lc_status_original"] == "unknown"
    assert cs.summary()["acute"].iloc[0] == 2


def test_concordance_uniprot_symbols_and_pairwise_shared_counts(tmp_path):
    conc = _load_script("concordance")
    hgnc = tmp_path / "hgnc.txt"
    pd.DataFrame({"symbol": ["IL6", "C3", "ALB"], "ensembl_gene_id": ["E1", "E2", "E3"], "uniprot_ids": ["P05231", "P01024", "P02768|Q00001"]}).to_csv(hgnc, sep="\t", index=False)
    sym = conc.feature_symbols(pd.Index(["P05231", "P01024_P02768", "Q00001", "Z99999"]), str(hgnc), "uniprot")
    assert list(sym) == ["IL6", "C3_ALB", "ALB", None]
    assert list(conc.feature_symbols(pd.Index(["E2"]), str(hgnc))) == ["C3"]
    eff = pd.DataFrame({"a": [0.1, 0.2, np.nan, 0.4], "b": [0.3, np.nan, 0.1, 0.2], "c": [np.nan, np.nan, 1.0, 2.0]})
    n = conc.shared_counts(eff)
    assert n.loc["a", "a"] == 3 and n.loc["a", "b"] == 2 and n.loc["a", "c"] == 1 and n.loc["b", "c"] == 2


def test_curation_helper_drops_samples_and_empty_features_and_schema_vocabulary(monkeypatch):
    import sys
    from pathlib import Path

    from chorus.adapters import Adapted
    from chorus.schema import CONTROL_TYPES, LAYERS

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "proteomics"))
    common = _load_script("proteomics/common")
    X = pd.DataFrame({"P1": [1.0, 2.0, 3.0], "P2": [np.nan, np.nan, 4.0]}, index=["r1", "r2", "r3"])
    a = Adapted(X, pd.DataFrame({"feature": ["P1", "P2"], "source_ids": ["p1", "p2"], "n_sources": 1}), "uniprot", "ms_proteomics")
    b = common.drop_samples(a, pd.Series([True, True, False], index=X.index))
    assert list(b.X.index) == ["r1", "r2"] and list(b.X.columns) == ["P1"] and list(b.features["feature"]) == ["P1"]
    assert "convalescent_unverified" in CONTROL_TYPES and "proteomics" in LAYERS


def _paired_views(n=60, seed=0):
    rng = np.random.default_rng(seed)
    y = pd.Series(np.r_[np.ones(n // 2), np.zeros(n - n // 2)], index=[f"s{i}" for i in range(n)])
    a = pd.DataFrame(rng.normal(size=(n, 30)), index=y.index, columns=[f"P{i}" for i in range(30)])
    a.iloc[:, :3] += 1.5 * y.values[:, None]
    b = pd.DataFrame(rng.normal(size=(n, 20)), index=y.index, columns=[f"P{i}" for i in range(20)])
    b.iloc[:, :3] += 1.0 * y.values[:, None]
    b.iloc[0, 5] = np.nan  # a missing value, filled inside the fold
    return {"a": a, "b": b}, y


def test_paired_cv_same_folds_for_every_model_and_agreement_tables():
    from chorus.eval.paired import PairedCV, effect_correlation, probability_correlation, summarise_folds

    views, y = _paired_views()
    cv = PairedCV(views, y, impute={"b": "min"}, k=5, n_repeats=2)
    F, P = cv.run(models=["layer:a", "layer:b", "early_fusion", "late_fusion"], verbose=False)
    assert set(F["model"]) == {"layer:a", "layer:b", "early_fusion", "late_fusion"} and F["error"].isna().all()
    # every model is scored on the same test folds
    folds = {m: g.sort_values(["repeat", "fold", "sample"])[["repeat", "fold", "sample"]].values.tolist() for m, g in P.groupby("model")}
    assert all(v == folds["layer:a"] for v in folds.values())
    assert len(P[P["model"] == "layer:a"]) == 2 * len(y)  # each subject once per repeat
    assert (cv.folds().groupby("repeat")["sample"].nunique() == len(y)).all()
    s = summarise_folds(F)
    assert (s["n_folds"] == 10).all() and (s["n_failed_folds"] == 0).all()
    pc = probability_correlation(P)
    assert len(pc) == 2 * 6 and pc["n"].eq(len(y)).all()
    E, C = effect_correlation(views, y)
    assert C.loc[0, "n_shared"] == 20 and E.loc["P0", "a"] > 0.5
    # a sample missing from a layer is refused
    with pytest.raises(ValueError):
        PairedCV({"a": views["a"], "b": views["b"].iloc[1:]}, y)


def test_metabolomics_adapter_and_write_cohort_other_root(tmp_path, monkeypatch):
    import sys
    from pathlib import Path

    from chorus.adapters import to_feature_space

    raw = pd.DataFrame({" lactate ": [1.0, 2.0], "glucose": [3.0, np.nan], "lactate": [3.0, 4.0]}, index=["x", "y"])
    a = to_feature_space("metabolomics", raw, prefix="m:")
    assert list(a.X.columns) == ["m:glucose", "m:lactate"] and a.X.loc["x", "m:lactate"] == 2.0 and a.feature_space == "metabolite_name"
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "proteomics"))
    common = _load_script("proteomics/common")
    monkeypatch.setattr(common, "META", tmp_path / "meta")
    meta = pd.DataFrame({"subject_id": ["x", "y"], "lc_status": ["LC", "COVID_recovered"], "case_definition": "study_specific", "sample_matrix": "plasma"}, index=a.X.index)
    common.write_cohort("M", a, meta, root=tmp_path / "metabolomics")
    written = pd.read_csv(tmp_path / "metabolomics" / "M" / "M_matrix.csv", index_col=0)
    assert written.index.name == "metabolite_name" and list(written.columns) == ["x", "y"]


def test_paired_samples_agree_across_layers_and_keep_the_earliest_draw():
    pl = _load_script("paired_layers")
    meta = pd.DataFrame(
        {
            "subject_id": ["p1", "p1", "p1", "h1", "p2"],
            "lc_status": ["LC", "LC", "LC", "never_infected", "COVID_recovered"],
            "sample_phase": ["acute", "acute", "post_acute", "not_applicable", "post_acute"],
            "timepoint": ["T2", "T1", "T3", "single_visit", "T3"],
        },
        index=["p1-T2", "p1-T1", "p1-T3", "h1", "p2-T3"],
    )
    views = {"a": pd.DataFrame({"f": range(5)}, index=meta.index), "b": pd.DataFrame({"g": range(4)}, index=meta.index[[0, 1, 2, 4]])}
    acute = pl.paired_samples(views, {"a": meta, "b": meta}, "acute", ["T1", "T2", "T3"])
    assert list(acute.index) == ["p1-T1"]  # earliest acute draw; h1 is not in layer b
    lc = pl.paired_samples(views, {"a": meta, "b": meta}, "lc", ["T1", "T2", "T3"])
    assert set(lc.index) == {"p1-T3", "p2-T3"} and lc.loc["p1-T3", "y"] == 1.0
    other = meta.assign(lc_status=["LC", "LC", "COVID_recovered", "never_infected", "COVID_recovered"])
    with pytest.raises(ValueError):
        pl.paired_samples(views, {"a": meta, "b": other}, "lc", ["T1", "T2", "T3"])


def test_paired_cv_mofa_and_diablo_when_available():
    from chorus.eval.paired import PairedCV, diablo_available, mofa_available

    if not (mofa_available() and diablo_available()):
        pytest.skip("mofapy2 or R mixOmics not available")
    views, y = _paired_views(n=40)
    views = {k: v.add_prefix(f"{k}_") if k == "b" else v for k, v in views.items()}
    F, P = PairedCV(views, y, k=5, n_splits=3, n_repeats=1, mofa_factors_n=3, diablo_keepx=5).run(models=["layer:a", "mofa", "diablo"], verbose=False)
    assert F["error"].isna().all(), F["error"].dropna().tolist()
    assert set(F["model"]) == {"layer:a", "mofa", "diablo"}
    d = P[P["model"] == "diablo"]
    assert d["p"].between(0, 1).all() and len(d) == len(y)


def test_programme_mapping_scoring_and_minimum_members(tmp_path):
    from chorus.integrate.programmes import common_programmes, map_gene_sets, programme_correlation, programme_members, read_gmt, score_programmes

    (tmp_path / "s.gmt").write_text("SET_A\turl\t1\t2\t3\t4\t5\t6\nSET_B\turl\t1\t2\t99\n")
    sets = read_gmt(tmp_path / "s.gmt")
    assert sets["SET_B"] == ["1", "2", "99"]
    hgnc = pd.DataFrame(
        {
            "entrez_id": [str(i) for i in range(1, 7)],
            "symbol": [f"G{i}" for i in range(1, 7)],
            "ensembl_gene_id": [f"ENSG{i}" for i in range(1, 7)],
            "uniprot_ids": ["P1|Q1", "P2", "P3", "P4", "P5", None],
        }
    )
    T = map_gene_sets(sets, hgnc, "test")
    assert set(T.loc[T["entrez_id"] == "1", "uniprot"]) == {"P1", "Q1"}  # one row per accession
    assert T.loc[T["entrez_id"] == "99", "symbol"].isna().all()  # not in HGNC: kept, unmapped
    # a multi-protein key is a member when any accession is; it counts once
    mem = programme_members(T, ["P1_P9", "P2", "P3", "Q7"], "uniprot")
    assert sorted(mem["SET_A"]) == ["P1_P9", "P2", "P3"] and mem["SET_B"] == ["P1_P9", "P2"]
    ens = programme_members(T, [f"ENSG{i}" for i in range(1, 7)], "ensembl_gene")
    assert len(ens["SET_A"]) == 6

    rng = np.random.default_rng(0)
    feats = [f"ENSG{i}" for i in range(1, 7)]
    X = pd.DataFrame(rng.normal(size=(20, 6)), columns=feats, index=[f"s{i}" for i in range(20)])
    X.iloc[10:] += 5.0  # study offset
    X.loc[X.index[10:], "ENSG6"] = np.nan  # not measured in study b
    st = pd.Series(["a"] * 10 + ["b"] * 10, index=X.index)
    S, cov = score_programmes(X, st, ens, min_members=6)
    assert S.loc[st == "a", "SET_A"].notna().all() and S.loc[st == "b", "SET_A"].isna().all()  # 5 members measured in b
    assert np.allclose(S.loc[st == "a", "SET_A"].mean(), 0)  # within cohort standardisation
    c = cov.set_index(["study_id", "programme"])
    assert c.loc[("b", "SET_A"), "n_measured"] == 5 and not c.loc[("b", "SET_A"), "scored"]
    # the scores of a cohort do not depend on the other cohorts
    Sa, _ = score_programmes(X[st == "a"], st[st == "a"], ens, min_members=6)
    assert np.allclose(Sa["SET_A"], S.loc[st == "a", "SET_A"])
    S5, _ = score_programmes(X, st, ens, min_members=5)
    assert common_programmes(S5, st, ["a", "b"]) == ["SET_A"]
    C = programme_correlation(S5, S5)
    assert np.allclose(C["spearman"], 1.0)


def test_per_cohort_features_and_acute_relabel_by_status():
    from chorus.harmonise.transcriptomics import expressed_per_cohort, mask_unexpressed, to_continuous
    from chorus.io import Cohort, CohortSet

    meta = lambda ids, st: pd.DataFrame({"subject_id": ids, "lc_status": st, "case_definition": "study_specific", "sample_matrix": "whole_blood"}, index=ids)  # noqa: E731
    Xa = pd.DataFrame({"g1": [100, 100, 100, 100], "g2": [0, 0, 0, 1]}, index=["a0", "a1", "a2", "a3"])
    Xb = pd.DataFrame({"g1": [100, 100, 100, 100], "g2": [100, 100, 100, 100]}, index=["b0", "b1", "b2", "b3"])
    cs = CohortSet([Cohort("A", "bulk_rna", Xa, meta(list(Xa.index), ["LC", "COVID_recovered", "never_infected", "unknown"]), is_counts=True), Cohort("B", "bulk_rna", Xb, meta(list(Xb.index), ["LC"] * 4), is_counts=True)])
    ex = expressed_per_cohort(cs)
    assert ex["A"] == {"g1"} and ex["B"] == {"g1", "g2"}
    out = mask_unexpressed(to_continuous(cs, ["g1", "g2"]), ex)
    assert out["A"].X["g2"].isna().all() and out["B"].X["g2"].notna().all()
    # a cohort without sample_phase: the listed lc_status values are acute, the rest unchanged
    r = cs.relabel_acute(acute_lc_status={"A": ["LC", "COVID_recovered"]})
    assert list(r["A"].meta["lc_status"]) == ["acute", "acute", "never_infected", "unknown"]
    assert list(r["B"].meta["lc_status"]) == ["LC"] * 4 and list(r["A"].meta["lc_status_original"])[:2] == ["LC", "COVID_recovered"]


def test_mosaic_schemes_keep_sensitivity_cohorts_out_of_training():
    mp = _load_script("mosaic_programmes")
    envs = {"rna": ["rna:a", "rna:b"], "protein": ["protein:p", "protein:q", "protein:r"], "sensitivity": ["protein:p_ms"]}
    lc = mp.schemes("lc", envs)
    assert all("protein:p_ms" not in train for _, train, _ in lc)
    assert [t for s, _, t in lc if s == "rna_to_protein"] == ["protein:p", "protein:q", "protein:r", "protein:p_ms"]
    loo = [(train, t) for s, train, t in lc if s == "leave_one_cohort_out"]
    assert len(loo) == 5 and all(t not in train and len(train) == 4 for train, t in loo)
    acute = mp.schemes("acute", {"rna": ["rna:a"], "protein": ["protein:p", "protein:q", "protein:r"], "sensitivity": []})
    assert acute[0] == ("protein_to_rna", ["protein:p", "protein:q", "protein:r"], "rna:a")
    assert acute[1] == ("rna_plus_other_protein", ["rna:a", "protein:q", "protein:r"], "protein:p")


def test_metabolite_ids_priority_cache_and_adapter(tmp_path):
    from chorus.adapters import to_feature_space
    from chorus.metabolites import match_names, metabolite_ids, pairwise_shared, strip_prefix

    assert strip_prefix("GCMS::Glutamine") == "Glutamine" and strip_prefix("PC 36:4") == "PC 36:4"
    refmet = pd.DataFrame(
        {"refmet_id": ["RM1", "RM2", "RM3"], "hmdb_id": ["HMDB0000641", None, None], "inchi_key": ["K1", "K2", None], "kegg_id": [None, None, None]},
        index=pd.Index(["Glutamine", "PC 36:4", "Rare"], name="refmet_name"),
    )
    calls = []

    def fake(batch):
        calls.extend(batch)
        lut = {"L-Glutamine": "Glutamine", "Gln": "Glutamine", "PC aa C36:4": "PC 36:4", "Rare one": "Rare"}
        return {q: lut.get(q, "") for q in batch}

    cache = tmp_path / "cache.csv"
    queries = ["L-Glutamine", "Gln", "PC aa C36:4", "Rare one", "Nothing"]
    m = match_names(queries, cache, fetch=fake, chunk=2, pause=0)
    assert m["Gln"] == "Glutamine" and m["Nothing"] == ""
    m2 = match_names(queries, cache, fetch=fake, chunk=2, pause=0)  # second call answered from the cache
    assert len(calls) == 5 and m2.equals(m)
    names = ["LC::L-Glutamine", "GC::Gln", "PC aa C36:4", "Rare one", "Nothing", "104.07_11.07", "Deposited only"]
    dep = pd.DataFrame({"hmdb_id": [None, "HMDB0000123", "HMDB0000999"], "inchi_key": [None, None, None]}, index=["PC aa C36:4", "Deposited only", "Rare one"])
    T = metabolite_ids(names, refmet, pd.concat([m, pd.Series({"Deposited only": ""})]), dep).set_index("name")
    assert T.loc["LC::L-Glutamine", "shared_id"] == "HMDB:HMDB0000641" and T.loc["GC::Gln", "shared_id"] == "HMDB:HMDB0000641"
    assert T.loc["PC aa C36:4", "shared_id"] == "RefMet:PC 36:4"  # RefMet has no HMDB, no deposited HMDB for it
    assert T.loc["Rare one", "id_source"] == "hmdb_deposited"  # RefMet matched without HMDB, author gave one
    assert T.loc["Deposited only", "shared_id"] == "HMDB:HMDB0000123" and T.loc["Nothing", "id_source"] == "unmatched"
    assert T.loc["104.07_11.07", "id_source"] == "unannotated"
    # adapter: unmapped names dropped; two names of one identifier -> the most complete is kept
    X = pd.DataFrame({"LC::L-Glutamine": [1.0, np.nan, 3.0], "GC::Gln": [2.0, 2.0, 2.0], "PC aa C36:4": [4.0, 5.0, 6.0], "Nothing": [1.0, 1.0, 1.0]}, index=["s1", "s2", "s3"])
    ad = to_feature_space("metabolomics_id", X, ids=T["shared_id"])
    assert list(ad.X.columns) == ["HMDB:HMDB0000641", "RefMet:PC 36:4"]
    assert ad.X["HMDB:HMDB0000641"].tolist() == [2.0, 2.0, 2.0]
    f = ad.features.set_index("feature")
    assert f.loc["HMDB:HMDB0000641", "kept_source"] == "GC::Gln" and f.loc["HMDB:HMDB0000641", "n_sources"] == 2
    assert ad.notes["names_unmapped"] == "1" and ad.notes["names_dropped_as_duplicate"] == "1"
    S = pairwise_shared({"a": ["x", "y"], "b": ["y", "z"], "c": ["y"]})
    assert S.loc["a", "b"] == 1 and S.loc["a", "a"] == 2 and S.loc["b", "c"] == 1


def test_read_mwtab_metabolite_block(tmp_path):
    bc = _load_script_path("metabolomics/build_cohorts")
    txt = "#MS_METABOLITE_DATA\nMS_METABOLITE_DATA:UNITS\tArea\nMS_METABOLITE_DATA_START\nSamples\tC1\tL1\tQC1\nFactors\ta\tb\tc\nGlutamine\t10\t20\t15\nCAR 10:0 \t1\t\t2\nMS_METABOLITE_DATA_END\n"
    p = tmp_path / "x.mwtab.txt"
    p.write_text(txt)
    D = bc.read_mwtab_data(p)
    assert list(D.index) == ["C1", "L1", "QC1"] and list(D.columns) == ["Glutamine", "CAR 10:0"]
    assert D.loc["L1", "Glutamine"] == 20 and np.isnan(D.loc["L1", "CAR 10:0"])


def test_summary_keeps_configurations_that_failed_in_every_fold():
    from chorus.eval import summarise

    r = pd.DataFrame([
        {"held_out": "A", "harmoniser": "identity", "selector": "variance", "classifier": "logreg", "auroc": 0.6, "auprc": 0.5, "brier": 0.2, "error": None},
        {"held_out": "B", "harmoniser": "identity", "selector": "variance", "classifier": "logreg", "auroc": 0.8, "auprc": 0.7, "brier": 0.1, "error": None},
        {"held_out": "A", "harmoniser": "chorus", "selector": "variance", "classifier": None, "auroc": np.nan, "auprc": np.nan, "brier": np.nan, "error": "only 4 null genes passed the screen"},
        {"held_out": "B", "harmoniser": "chorus", "selector": "variance", "classifier": None, "auroc": np.nan, "auprc": np.nan, "brier": np.nan, "error": "only 6 null genes passed the screen"},
    ])
    s = summarise(r)
    assert s.loc[("identity", "variance", "logreg"), "mean"] == pytest.approx(0.7)
    assert s.loc[("identity", "variance", "logreg"), "n_errors"] == 0
    row = s.loc[("chorus", "variance", "not_fitted")]
    assert row["n_studies"] == 0 and row["n_errors"] == 2 and np.isnan(row["mean"])
    assert "only 4" in row["error"] and "only 6" in row["error"]
    assert list(s.index)[-1] == ("chorus", "variance", "not_fitted")  # failed configurations sort last


def test_salmon_prefix_download_keeps_equal_complete_mates(tmp_path, monkeypatch):
    import gzip

    q = _load_script_path("salmon_quant_ena")
    rng = np.random.default_rng(0)
    src = {}
    for mate in (1, 2):  # random reads, so that the compressed prefix is proportional to the reads it holds
        seqs = ["".join(rng.choice(list("ACGT"), 100)) for _ in range(20000)]
        src[mate] = gzip.compress("".join(f"@R.{i} x/{mate}\n{s}\n+\n{'F' * 100}\n" for i, s in enumerate(seqs, 1)).encode())

    def fake_download(url, dest, size, segments):  # leading bytes of the file, as the byte range download would give
        dest.write_bytes(src[1 if url.endswith("_1.fastq.gz") else 2][:size])

    monkeypatch.setattr(q, "download", fake_download)
    links = "h/RUN.fastq.gz;h/RUN_1.fastq.gz;h/RUN_2.fastq.gz"  # ENA lists an unpaired file next to the mates
    sizes = f"10;{len(src[1])};{len(src[2])}"
    assert q._mates(links.split(";"), sizes.split(";"))[0] == ["h/RUN_1.fastq.gz", "h/RUN_2.fastq.gz"]
    row = SimpleNamespace(run_accession="RUN", sample_alias="S1", read_count=20000, fastq_ftp=links, fastq_bytes=sizes, fastq_md5="a;b;c")
    files = q.fetch_prefix(row, tmp_path, 4, 5000, tmp_path / "prefix_reads.tsv")
    n = [sum(1 for _ in open(f)) for f in files]
    assert n[0] == n[1] == 4 * 5000
    assert open(files[0]).readlines()[-4].split()[0] == open(files[1]).readlines()[-4].split()[0]
    log = pd.read_csv(tmp_path / "prefix_reads.tsv", sep="\t")
    assert log.loc[0, "read_pairs_used"] == 5000 and log.loc[0, "read_pairs_in_run"] == 20000


def test_salmon_quantify_only_runs(tmp_path, monkeypatch):
    q = _load_script_path("salmon_quant_ena")
    runs = pd.DataFrame({"run_accession": ["A", "B", "C"], "sample_alias": ["sa", "sb", "sc"]})
    monkeypatch.setattr(q, "filereport", lambda project: runs)
    fetched = []

    def fake_fetch(row, tmp, segments):  # record the run and fail it, so that salmon is never called
        fetched.append(row.run_accession)
        raise RuntimeError("not downloaded")

    monkeypatch.setattr(q, "fetch", fake_fetch)
    q.quantify("P", tmp_path / "idx", tmp_path / "out", 1, 1, only_runs=["C", "A"])
    assert fetched == ["A", "C"]
    with pytest.raises(ValueError, match="runs not in P"):
        q.quantify("P", tmp_path / "idx", tmp_path / "out", 1, 1, only_runs=["D"])


def test_series_matrix_samples_and_gse297499_column_join(tmp_path):
    import gzip

    b = _load_script_path("build_rna_cohorts_2026_11")
    text = ('!Series_title\t"x"\n!Sample_title\t"2057_V2_PBMC"\t"2014R_PBMC"\t"2023_CSF"\n'
            '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\n'
            '!Sample_characteristics_ch1\t"tissue: PBMC"\t"tissue: PBMC"\t"tissue: Cerebrospinal Fluid"\n'
            '!Sample_characteristics_ch1\t"treatment: NP-Long COVID"\t"treatment: Control"\t"treatment: Control"\n')
    p = tmp_path / "sm.txt.gz"
    p.write_bytes(gzip.compress(text.encode()))
    g = b.series_matrix_samples(p)
    assert list(g["gsm"]) == ["GSM1", "GSM2", "GSM3"] and list(g["ch_treatment"]) == ["NP-Long COVID", "Control", "Control"]
    cols = g["title"].str.replace("_", "", regex=False).str.replace("PBMC$", "_PBMC", regex=True)
    assert list(cols[:2]) == ["2057V2_PBMC", "2014R_PBMC"]


def test_keep_timepoints_numeric_and_text():
    from chorus.io import Cohort, CohortSet

    X = pd.DataFrame(np.ones((4, 2)), index=list("abcd"), columns=["f1", "f2"])
    meta = pd.DataFrame({"sample_id": list("abcd"), "subject_id": list("abcd"), "study_id": "S", "lc_status": "LC", "timepoint": [0.0, 1.0, np.nan, 0.0]})
    cs = CohortSet([Cohort("S", "bulk_rna", X, meta.set_index("sample_id", drop=False))])
    assert list(cs.keep_timepoints([0])["S"].X.index) == ["a", "d"]
    meta["timepoint"] = ["T1", "T2", "T1", "acute"]
    cs = CohortSet([Cohort("S", "bulk_rna", X, meta.set_index("sample_id", drop=False))])
    assert list(cs.keep_timepoints(["T1"])["S"].X.index) == ["a", "c"]


def test_paired_layers_early_contrast():
    pl = _load_script_path("paired_layers")
    meta = pd.DataFrame({"lc_status": ["LC", "COVID_recovered", "LC", "never_infected", "COVID_recovered"], "sample_phase": ["acute", "acute", "post_acute", "not_applicable", "post_acute"]})
    assert pl.contrast_label(meta, "early").tolist()[:2] == [1.0, 0.0] and pl.contrast_label(meta, "early")[2:].isna().all()
    lc = pl.contrast_label(meta, "lc")
    assert lc.isna()[:2].all() and lc[2] == 1.0 and lc[4] == 0.0


def test_concordance_vs_transfer_statistics():
    cv = _load_script_path("concordance_vs_transfer")
    rng = np.random.default_rng(0)
    rows = []
    for b in ("A", "B"):
        for c in range(4):
            for h in cv.LABEL_FREE[:4]:
                m = rng.normal()
                rows.append({"block": b, "held_out": f"{b}{c}", "harmoniser": h, "effect_spearman": m, "sign_agreement": m, "auroc": 0.5 + 0.1 * m})
    df = pd.DataFrame(rows)
    r = cv.correlation(df, "effect_spearman", n_boot=50, n_perm=200)
    assert r["spearman"] == pytest.approx(1.0) and r["n_cohorts"] == 8 and r["p_permutation"] < 0.01
    v = cv.variance_shares(df).set_index("term")
    assert set(v.index) == {"held_out_cohort", "harmoniser", "residual"} and (v["share"] >= 0).all()


def test_clustering_consensus_separates_groups():
    cr = _load_script_path("clustering_readiness")
    rng = np.random.default_rng(0)
    P = np.vstack([rng.normal(0, 0.1, (20, 3)), rng.normal(5, 0.1, (20, 3))])
    labels, pac = cr.consensus(P, 2, 20, 0.8, 0)
    from sklearn.metrics import adjusted_rand_score

    assert adjusted_rand_score(np.repeat([0, 1], 20), labels) == 1.0 and pac == 0.0


def test_acute_persistence_cohort_stats_and_days():
    ap = _load_script_path("acute_persistence")
    rng = np.random.default_rng(0)
    y = pd.Series(np.repeat([1, 0], 30))
    days = pd.Series(rng.uniform(90, 400, 60))
    score = pd.Series(1.0 * y + 0.01 * days + rng.normal(0, 1, 60))
    st = ap.cohort_stats(score, y, days, n_boot=50, seed=0)
    assert st["n_lc"] == 30 and st["g"] > 0 and st["g_lo"] < st["g"] < st["g_hi"] and st["days_coef"] == pytest.approx(0.01, abs=0.01)
    assert np.isnan(ap.cohort_stats(score, y, pd.Series(np.nan, index=y.index), n_boot=10, seed=0)["g_adjusted"])


def test_schema_accepts_pulmonary_sequelae_case_definition():
    meta = pd.DataFrame({"sample_id": ["a", "b"], "subject_id": ["a", "b"], "study_id": "S", "lc_status": ["LC", "COVID_recovered"],
                         "case_definition": ["pulmonary_sequelae", "not_applicable"], "sample_matrix": "whole_blood", "omics_layer": "bulk_rna"})
    assert validate_metadata(meta, ClinicalSchema(), strict=True) == []
    with pytest.raises(ValueError):
        validate_metadata(meta.assign(case_definition="lung"), ClinicalSchema(), strict=True)


def test_nanostring_qc_flags():
    b = _load_script_path("build_rna_cohorts_2026_11")
    pos = {k: [v * 100, v * 100] for k, v in b.NANOSTRING_POS_FM.items()}
    raw = pd.DataFrame({**pos, **{f"NEG_{c}": [10, 10] for c in "ABCDEFGH"}, "HK1": [1000, 1000], "HK2": [500, 500]}, index=["s1", "s2"])
    raw.loc["s2", "NEG_A"] = 400  # negatives above POS_E (50) in s2
    q = b.nanostring_qc(raw, ["HK1", "HK2"])
    assert not q.loc["s1", "flag"] and q.loc["s2", "flag"] and "POS_E" in q.loc["s2", "flag_reason"]
    assert q["pos_scaling_factor"].round(6).tolist() == [1.0, 1.0] and q["pos_linearity_r2"].min() > 0.99
