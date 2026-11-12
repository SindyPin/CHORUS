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