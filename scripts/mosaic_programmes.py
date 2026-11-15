"""Test 3, mosaic: RNA and protein cohorts bridged at the programme level (docs/mosaic.md).

Every cohort is scored per programme (MSigDB Hallmark and Reactome, chorus.integrate.programmes)
from its own layer's features. The environment of a cohort is its cohort and layer
(``rna:GSE224615``, ``protein:Su2022``). Per contrast, in --out/<contrast>/:

  cohorts.csv               per environment: layer, cohort, samples per class, features measured,
                            programmes scored, role (primary or sensitivity)
  coverage.csv              per environment and programme: members in the layer's feature space,
                            members measured, fraction, scored (at least --min-members)
  effects.csv, effects_var.csv   Hedges g (positive minus negative) of every programme in every
                            environment, and its variance
  concordance_<set>.csv, concordance_<set>_p.csv, n_shared_<set>.csv   Spearman correlation of the
                            programme effects between environments over the programmes scored in
                            both, for set = all, hallmark, reactome
  partition.csv             invariance partition of the programmes with the environments as
                            environments (InvarianceAnchored.partition, no cap on the null set) over
                            the programmes scored in every environment of the analysis (the others:
                            set not_partitioned, with n_env_scored): meta-analysis table and set,
                            for the primary environments and, when the contrast has one, for the
                            sensitivity environment in place of the cohort it repeats (``analysis``)
  loso.csv, loso_predictions.csv   cross layer leave one cohort out in programme space: per scheme,
                            training environments, test environment, harmoniser and classifier,
                            the programmes used (scored in every training and the test cohort),
                            metrics on the test cohort with a bootstrap AUROC interval
  bridge_cervia.csv         (lc only) per programme, correlation across the same sera of the
                            SomaScan and MS scores, with each layer's own measured members
                            (``own``) and with the members measured in both layers (``shared``)

The scores themselves (scores.csv.gz) are written next to these files and are not committed.
Programme membership: --out/genesets.csv (per programme, genes, mapped to Ensembl and UniProt).

  python scripts/mosaic_programmes.py --out results/mosaic_programmes
  python scripts/mosaic_programmes.py --out results/mosaic_programmes_v4 --contrasts lc_v4
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from concordance import shared_counts, spearman_matrix  # noqa: E402

from chorus.cli import load_config  # noqa: E402
from chorus.eval.metrics import binary_metrics, bootstrap_ci  # noqa: E402
from chorus.features import per_study_effects  # noqa: E402
from chorus.harmonise import InvarianceAnchored, make_harmoniser  # noqa: E402
from chorus.integrate.programmes import common_programmes, map_gene_sets, programme_correlation, programme_members, read_gmt, score_programmes  # noqa: E402
from chorus.models import make_classifier  # noqa: E402
from chorus.schema import binary_label  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
HGNC = ROOT / "data/ref/hgnc_complete_set.txt"
GENESETS = {  # collection -> GMT with Entrez ids (MSigDB 2026.2, docs/provenance.md)
    "hallmark": ROOT / "data/ref/msigdb/h.all.v2026.2.Hs.entrez.gmt",
    "reactome": ROOT / "data/ref/msigdb/c2.cp.reactome.v2026.2.Hs.entrez.gmt",
}
HARMONISERS = ["identity", "center_scale"]
CLASSIFIER = "logreg"

# layer -> config; sensitivity: environments reported but never used for training or in the primary
# partition, with the cohort they repeat (Cervia2024_MS is the same sera as Cervia2024_SomaScan)
CONTRASTS = {
    "acute": {
        "configs": {"rna": "configs/mosaic_acute_rna.yaml", "protein": "configs/mosaic_acute_protein.yaml"},
        "sensitivity": {},
    },
    "lc": {
        "configs": {"rna": "configs/mosaic_lc_rna.yaml", "protein": "configs/concordance_proteomics_all.yaml"},
        "sensitivity": {"protein:Cervia2024_MS": "protein:Cervia2024_SomaScan"},
    },
    # version 4 cohort set (16 November 2026): the lc contrast with PRJNA1149196_SA, PRJNA1149196_QLD, GSE297499,
    # Fineschi2026 and Chowdhury2026 added; run with --contrasts lc_v4
    "lc_v4": {
        "configs": {"rna": "configs/mosaic_lc_rna_v4.yaml", "protein": "configs/concordance_proteomics_all_v4.yaml"},
        "sensitivity": {"protein:Cervia2024_MS": "protein:Cervia2024_SomaScan"},
    },
}


def gene_set_table() -> pd.DataFrame:
    hgnc = pd.read_csv(HGNC, sep="\t", dtype=str, usecols=["entrez_id", "symbol", "ensembl_gene_id", "uniprot_ids"])
    return pd.concat([map_gene_sets(read_gmt(p), hgnc, c) for c, p in GENESETS.items()], ignore_index=True)


def gene_set_summary(T: pd.DataFrame) -> pd.DataFrame:
    g = T.groupby(["collection", "programme"], sort=False)
    return pd.DataFrame(
        {
            "n_genes": g["entrez_id"].nunique(),
            "n_in_hgnc": g["symbol"].nunique(),
            "n_ensembl": g["ensembl_gene_id"].nunique(),
            "n_genes_with_uniprot": T[T["uniprot"].notna()].groupby(["collection", "programme"], sort=False)["entrez_id"].nunique(),
            "n_uniprot": g["uniprot"].nunique(),
        }
    ).fillna(0).astype(int).reset_index()


def load_layer(config: str, layer: str, T: pd.DataFrame, min_members: int) -> dict:
    cs, lc, cfg = load_config(str(ROOT / config))
    X, M = cs.pooled(cfg["filtered_features"])
    y = binary_label(M, lc.positive, lc.negative)
    X, M, y = X[y.notna()], M[y.notna()], y[y.notna()]
    env = (layer + ":" + M["study_id"].astype(str)).rename("env")
    members = programme_members(T, X.columns, cs.cohorts[0].feature_space)
    S, cov = score_programmes(X, env, members, min_members)
    cov.insert(0, "layer", layer)
    return {"X": X, "M": M, "y": y, "env": env, "S": S, "coverage": cov.rename(columns={"study_id": "env"}), "members": members, "space": cs.cohorts[0].feature_space, "positive": lc.positive, "negative": lc.negative}


def schemes(contrast: str, envs: dict) -> list:
    """(scheme, training environments, test environment) of the cross layer evaluation."""
    rna, prot = envs["rna"], envs["protein"]
    out = []
    if contrast == "acute":
        out += [("protein_to_rna", prot, t) for t in rna]
        out += [("rna_plus_other_protein", rna + [p for p in prot if p != t], t) for t in prot]
    else:
        out += [("rna_to_protein", rna, t) for t in prot + envs["sensitivity"]]
        out += [("protein_to_rna", prot, t) for t in rna]
        both = rna + prot
        out += [("leave_one_cohort_out", [e for e in both if e != t], t) for t in both]
    return out


def evaluate(S, y, env, train, test, n_boot=1000, min_test_per_class=5, seed=0) -> tuple:
    rows, preds = [], []
    progs = common_programmes(S, env, train + [test])
    tr, te = env.isin(set(train)), env == test
    base = {"train": ";".join(train), "test": test, "n_train": int(tr.sum()), "n_test_pos": int((y[te] == 1).sum()), "n_test_neg": int((y[te] == 0).sum()), "n_programmes": len(progs)}
    for h in HARMONISERS:
        row = {**base, "harmoniser": h, "classifier": CLASSIFIER}
        if not progs:
            rows.append({**row, "skipped": "no programme scored in every cohort"})
            continue
        if min(row["n_test_pos"], row["n_test_neg"]) < min_test_per_class:
            rows.append({**row, "skipped": f"fewer than {min_test_per_class} test samples in a class"})
            continue
        H = make_harmoniser(h).fit(S.loc[tr, progs], env[tr])
        Ztr, Zte = H.transform(S.loc[tr, progs], env[tr]), H.transform(S.loc[te, progs], env[te])
        clf = make_classifier(CLASSIFIER, seed).fit(Ztr, y[tr])
        p = clf.predict_proba(Zte)[:, 1]
        est, lo, hi = bootstrap_ci(y[te].values, p, roc_auc_score, n_boot=n_boot, seed=seed)
        rows.append({**row, **binary_metrics(y[te].values, p), "auroc_lo": lo, "auroc_hi": hi, "skipped": ""})
        preds.append(pd.DataFrame({"sample": S.index[te], "test": test, "harmoniser": h, "classifier": CLASSIFIER, "y": y[te].values, "p": p}))
    return rows, preds


def run_contrast(contrast: str, out: Path, T: pd.DataFrame, min_members: int = 5, n_boot: int = 1000) -> dict:
    spec = CONTRASTS[contrast]
    d = out / contrast
    d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    L = {layer: load_layer(cfg, layer, T, min_members) for layer, cfg in spec["configs"].items()}
    S = pd.concat([L[l]["S"] for l in L])
    y = pd.concat([L[l]["y"] for l in L])
    env = pd.concat([L[l]["env"] for l in L])
    S.to_csv(d / "scores.csv.gz")
    sens = spec["sensitivity"]
    primary = [e for e in pd.unique(env) if e not in sens]
    order = primary + [e for e in pd.unique(env) if e in sens]
    envs = {"rna": [e for e in primary if e.startswith("rna:")], "protein": [e for e in primary if e.startswith("protein:")], "sensitivity": list(sens)}
    collection = T.drop_duplicates("programme").set_index("programme")["collection"]

    cov = pd.concat([L[l]["coverage"] for l in L], ignore_index=True)
    cov.insert(3, "collection", cov["programme"].map(collection))
    cov.to_csv(d / "coverage.csv", index=False)
    rows = []
    for e in order:
        layer = e.split(":")[0]
        m = env == e
        X = L[layer]["X"][L[layer]["env"] == e]
        c = cov[cov["env"] == e]
        rows.append(
            {
                "env": e, "layer": layer, "cohort": e.split(":", 1)[1], "role": "sensitivity" if e in sens else "primary",
                "n_positive": int((y[m] == 1).sum()), "n_negative": int((y[m] == 0).sum()),
                "n_features_measured": int(X.notna().any(axis=0).sum()),
                **{f"n_scored_{k}": int(c.loc[c["collection"] == k, "scored"].sum()) for k in GENESETS},
                **{f"n_programmes_{k}": int((c["collection"] == k).sum()) for k in GENESETS},
            }
        )
    pd.DataFrame(rows).to_csv(d / "cohorts.csv", index=False)

    eff = per_study_effects(S, y, env)
    G = pd.DataFrame({e: eff[e]["g"] for e in order if e in eff})
    V = pd.DataFrame({e: eff[e]["var"] for e in order if e in eff})
    G.insert(0, "collection", collection.reindex(G.index).values)
    V.insert(0, "collection", collection.reindex(V.index).values)
    G.rename_axis("programme").to_csv(d / "effects.csv")
    V.rename_axis("programme").to_csv(d / "effects_var.csv")
    g = G.drop(columns="collection")
    for name, sel in [("all", g.index), *[(k, g.index[G["collection"] == k]) for k in GENESETS]]:
        r, p = spearman_matrix(g.loc[sel])
        r.to_csv(d / f"concordance_{name}.csv")
        p.to_csv(d / f"concordance_{name}_p.csv")
        shared_counts(g.loc[sel]).to_csv(d / f"n_shared_{name}.csv")

    parts = []
    analyses = {"primary": primary}
    for s_env, repl in sens.items():
        analyses[f"sensitivity_{s_env.split(':', 1)[1]}"] = [s_env if e == repl else e for e in primary]
    for name, use in analyses.items():
        m = env.isin(set(use))
        n_scored = S[m].notna().groupby(env[m]).any().sum(axis=0)
        cols = common_programmes(S, env, use)  # the operator assumes a complete matrix
        P = InvarianceAnchored(n_controls=10**9)
        sets = P.partition(S.loc[m, cols], y[m], env[m])
        tab = P.table_.copy()
        tab["set"] = "none"
        for k in ("invariant", "null", "unstable"):
            tab.loc[sets[k], "set"] = k
        tab["max_abs_g"] = P.per_study_g_.abs().max(axis=1)
        rest = pd.DataFrame({"set": "not_partitioned"}, index=pd.Index([c for c in S.columns if c not in set(cols)]))
        tab = pd.concat([tab, rest])
        tab.insert(0, "analysis", name)
        tab.insert(1, "collection", collection.reindex(tab.index).values)
        tab["n_env_scored"] = n_scored.reindex(tab.index).astype(int).values
        tab["n_env"] = len(use)
        tab["environments"] = ";".join(use)
        parts.append(tab.rename_axis("programme").reset_index())
    pd.concat(parts, ignore_index=True).to_csv(d / "partition.csv", index=False)

    res, prs = [], []
    for scheme, train, test in schemes(contrast, envs):
        r, p = evaluate(S, y, env, train, test, n_boot=n_boot)
        res += [{"scheme": scheme, "test_role": "sensitivity" if test in sens else "primary", **x} for x in r]
        prs += [x.assign(scheme=scheme) for x in p]
    R = pd.DataFrame(res)
    R.to_csv(d / "loso.csv", index=False)
    if prs:
        pd.concat(prs, ignore_index=True)[["scheme", "test", "harmoniser", "classifier", "sample", "y", "p"]].to_csv(d / "loso_predictions.csv", index=False)
    info = {"contrast": contrast, "positive": L["rna"]["positive"], "negative": L["rna"]["negative"], "configs": spec["configs"], "environments": order, "sensitivity": sens, "min_members": min_members}

    if contrast in ("lc", "lc_v4") and {"protein:Cervia2024_SomaScan", "protein:Cervia2024_MS"} <= set(order):
        B = bridge_cervia(L["protein"], collection, min_members)
        B.to_csv(d / "bridge_cervia.csv", index=False)
        info["bridge_n_programmes"] = B.groupby("member_set")["programme"].size().to_dict()
    info["seconds"] = round(time.time() - t0, 1)
    print(R.drop(columns=["train"], errors="ignore").to_string(index=False), flush=True)
    return info


def bridge_cervia(P: dict, collection: pd.Series, min_members: int) -> pd.DataFrame:
    """Correlation per programme of the SomaScan and MS scores of the same sera."""
    a, b = "protein:Cervia2024_SomaScan", "protein:Cervia2024_MS"
    env, X = P["env"], P["X"]
    sid = lambda idx: [i.split("::", 1)[1] for i in idx]  # noqa: E731
    Xa, Xb = X[env == a], X[env == b]
    Xa.index, Xb.index = sid(Xa.index), sid(Xb.index)
    both = Xa.index.intersection(Xb.index)
    Xa, Xb = Xa.loc[both], Xb.loc[both]
    meas = lambda Z: set(Z.columns[Z.notna().any(axis=0).values & (Z.std(axis=0) > 0).values])  # noqa: E731
    ma, mb = meas(Xa), meas(Xb)
    out = []
    for member_set in ("own", "shared"):
        keep = (lambda f, own: f in own) if member_set == "own" else (lambda f, own: f in ma and f in mb)
        mem_a = {p: [f for f in fs if keep(f, ma)] for p, fs in P["members"].items()}
        mem_b = {p: [f for f in fs if keep(f, mb)] for p, fs in P["members"].items()}
        Sa, _ = score_programmes(Xa, pd.Series("a", index=both), mem_a, min_members)
        Sb, _ = score_programmes(Xb, pd.Series("b", index=both), mem_b, min_members)
        C = programme_correlation(Sa, Sb)
        C.insert(0, "member_set", member_set)
        C.insert(2, "collection", C["programme"].map(collection))
        C["n_members_somascan"] = C["programme"].map({p: len(v) for p, v in mem_a.items()})
        C["n_members_ms"] = C["programme"].map({p: len(v) for p, v in mem_b.items()})
        C["n_members_shared"] = C["programme"].map({p: len(set(mem_a[p]) & set(mem_b[p])) for p in mem_a})
        out.append(C)
    return pd.concat(out, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--contrasts", default="acute,lc", help="run order; the positive control first")
    ap.add_argument("--min-members", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=1000)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    T = gene_set_table()
    gene_set_summary(T).to_csv(out / "genesets.csv", index=False)
    info = {"genesets": {k: str(p.relative_to(ROOT)) for k, p in GENESETS.items()}, "harmonisers": HARMONISERS, "classifier": CLASSIFIER, "contrasts": {}}
    for c in a.contrasts.split(","):
        info["contrasts"][c] = run_contrast(c, out, T, a.min_members, a.n_boot)
        (out / "run_info.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
