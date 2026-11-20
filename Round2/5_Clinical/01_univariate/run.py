"""Step 01: each clinical variable alone (Amendment 5, C).

For each task, unit and variable kept in that unit (variables.csv): the AUROC of the variable alone and the odds ratio
of LC with its 95 percent interval. Where the same term is estimated in at least three units, a random effects pooled
odds ratio (DerSimonian and Laird) and heterogeneity (I squared). The timing variable is run the same way and marked
as a sensitivity analysis.

  numeric       AUROC of the raw value (higher value scored as more likely LC); odds ratio per unit from a logistic
                regression (age per 10 years, BMI per 5 kg/m2, WHO scale per point, counts per condition, time per 30 days)
  binary        AUROC of the 0/1 value; odds ratio from the 2 by 2 table (0.5 added to every cell when one is zero)
  categorical   each category against all others, as a binary term "variable=category"

Complete cases per variable; no imputation here. Units with fewer than 2 per class among the recorded are skipped.

  python Round2/5_Clinical/01_univariate/run.py
"""

import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import clinical as Q  # noqa: E402

OUT = Path(__file__).resolve().parent / "output"
SCALE = {"age": 10.0, "bmi": 5.0, "who_severity": 1.0, "comorbidity_count": 1.0, "n_preexisting_conditions": 1.0,
         "days_since_infection": 30.0, "days_since_infection_at_sample": 30.0}


def auroc(x, y):
    from sklearn.metrics import roc_auc_score

    return roc_auc_score(y, x) if len(set(y)) == 2 and len(set(x)) > 1 else np.nan


def or_table(x, y):
    a, b = int(((x == 1) & (y == 1)).sum()), int(((x == 1) & (y == 0)).sum())
    c, d = int(((x == 0) & (y == 1)).sum()), int(((x == 0) & (y == 0)).sum())
    note = ""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        note = "0.5 added to every cell (a zero cell)"
    lor = np.log(a * d / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return lor, se, note


def or_logit(x, y):
    import statsmodels.api as sm

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = sm.Logit(y, sm.add_constant(x)).fit(disp=0)
    return float(f.params[1]), float(f.bse[1]), ""


def pooled(lor, se):
    """DerSimonian and Laird random effects on log odds ratios; I squared from Cochran's Q."""
    w = 1 / se ** 2
    fixed = (w * lor).sum() / w.sum()
    Qs = (w * (lor - fixed) ** 2).sum()
    k = len(lor)
    tau2 = max(0.0, (Qs - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
    wr = 1 / (se ** 2 + tau2)
    m = (wr * lor).sum() / wr.sum()
    s = np.sqrt(1 / wr.sum())
    i2 = max(0.0, (Qs - (k - 1)) / Qs) if Qs > 0 else 0.0
    return m, s, tau2, i2, Qs, stats.chi2.sf(Qs, k - 1)


def terms(R, v):
    """(term, x, unit label, analysis) for one variable of one unit."""
    t = Q.V[v]["type"]
    x = Q.typed(v, R[v])
    if t == "categorical":
        for c in sorted(x.dropna().unique()):
            yield f"{v}={c}", (x == c).astype(float).where(x.notna()), "category vs all others", "binary"
    elif t == "binary":
        yield v, x, Q.V[v]["unit"], "binary"
    else:
        sd = x.std() if Q.V[v]["unit"] == "per SD" else SCALE.get(v, 1.0)
        yield v, x / sd, Q.V[v]["unit"], "numeric"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    U = pd.read_csv(Q.HERE / "output/units.csv").fillna("")
    D = pd.read_csv(Q.HERE / "variables.csv")
    rows = []
    for u in U.itertuples():
        with open(Q.CACHE / f"{u.task}__{u.unit}.pkl", "rb") as f:
            R = pickle.load(f)["R"]
        vs = [(v, False) for v in Q.kept_vars(D, u.unit, u.task)]
        if u.timing_variable:
            vs.append((u.timing_variable, True))
        for v, sens in vs:
            for term, x, unit, kind in terms(R, v):
                ok = x.notna()
                xx, yy = x[ok].to_numpy(dtype=float), R.loc[ok, "y"].to_numpy(dtype=int)
                n1, n0 = int((yy == 1).sum()), int((yy == 0).sum())
                r = dict(task=u.task, unit=u.unit, variable=v, term=term, sensitivity_timing=sens, or_unit=unit, n_LC=n1, n_recovered=n0)
                if n1 < 2 or n0 < 2 or len(set(xx)) < 2:
                    rows.append({**r, "note": "not estimable: fewer than 2 per class or no variation among the recorded"})
                    continue
                a = auroc(xx, yy)
                try:
                    lor, se, note = or_table(xx, yy) if kind == "binary" else or_logit(xx, yy)
                except Exception as e:  # noqa: BLE001 (separation and the like: recorded)
                    lor, se, note = np.nan, np.nan, f"odds ratio not estimable: {type(e).__name__}"
                if kind == "binary":
                    r.update(prevalence_LC=round(xx[yy == 1].mean(), 3), prevalence_recovered=round(xx[yy == 0].mean(), 3))
                else:
                    raw = Q.typed(v, R[v])[ok].to_numpy(dtype=float)
                    r.update(median_LC=float(np.median(raw[yy == 1])), median_recovered=float(np.median(raw[yy == 0])))
                rows.append({**r, "AUROC": a, "log_OR": lor, "se_log_OR": se, "OR": np.exp(lor),
                             "OR_low": np.exp(lor - 1.96 * se), "OR_high": np.exp(lor + 1.96 * se), "note": note})
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "by_cohort.csv", index=False)
    P = []
    ok = T.dropna(subset=["log_OR", "se_log_OR"])
    ok = ok[np.isfinite(ok["se_log_OR"]) & (ok["se_log_OR"] > 0)]
    for (task, term, sens), g in ok.groupby(["task", "term", "sensitivity_timing"]):
        if len(g) < 3:
            continue
        m, s, tau2, i2, Qs, pq = pooled(g["log_OR"].to_numpy(), g["se_log_OR"].to_numpy())
        P.append(dict(task=task, term=term, sensitivity_timing=sens, units=len(g), unit_names="; ".join(g["unit"]), or_unit=g["or_unit"].iloc[0],
                      n_LC=int(g["n_LC"].sum()), n_recovered=int(g["n_recovered"].sum()), pooled_OR=np.exp(m),
                      pooled_OR_low=np.exp(m - 1.96 * s), pooled_OR_high=np.exp(m + 1.96 * s), tau2=tau2, I2=i2, Q=Qs, Q_p=pq,
                      AUROC_median=g["AUROC"].median(), AUROC_min=g["AUROC"].min(), AUROC_max=g["AUROC"].max(),
                      interval_excludes_1=bool(np.exp(m - 1.96 * s) > 1 or np.exp(m + 1.96 * s) < 1)))
    P = pd.DataFrame(P)
    P.to_csv(OUT / "pooled.csv", index=False)
    pd.set_option("display.width", 250)
    if len(P):
        print(P[["task", "term", "units", "pooled_OR", "pooled_OR_low", "pooled_OR_high", "I2", "AUROC_median"]].round(3).to_string(index=False))
    print(f"{len(T)} unit-level estimates; {len(P)} pooled terms")


if __name__ == "__main__":
    main()
