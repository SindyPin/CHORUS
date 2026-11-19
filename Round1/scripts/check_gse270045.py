"""GSE270045 label check: sex from RPS4Y1 and Hedges g (LC minus recovered) against Shahbaz 2025.

Writes the per gene effects, the comparison with the genes named in the paper and the top 50 genes
by |g| to OUT_DIR. See docs/gse270045_check.md.

  python scripts/check_gse270045.py OUT_DIR
"""
import sys
import numpy as np, pandas as pd
S = sys.argv[1]
raw = pd.read_csv("data/counts/GSE270045/GSE270045_LC_counts.tsv.gz", sep="\t", index_col=0)
meta = pd.read_csv("data/meta/GSE270045_meta.csv").set_index("sample_id")
cpm = raw / raw.sum() * 1e6
sex = pd.Series(np.where(cpm.loc["RPS4Y1"] > 20, "male", "female"), index=raw.columns)
keep = (cpm >= 1).mean(axis=1) >= 0.5
L = np.log2(cpm[keep] + 0.5)
y = meta.loc[L.columns, "lc_status"]

def hedges_g(L, a, b):
    xa, xb = L.loc[:, a], L.loc[:, b]
    na, nb = a.sum(), b.sum()
    sp = np.sqrt(((na - 1) * xa.var(axis=1) + (nb - 1) * xb.var(axis=1)) / (na + nb - 2))
    d = (xa.mean(axis=1) - xb.mean(axis=1)) / sp.replace(0, np.nan)
    return d * (1 - 3 / (4 * (na + nb) - 9))

lc, rec = (y == "LC").values, (y == "COVID_recovered").values
fem, mal = (sex == "female").values, (sex == "male").values
g_all = hedges_g(L, lc, rec)
g_f = hedges_g(L, lc & fem, rec & fem)
g_m = hedges_g(L, lc & mal, rec & mal)
print("genes tested", keep.sum(), "| groups: LC F/M", (lc & fem).sum(), (lc & mal).sum(), " R F/M", (rec & fem).sum(), (rec & mal).sum())

ref = {  # direction in Shahbaz 2025 main text (LC minus recovered), by sex stratum
 "F": {"up": ["ZNF469","BRINP2","FEZF2","HOXC12","SKOR2","ESR1","ESR2","ESRRG","AR","NR3C1","CD69","CD302","CD226","CD200R1","MEIS2"],
       "down": ["WDR62","MUC12","SYCE1","SYCE1L","CD3E","CD14","CD7","CD320","CD101","CD82","CD81","CD6","CD68","CXCR3","IL32","IL15RA","IL27","IL2RB","IL12RB1","IL11RA"]},
 "M": {"up": ["PALB2","MEIS2","SCAF8","CLEC4D","LYZ","TSNAX-DISC1","CD36","CD48","CD164","TGFBR1","TFRC","ARTN","OR10H1","EPSTI1","CASP5"],
       "down": ["MCC","DNAAF2","HLA-DRB5","NECTIN2","GMPR","CD7","CD82","CD320","CD14","CD79A","IL27RA","IL2RA","IL11RA","CXCR3"]}}
rows = []
for stratum, gser in [("F", g_f), ("M", g_m)]:
    for direction, genes in ref[stratum].items():
        for gname in genes:
            v = gser.get(gname, np.nan) if gname in gser.index else (hedges_g(np.log2(cpm.loc[[gname]] + 0.5), lc & (fem if stratum == "F" else mal), rec & (fem if stratum == "F" else mal)).iloc[0] if gname in cpm.index else np.nan)
            rows.append({"stratum": "LCF vs RF" if stratum == "F" else "LCM vs RM", "gene": gname, "paper": direction, "g": round(float(v), 2) if pd.notna(v) else np.nan,
                         "agrees": (np.sign(v) == (1 if direction == "up" else -1)) if pd.notna(v) else np.nan})
ref_df = pd.DataFrame(rows)
ref_df.to_csv(f"{S}/gse270045_paper_genes.csv", index=False)
for st, d in ref_df.groupby("stratum"):
    d2 = d.dropna(subset=["g"])
    print(f"{st}: {int(d2.agrees.sum())} of {len(d2)} paper genes agree in sign ({len(d) - len(d2)} absent)")
print(ref_df.to_string(index=False))

top = g_all.abs().sort_values(ascending=False).head(50).index
t50 = pd.DataFrame({"g_all": g_all[top].round(2), "g_female": g_f[top].round(2), "g_male": g_m[top].round(2)})
t50.to_csv(f"{S}/gse270045_top50.csv")
named = set(sum([v for s in ref.values() for v in s.values()], []))
print("\nTop 50 by |g| overall: up", (t50.g_all > 0).sum(), "down", (t50.g_all < 0).sum(), "| named in paper:", sorted(set(top) & named))
print("Top 50: sign agrees between female and male strata:", int((np.sign(t50.g_female) == np.sign(t50.g_male)).sum()), "of 50")
print(t50.head(50).to_string())
pd.DataFrame({"g_all": g_all, "g_female": g_f, "g_male": g_m}).to_csv(f"{S}/gse270045_g.csv")
