"""Build the Captur2022 targeted MS proteomics cohort (Captur G, ..., Heywood WE, UK COVIDsortium
Investigators 2022 eBioMedicine 85:104293, "Plasma proteomic signature predicts who will get
persistent symptoms following SARS-CoV-2 infection", doi:10.1016/j.ebiom.2022.104293, PMID 36182629,
PMC9515404; PRIDE PXD036590, CC0).

Raw files (data/proteomics/Captur2022/raw):
  Plate_<k>_Injection<j>_TransitionResults_QuantQuali_YE_IS.csv   Skyline transition results of the
      targeted MRM assay, plates 1 to 4. Injection1 and Injection2 are NOT technical replicates:
      they are two different MRM methods (replicate names "Barts HCW Inflam 1" / "Inflam 2"), with
      almost disjoint peptide lists (only ALDOA, IGHA2 TPL, NRP1 ATH, C3 SSL and THBS1 TIV are in
      both). Every plasma sample was injected once with each method (suffix _A = Injection1, _B =
      Injection2). Each peptide has one quantifier and one qualifier transition (column
      "Quantitative"). The internal standard is a spiked yeast enolase peptide (GNPTVEVELTTEK,
      "YE_GP" in Injection1, "YeastEno_GNP" in Injection2); ALDOA ALQASALK is the only peptide
      with a stable isotope labelled (heavy Leu) partner.
  Supplementary_Data_1.xlsx   pathfindR KEGG enrichment output only (no protein quantities, no
      sample or clinical data); not used.

Quantity per peptide and sample: sum of the areas of the transitions Skyline flags Quantitative
(the quantifier transition), divided by the same quantity of the yeast enolase GNPTVEVELTTEK
internal standard measured in the same injection (the paper: "Integrated peak areas were expressed
as a ratio to internal standard"). For ALDOA the ratio is light / heavy ALQASALK. Missing (#N/A) or
zero areas are missing. Peptides go to the ms_proteomics adapter as rows keyed by the protein's
UniProt accession, with ``log2=True``; the adapter averages all peptides of one protein (both
injections) on the log2 scale, so each matrix value is the mean log2 peptide/IS ratio.
The authors' z scoring across analytes is not applied.

Protein labels are UniProt entry-name stems, gene symbols or FASTA headers. They are mapped to gene
symbols here (LABEL_GENE) and to UniProt with the HGNC complete set (data/ref/hgnc_complete_set.txt);
every peptide sequence was checked against the UniProt sequence of the mapped entry when the map
was written (CDKL5 GTSMCPTLQVR matches a CDKL5 isoform). Excluded: the yeast enolase standards, and
the APOE allele specific peptides (APOE_E2, APOE_E3, APOE_E4), whose signal depends on the APOE
genotype; APOE is quantified from AATVGSLAGQPLQER and LGPLVEQGR.

Samples: replicate names "BHC<participant>_<tube barcode> _A|_B". Calibration curves, pooled
plasma QCs and washes are dropped. The two injections of a sample are joined by plate and tube
barcode. For 10 barcodes the participant number differs between the injections (e.g. BHC141 in
Injection1 vs BHC142 in Injection2 for barcode 1208861964); the Injection1 number is used (in each
case the Injection2 list gives two consecutive barcodes the same participant), and both are kept in
``notes`` / ``bhc_injection2``. One plate 3 well is "BHC81_1217889158" in Injection1 and
"BHC77_1217889156" in Injection2 (barcode and participant both differ); it is dropped.

Labels: the deposited files carry no timepoint, infection status, symptom or demographic
information (the authors' clinical data are available only by application to the COVIDsortium
data access portal). So lc_status and sample_phase are ``unknown`` for every sample. The column
``n_samples_subject`` gives the number of samples of the participant in this dataset (uninfected
HCW contributed one baseline sample, infected HCW serial samples); it is not used as a label.

No vendor QC flags exist in the Skyline export, so no values are masked; no normalisation,
imputation or cross-plate correction is applied.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, raw_dir, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Captur2022"
HGNC = ROOT / "data" / "ref" / "hgnc_complete_set.txt"
PLATFORM = "Targeted LC-MRM (custom ~91 protein inflammation/neuroinflammation panel, two MRM methods; Skyline); peptide area / yeast enolase internal standard"
IS_LABELS = {"YE_GP", "YeastEno_GNP"}  # yeast enolase GNPTVEVELTTEK, the internal standard
EXCLUDE = {"YE AVD", "YE SIV", "YEAST ENOLASE AVD", "YEAST ENOLASE SIV", "APOE_E2", "APOE_E3", "APOE_E4"} | IS_LABELS
SAMPLE_RE = re.compile(r"^BHC(\d+)_(\d+)\s*_+([AB])$")
FILE_RE = re.compile(r"Plate_(\d)_Injection_?(\d)_")

# label stem (text before "_", "-" or " ") -> HGNC symbol
LABEL_GENE = {
    "A1AT": "SERPINA1", "A2AP": "SERPINF2", "A2MG": "A2M", "AACT": "SERPINA3", "ABCA7": "ABCA7",
    "ADIPOQ": "ADIPOQ", "ALDOA": "ALDOA", "GSTO1": "GSTO1", "APOA4": "APOA4", "APOE": "APOE",
    "APP": "APP", "C3": "C3", "TSP1": "THBS1", "B4GT1": "B4GALT1", "BCHE": "BCHE", "BST1": "BST1",
    "CAPN2": "CAPN2", "CCL4": "CCL4", "CD5L": "CD5L", "CDKL5": "CDKL5", "CO6A3": "COL6A3",
    "CPSM": "CPS1", "CP": "CP", "CRLD1": "CRISPLD1", "CRP": "CRP", "CSF1R": "CSF1R", "CYCS": "CYCS",
    "CYTC": "CST3", "CATB": "CTSB", "C4": "C4A", "DCXR": "DCXR", "DKK3": "DKK3", "DYH12": "DNAH12",
    "DYH1": "DNAH1", "FABP5": "FABP5", "FGA": "FGA", "FGF21": "FGF21", "FGF2": "FGF2", "GM2A": "GM2A",
    "SAP3": "GM2A", "GOLM1": "GOLM1", "GOT1": "GOT1", "AATC": "GOT1", "GRN": "GRN", "GUAD": "GDA",
    "HBE1": "HBE1", "HPX": "HPX", "HS71L": "HSPA1L", "HS90B": "HSP90AB1", "HSCB": "HSCB",
    "HSP7C": "HSPA8", "HSPA5": "HSPA5", "IC1": "SERPING1", "ICAM1": "ICAM1", "IGHA2": "IGHA2",
    "IGHE": "IGHE", "IL4": "IL4", "IL7": "IL7", "ITIH2": "ITIH2", "LAMP1": "LAMP1", "LMO7": "LMO7",
    "LYAM2": "SELE", "MASP2": "MASP2", "MIEN1": "MIEN1", "MMP3": "MMP3", "MOES": "MSN",
    "MYO10": "MYO10", "NCAM1": "NCAM1", "NCAM2": "NCAM2", "NEFH": "NEFH", "NEFL": "NEFL",
    "NEFM": "NEFM", "NPC2": "NPC2", "NRP1": "NRP1", "OLR1": "OLR1", "PEDF": "SERPINF1", "PGK1": "PGK1",
    "PKM": "PKM", "PLD3": "PLD3", "PPP3CB": "PPP3CB", "PRDX3": "PRDX3", "PRG4": "PRG4",
    "PTGDS": "PTGDS", "QSOX1": "QSOX1", "REG1B": "REG1B", "RSF1": "RSF1", "S100A9": "S100A9",
    "SAA1": "SAA1", "SCRG1": "SCRG1", "SIR2": "SIRT2", "SOD3": "SOD3", "SPP2": "SPP2",
    "TCHHL1": "TCHHL1", "TFEB": "TFEB", "TNNT3": "TNNT3", "TRAP1": "TRAP1", "TUBA4A": "TUBA4A",
    "UH1BL": "BLTP3B", "VCAM1": "VCAM1", "VCAM": "VCAM1", "HEXA": "HEXA", "CATH": "CTSH",
    "GALNS": "GALNS", "PEBP1": "PEBP1",
}
LABEL_ALIAS = {"AP GSTO1": "GSTO1", "AP_C3_SSL": "C3", "AP_TSP1_TIV": "TSP1", "Cathepsin B": "CATB", "CompC4-VTA": "C4"}


def label_stem(label: str) -> str:
    if label in LABEL_ALIAS:
        return LABEL_ALIAS[label]
    if label.startswith("sp|"):  # sp|P07858|CATB_HUMAN
        return label.split("|")[2].split("_")[0]
    return re.split(r"[_\- ]", label)[0]


def hgnc_uniprot() -> dict:
    h = pd.read_csv(HGNC, sep="\t", usecols=["symbol", "uniprot_ids"], dtype=str)
    return h.dropna().set_index("symbol")["uniprot_ids"].to_dict()


def read_injection(path: Path) -> pd.DataFrame:
    """Long table: one row per sample replicate and peptide (Protein label, Peptide, precursor),
    quantifier area summed over Quantitative transitions."""
    d = pd.read_csv(path, na_values=["#N/A"], low_memory=False)
    d = d[d["Quantitative"].astype(str) == "True"]
    d = d[d["Replicate"].str.match(SAMPLE_RE)]
    d["Area"] = pd.to_numeric(d["Area"], errors="coerce")
    g = d.groupby(["Replicate", "Protein", "Peptide", "Precursor Mz"], sort=False)["Area"].sum(min_count=1).reset_index()
    return g


def peptide_ratios(g: pd.DataFrame) -> pd.DataFrame:
    """Replicates x peptide rows: area / IS area; ALDOA light / heavy."""
    g = g.copy()
    g.loc[g["Area"] <= 0, "Area"] = np.nan
    is_area = g[g["Protein"].isin(IS_LABELS)].groupby("Replicate")["Area"].first()
    rows = []
    # ALDOA: two precursors of ALQASALK, the heavier one is the labelled standard
    al = g[g["Protein"] == "ALDOA_ALQ"]
    if len(al):
        mz = sorted(al["Precursor Mz"].unique())
        light = al[al["Precursor Mz"] == mz[0]].set_index("Replicate")["Area"]
        heavy = al[al["Precursor Mz"] == mz[-1]].set_index("Replicate")["Area"]
        rows.append(("ALDOA_ALQ", "ALQASALK", light / heavy.reindex(light.index)))
    for (prot, pep), s in g[~g["Protein"].isin(EXCLUDE | {"ALDOA_ALQ"})].groupby(["Protein", "Peptide"], sort=False):
        s = s.groupby("Replicate")["Area"].sum(min_count=1)
        rows.append((prot, pep, s / is_area.reindex(s.index)))
    out = pd.DataFrame({f"{p}::{q}": v for p, q, v in rows})
    return out, is_area


def build() -> list:
    raw = raw_dir(STUDY)
    up = hgnc_uniprot()
    blocks, samples = [], []
    for path in sorted(raw.glob("Plate_*_TransitionResults_QuantQuali_YE_IS.csv")):
        m = FILE_RE.search(path.name)
        plate, inj = int(m.group(1)), int(m.group(2))
        R, is_area = peptide_ratios(read_injection(path))
        for rep in R.index:
            mm = SAMPLE_RE.match(rep)
            samples.append({"plate": plate, "inj": inj, "replicate": rep, "bhc": int(mm.group(1)), "barcode": mm.group(2), "is_area": is_area.get(rep, np.nan)})
        R.index = pd.MultiIndex.from_arrays([[plate] * len(R), [inj] * len(R), list(R.index)], names=["plate", "inj", "replicate"])
        blocks.append(R)
    S = pd.DataFrame(samples)

    # join the two injections of a sample by plate and barcode
    w = S.pivot_table(index=["plate", "barcode"], columns="inj", values="bhc", aggfunc="first")
    unmatched = w[w.isna().any(axis=1)]
    w = w.dropna().astype(int)
    rep = S.set_index(["plate", "barcode", "inj"])["replicate"]
    isa = S.set_index(["plate", "barcode", "inj"])["is_area"]

    cols, meta_rows = {}, []
    by_inj = {i: pd.concat([b for b in blocks if b.index.get_level_values("inj")[0] == i]) for i in (1, 2)}
    for (plate, bc), r in w.iterrows():
        b1, b2 = int(r[1]), int(r[2])
        sid = f"BHC{b1}_{bc}"
        v1 = by_inj[1].loc[(plate, 1, rep[(plate, bc, 1)])]
        v2 = by_inj[2].loc[(plate, 2, rep[(plate, bc, 2)])]
        v = pd.concat([v1.dropna().rename(lambda c: f"I1::{c}"), v2.dropna().rename(lambda c: f"I2::{c}")])
        cols[sid] = v
        note = []
        if b1 != b2:
            note.append(f"participant number differs between injections (Injection1 BHC{b1}, Injection2 BHC{b2}); Injection1 used")
        meta_rows.append({
            "sample_id": sid, "subject_id": f"{STUDY}_BHC{b1}", "bhc_injection2": f"BHC{b2}", "tube_barcode": bc, "plate": plate,
            "replicate_injection1": rep[(plate, bc, 1)], "replicate_injection2": rep[(plate, bc, 2)],
            "is_area_injection1": isa[(plate, bc, 1)], "is_area_injection2": isa[(plate, bc, 2)], "notes": "; ".join(note),
        })
    P = pd.DataFrame(cols)  # peptide rows x samples
    P.index.name = "peptide_key"

    # peptide rows -> UniProt accession of the protein
    labels = P.index.str.split("::").str[1]
    stems = labels.map(label_stem)
    genes = stems.map(LABEL_GENE)
    if genes.isna().any():
        raise ValueError(f"unmapped labels: {sorted(set(labels[genes.isna()]))}")
    acc = genes.map(up)
    if acc.isna().any():
        raise ValueError(f"genes without UniProt in HGNC: {sorted(set(genes[acc.isna()]))}")
    if acc.str.contains(r"\|").any():
        raise ValueError(f"genes with several UniProt ids: {sorted(set(genes[acc.str.contains('[|]')]))}")
    table = P.reset_index(drop=False)
    table.insert(0, "protein_id", acc.values)
    table["protein_id"] = table["protein_id"] + " " + table["peptide_key"]  # leading accession first
    sample_cols = list(P.columns)
    adapted = to_feature_space("ms_proteomics", table, id_col="protein_id", sample_cols=sample_cols, log2=True)

    meta = pd.DataFrame(meta_rows).set_index("sample_id").loc[adapted.X.index]
    nsub = meta.groupby("subject_id").size()
    meta["n_samples_subject"] = meta["subject_id"].map(nsub).values
    meta["lc_status"] = "unknown"
    meta["case_definition"] = "study_specific"
    meta["sample_matrix"] = "plasma"
    meta["platform"] = PLATFORM
    meta["sample_phase"] = "unknown"
    meta["time_anchor"] = "unknown"
    meta["acute_severity"] = "unknown"
    meta["hospitalised"] = False  # "None of the SARS-CoV-2 infected individuals required hospitalization"
    meta["variant_era"] = "pre_alpha"  # first UK wave, recruitment from 23 March 2020
    meta["vaccinated_before_infection"] = False  # "before the availability of vaccines"
    meta["collection_country"] = "United Kingdom"
    base = "no timepoint, infection status or symptom outcome in the deposited files"
    meta["notes"] = [base if not n else f"{base}; {n}" for n in meta["notes"]]
    first = [c for c in ["subject_id", "lc_status", "case_definition", "sample_matrix", "platform", "control_type", "sample_phase", "timepoint",
                         "time_anchor", "acute_severity", "hospitalised", "vaccinated_before_infection", "variant_era", "collection_country", "notes"] if c in meta.columns]
    meta = meta[first + [c for c in meta.columns if c not in first]]

    extra = {
        "n_peptide_rows": int(P.shape[0]),
        "n_dropped_unmatched_wells": int(len(unmatched)),
        "n_participant_id_conflicts": int((w[1] != w[2]).sum()),
        "n_subjects_multi_sample": int((nsub > 1).sum()),
        "n_samples_of_multi_sample_subjects": int(nsub[nsub > 1].sum()),
        "adapter_notes": adapted.notes,
    }
    return [write_cohort(STUDY, adapted, meta, extra)]


if __name__ == "__main__":
    for s in build():
        print(s)
