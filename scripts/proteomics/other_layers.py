"""Other omics layers measured on the same people, per proteomics cohort, from the papers and data
records (see docs/proteomics_cohorts.md). Used by build_all.py for the cohort table."""

OTHER_LAYERS = {
    "Buggert2025": "no (flow cytometry, serology on request only)",
    "Cervia2024_SomaScan": "yes: MS proteomics (Cervia2024_MS, same sera at 6 months); monocyte scRNA-seq, PBMC flow, VirScan (Mendeley 10.17632/dvf6yvrg4x.1, subsets); Olink Inflammation (Ceglarek 2025)",
    "Cervia2024_MS": "yes: SomaScan (Cervia2024_SomaScan, same sera)",
    "Su2022": "yes: plasma metabolomics (Table S2.2), PBMC scRNA-seq/CITE-seq (ArrayExpress E-MTAB-10129); TCR-seq without accession",
    "Kovarik2023_plasma": "yes, per sample in the paper's supplementary tables: eicosanoids, Biocrates metabolomics, triglycerides; no repository accession",
    "Captur2022": "none public",
    "Chowdhury2024": "no (glutathione, 8-OHdG, ELISA, flow cytometry; PCC participants from the BQC19 biobank, no omics link given)",
    "Zoodsma2022_Hannover": "no omics (flow counts, serology in the annotation; Bonifacius 2021)",
    "Zoodsma2022_Breda": "none public",
    "Zoodsma2022_Nijmegen": "no omics (flow cytometry; Janssen 2021)",
}
