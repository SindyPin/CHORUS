# Data inventory: transcriptomic and proteomic cohorts

Queried on 2026-11-09 from NCBI E-utilities and GEO SOFT (series, samples, platforms, BioProject, SRA), the ENA
portal `filereport` API (read runs, FASTQ availability), the EGA metadata API and PRIDE `ws/archive/v2`
(with ProteomeXchange as fallback). Only metadata were retrieved; no sequencing data were downloaded. The full
record for each accession is in `docs/data_inventory.csv`. The `provenance` column says how the matrix will be
built: `requantify` (public reads, nf-core/rnaseq), `author_matrix` (reads not public, author counts used) or
`author_table` (affinity proteomics, author NPX table).

## Transcriptomics (`configs/lc_transcriptomics_blood.yaml` plus candidates)

| study_id | matrix | n (GEO) | LC / control labels in GEO | raw reads | author matrix | provenance | status |
|---|---|---|---|---|---|---|---|
| GSE215865 | whole blood | 1392 libraries | none (acute cohort; LC labels in Synapse syn35874390) | yes | raw counts | requantify | usable once Synapse labels obtained |
| GSE275334 | PBMC | 47 | 15 LC / 18 healthy / 14 ME/CFS | **no** (NanoString) | RCC + normalised xlsx | author_matrix | **does not match config; drop from bulk_rna** |
| GSE267625 | whole blood | 111 | none (3 mo 79, 12 mo 32) | yes | csv (type not verified) | requantify | needs labels from P4O2 / paper |
| GSE169687 | whole blood | 152 | 14 healthy, convalescents 12–24 wpi | yes | RAW.tar | requantify | sensitivity contrast only |
| GSE270045 | whole blood | 36 | 19 LC / 17 **healthy** controls | yes | counts tsv | requantify | check controls were infected |
| GSE226260 | PBMC + whole blood | 331 | WB: 52 PASC / 51 recovered; PBMC: 72 PASC / 14 NOPASC; 142 acute | yes | raw counts | requantify | include; WB subset fits primary |
| GSE251849 | PBMC | 23 | 11 LC / 5 convalescent / 7 healthy | yes | counts txt | requantify | include (PBMC set) |
| PRJNA1184005 | PBMC | 14 runs | not in GEO (BioSample) | yes | none | requantify | include (PBMC set) |
| GSE224615 | PBMC | 36 | 23 LC / 13 non-LC | yes | **none** (DEG table only) | requantify | include; must requantify |
| GSE222253 (new) | whole blood | 35 | endotype only; symptom labels at BQC19 | yes | RAW.tar | requantify | conditional candidate |
| GSE286325 (new) | PBMC scRNA | 8 runs | 44 post-COVID / 44 convalescent (multiplexed) | **no** (0-read runs) | RAW.tar | author_matrix | pseudobulk only |

Candidates checked and not added: GSE152418 (acute COVID-19 vs healthy PBMC, no Long COVID), GSE265753 (scRNA of
long COVID patients in a herbal medicine trial, two healthy controls, no recovered controls), GSE235938 (scRNA,
two PASC and one control). GSE226260 and GSE224615 were already in the config and meet the criterion.

## Proteomics (`configs/lc_proteomics_plasma.yaml`)

| study_id | accession in config | finding |
|---|---|---|
| Zoodsma2022 | GitHub | repository name only, nothing to query |
| Buggert2025 | Zenodo | repository name only, nothing to query |
| Su2022 | Mendeley | repository name only, nothing to query |
| EatonFitch2024 | GSE275334 | series holds NanoString PBMC RNA, **no Olink data** |
| Vijayakumar2022 | Mendeley | repository name only, nothing to query |
| Liew2024 | EGAD00001009745 | **wrong accession**: EGA returns a hepatocellular carcinoma dataset |

The config has no PRIDE accessions. As a side check, the MS cohorts named in the project brief (T2) resolve as
follows: PXD036590 (Captur, SRM/MRM plasma, 43 files) and PXD036597 (Pushalkar, SWATH serum, 75 files) and
PXD047623 (Chowdhury, DIA plasma, 55 files) are in PRIDE; PXD026442 (Chen) and PXD044441 (Li 2024) exist in
ProteomeXchange but not PRIDE (another member repository, likely iProX); PXD052898 (Song) was not found in either.

## Blockers

1. **GSE275334 is not the Eaton-Fitch cohort described in the configs.** It is a 47-sample NanoString panel with
   no RNA-seq reads, no recovered controls and no Olink. Both configs need the correct accession (or the cohort removed),
   and the RNA-Olink bridge cohort for mosaic integration currently has no public source.
2. **Liew2024 EGAD00001009745 is an unrelated EGA dataset.** The PHOSP-COVID proteomics accession and access route
   must be identified before applying for access.
3. **Labels are missing from GEO for several core cohorts:** GSE215865 (Synapse syn35874390), GSE267625 (P4O2 data
   committee), GSE169687 (LC status of convalescents), GSE222253 (BQC19), PRJNA1184005 (BioSample attributes, not yet
   checked). Metadata curation (CC3) depends on these.
4. **Control definitions disagree with the config:** GSE270045 controls are labelled healthy, not COVID recovered,
   and GSE169687 has only never-infected controls. Several config sample counts (GSE226260, GSE224615, GSE215865)
   do not match GEO.
5. **Raw reads not public:** GSE275334 (NanoString, no reads) and GSE286325 (runs registered with zero reads) can
   only use author matrices (`author_matrix`). GSE224615 has public reads but no author matrix, so it has to be re-quantified.
6. **Proteomics cohorts on GitHub, Zenodo and Mendeley** have no record identifiers in the config, so their
   availability and sample counts are unverified.
7. **Matrix heterogeneity:** with GSE275334 removed, only GSE215865 (acute), GSE267625, GSE169687, GSE270045, the
   GSE226260 whole blood subset and GSE222253 are whole blood; the rest are PBMC.
