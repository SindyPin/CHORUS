# CHORUS project brief: harmonising Long COVID omics cohorts for externally validated machine learning

Prepared 4 November 2026 for Sindy Pinero, updated 5 November 2026. Working title and package name CHORUS (Cross cohort Harmonisation of Omics for Robust and Unified early detection Signatures). PyPI name `choruslc` is free; `cohere` and `chorus` are taken. Rename freely before the first public release.

Contents

1. Verdict in one page
2. Conference fit (ISMB/ECCB 2027)
3. Literature review synthesis and novelty check
4. The computational framework
5. Feasibility, data inventory and compute plan
6. Ten week plan to the 14 January 2027 deadline
7. Prompts for Claude Code and Claude Science
8. Repository, package and GitHub
9. Paper plan and LaTeX template
10. Reviewer audit checklist (run before submission)
11. Reusable master prompt for future papers

---

## 1. Verdict in one page

The idea is novel and feasible for the ISMB/ECCB 2027 proceedings if it is framed as a validated framework plus benchmark rather than as a new batch correction method.

What exists. Pooled re-analyses of Long COVID transcriptomes exist (Maison 2025, 21 series, all 23 LC samples from one study; Avdiu 2025 thesis, 10 series, same confounding) but cannot measure cross study generalisation. One lineage does multi dataset Long COVID transcriptomic ML with cross dataset testing (Lopez-Rincon 2025, Varga 2026) but applies no batch harmonisation and rests on a 33 sample discovery set. One cross cohort proteomic signature exists (Maestri and Kwon 2025, nine public Olink cohorts, summary AUROC 0.78 discovery, 0.76 validation) but produces a score, not a harmonised resource or a locked classifier. Two independent replication attempts failed (Farztdinov 2024 complement; Gao 2025 two country Olink). The largest single consortium resources (PHOSP-COVID n 657, Gutenberg n 495, IMPACC n 513) all stop at internal cross validation.

What does not exist. A harmonised, sample level, multi cohort Long COVID omics resource with a leave one study out benchmark; any Long COVID paper that fits the harmoniser inside the training fold; any omics paper that treats observational cohorts as environments for invariant feature selection and tests the transfer benefit; any framework coupling per layer harmonisation, mosaic multi omics integration and cross study predictive validation with harmonised clinical covariates.

The defensible claim. "First framework to treat Long COVID cohorts as environments, harmonise them without leakage, and show under strict leave one study out validation whether invariant and causally supported features transfer better than variance selected features." Claim the combination and the evaluation, never the components (ComBat, ICP, TabPFN and causal batch correction all exist).

The key risk. Nastl and Hardt (NeurIPS 2024) found that causal feature subsets did not generalise better on 16 tabular benchmarks. The paper must be framed as a test, must include all features and differential expression baselines under identical harmonisation, and must report the within study versus cross study gap. A negative or mixed result on the invariance question is still publishable because the harmonisation benchmark and resource stand on their own.

Where your previous models sit. MRCONTROL, TACO and SPLIT all use one cohort (GSE215865, Mount Sinai). TACO reached ROC AUC 0.675 within cohort. CHORUS is the external validation step all three papers list as future work, and the harmonised resource feeds SPLIT and MRCONTROL afterwards.

---

## 2. Conference fit (ISMB/ECCB 2027)

Source: iscb.org/ismbeccb2027 pages fetched 4 November 2026.

- Dates: 18 to 22 July 2027, Bella Center Copenhagen, hybrid.
- Proceedings deadline: 14 January 2027, 23:59 AoE, no extensions. Submission through EasyChair (ismbeccb2027). Typical cycle from earlier years: conditional acceptance early March, revised papers mid March, final notification early April. Payment system opens 5 April 2027; fee USD 1052 due 6 May 2027.
- Format: OUP Bioinformatics template only, nine pages maximum including references, structured abstract of at most 250 words with headings Motivation, Results, Availability and Implementation, Contact, Supplementary Information. Alt text on figures. References in Oxford SCIMED style as a .bib. Font or spacing changes are grounds for rejection. Supplementary material discouraged.
- Not double blind: author information required. Your GitHub does not need anonymising.
- Code: a public, documented repository must be linked in the abstract at submission; missing code or reproducibility is grounds for rejection.
- Review: three reviews, two rounds (accept or minor; major changes with response letter; reject). Highly ranked non accepted papers may be offered transfer to Bioinformatics Advances.
- Published in an open access section of Bioinformatics; 20 minute talk; in person strongly preferred.

Area choice. First choice Biomedical Informatics; second choice General Computational Biology (Systems Biology and Networks also fits the mosaic and causal parts). COSIs to tick (up to three): MLCSB (Machine Learning in Computational and Systems Biology), TransMed (Translational Medicine Informatics), and SysMod or BioInfo-Core. Confirm the 2027 COSI list on the submission form.

Fallbacks. Abstract track for talks and posters opens about late January with a mid April deadline; the LC harmonised resource would also fit the Tech Track (software) in May.

---

## 3. Literature review synthesis and novelty check

Six parallel reviews were run on 4 November 2026 (Long COVID multi cohort omics; transcriptomic cross study harmonisation; proteomics and metabolomics harmonisation; mosaic multi omics integration; causal and invariant learning across cohorts; clinical metadata harmonisation). Full outputs are in `docs/literature/`. Headline findings and the references to cite:

Long COVID multi cohort work. Maison 2025 (Mitochondrion 85:102072); Varga 2026 (Front Syst Biol 6:1873899) and Lopez-Rincon 2025 (medRxiv 2025.02.07.25321856); Eroglu and Firat 2026 (Zenodo 21398550; ComBat on all six GEO cohorts then external AUC 0.984, a leakage pattern, cite with caution); Maestri and Kwon 2025 STOP-PASC (medRxiv 2025.12.04.25341650, nine cohort Olink meta-analysis); Farztdinov 2024 (medRxiv 2024.03.14.24304224); Gao 2025 (Nat Immunol 26:692); Aid 2025 (Nat Immunol 27:61, GSE226260, independent validation cohort); Gabernet 2025 (JCI 193698, IMPACC SPEAR versus MOFA); Liew 2024 (Nat Immunol, PHOSP-COVID); Ten Cate 2025 (Research Square, Gutenberg proteomic atlas, LOOCV only); Klein 2023 (Nature, MY-LC); Lammi 2025 (Nat Genet, 24 cohorts). No Long COVID paper performs leave one study out validation.

Harmonisation methods for prediction. Addon adjustment (Hornung 2017 Bioinformatics, bapred; Parker 2014 PeerJ frozen SVA); dangers of ComBat with group covariate under imbalance (Nygaard 2016; Zindler 2020); batch confounded with class biases CV (Soneson 2014); no pipeline generalised to all external sets (Van 2024 BMC Bioinformatics 25:181); supervised correction leaks labels (Sun 2024 bioRxiv; Hui 2024 Brief Bioinform); FSQN (Franks 2018); TDM (Thompson 2016); gene pair methods (Geman 2004; Wu 2025 review); MINT (Rohart 2017); cross study learning and merging versus ensembling (Patil and Parmigiani 2018 PNAS; Guan 2019); cross study validation (Bernau 2014); bioLeak (Korkmaz 2026 arXiv 2604.10965, leakage audit with 0.2 AUC inflation); uniform reprocessing necessary but not sufficient (Wang 2018 Sci Data; Arora 2020); whole blood composition (O'Connell 2023; Gautam 2019); batch effects on multi study classifiers (Raab, Johnson, Piccolo 2026 bioRxiv, supervised ComBat can invert class signal).

Proteomics and metabolomics. OlinkAnalyze bridging vignettes (no bridging without shared samples); Explore HT validation (median r 0.88 to 3072); Rooney 2025 Clin Chem (HT precision); Heylen 2024 Brief Bioinform (SPOC; rank INT unreliable across cohorts); SCALLOP analysis plan (effect size pooling); Pietzner 2021, Katz 2022, Eldjarn 2023, Wang 2025, Barente 2026 (Olink versus SomaScan; model transfer loses performance; correlation cutoff 0.3 to 0.4); cpiVAE 2025; quantms (Dai 2024 Nat Methods); HarmonizR (Voss 2022); proBatch (Cuklina 2021); protein level correction and ratio to reference (Chen 2025 Nat Commun; Yu 2023 Genome Biol Quartet); round robin 71 reproducible plasma proteins (Kardell 2025); Sissala 2025 Olink versus MS median rho 0.59; RefMet, metLinkR (Patt 2025), MMS (Villalba 2023).

Mosaic integration. MOFA+ groups (Argelaguet 2020); MOTL transfer (Hirst 2025 Genome Biol); StabMap (Ghazanfar 2024 Nat Biotechnol); MIDAS (He 2024); Multigrate; SC5 (Arriola 2024, cohort as domain); SCMMIB benchmark 2025; MSB late fusion stacking (Boussena 2026); MAGNET, MOIRA, MIMIR, Flex MoE for missing modalities; fedRBE and Fed ComBat for federated settings; SurvBoard (Wissel 2025); Scicluna 2025 Nat Med (consensus sepsis subtypes validated in external cohorts). No benchmark compares mosaic methods on cross cohort generalisation of disease classifiers.

Causal and invariant learning. Peters 2016 ICP; Rojas-Carulla 2018; Arjovsky 2019 IRM; Rothenhausler 2021 anchor regression; Pfister 2021 stabilised regression (proteomics application); Wu 2025 Bayesian invariance (6169 features); Londschien 2026 (anchor regression with ICU database as environment, 1 to 3 percent OOD gains); Nastl and Hardt 2024 (negative result); Perez 2025 (domain invariant single cell patient prediction); Ng 2025 InvGraphCox (cohorts as environments, cancer survival); Yu, Liu, Li, Ding, Le 2019 TPAMI MCFS (multi source causal feature selection, Thuc is a coauthor, natural to cite and possibly reuse); Bridgeford 2025 (causal perspective on batch effects, causalBatch); Wang 2021 flow based causal harmonisation; Mooij 2020 JCI; Huang 2020 CD-NOD; Jablonski 2022 DCE; TabPFN under shift (Cheng 2025; Loza 2026).

Clinical metadata. WHO 2021, NICE NG188, CDC, NASEM 2024 (Chu 2025 JGIM), RECOVER index 2023 and 2024 update (Geng 2024 JAMA); definition changes prevalence and risk factors (Lovaglio 2025; Wisk 2025 JAMA Netw Open); Foulkes 2026 arXiv (statistical considerations, pseudo labels, left censoring, variant era confounding); RECOVER CDEs (Pike Welch 2025); ISARIC follow up protocol (Sigfrid 2021); ORCHESTRA harmonisation and FHIR/OMOP mapping (Rinaldi 2022, 2024); PC-COS core outcomes (Munblit 2022); N3C phenotype (Pfaff 2022); LLM assisted metadata curation (Sundaram and Musen 2024, 2025; multi agent GEO curator 2025; Pavlidis 2026; Verbitsky 2025; MetaHQ 2026); doppelgangR for duplicates. No published harmonised metadata schema exists for Long COVID omics reuse and no LLM curation effort covers disease specific fields such as case definition or time since infection.

Novelty verdict. Novel as a combination and as a validated resource. Not novel as a batch correction method, as causal feature selection, or as TabPFN under shift. The Introduction should say exactly this.

---

## 4. The computational framework

Five stages, all implemented in the package (`src/chorus`).

Stage 1, clinical harmonisation (`chorus.schema`). Minimal schema with controlled vocabularies: sample and subject identifiers, study, LC status (LC, COVID recovered, never infected, disease comparator), case definition used by the study (WHO 2021, NICE NG188, CDC 4 week, NASEM 2024, RECOVER index, self report, study specific), sample matrix, omics layer; optional days since infection with time anchor, acute severity in four WHO bands, hospitalisation, age, sex, BMI, vaccination before infection, variant era, comorbidities, symptom domains (PC-COS), platform, overlapping study. LLM drafted, human verified curation with provenance per field (`scripts/hpc/00_curate_geo_metadata.py`). Primary contrast LC versus COVID recovered. Sensitivity analysis on WHO compatible cases. Subject de-duplication across series.

Stage 2, per layer technical harmonisation (`chorus.harmonise`). Uniform preprocessing per layer: bulk RNA re-quantified from FASTQ with nf-core/rnaseq and Salmon against one GENCODE release, log CPM, gene filter applied per study and intersected, analysis within sample matrix; Olink restricted to UniProt assays shared across panel generations, within study standardisation or ranks only; MS re-searched with DIA-NN or FragPipe against one FASTA, protein level median centring, missingness indicators kept; metabolomics mapped to RefMet with MSI level recorded. Harmonisation operators with fit on training studies and label free addon transform on the held out study: identity, per study centring and scaling, quantile to training reference, feature specific quantile normalisation, within sample rank, top scoring pairs, addon ComBat with protected covariates. Transductive ComBat (inmoose) included only as a flagged comparator.

Stage 3, cohort as environment feature selection (`chorus.features`). Per study Hedges g, random effects meta-analysis (DerSimonian and Laird), invariance screen on pooled p, I squared, Cochran Q and sign consistency, ranking by pooled z. Baselines: most variable genes after per study centring (the TACO baseline), meta-analytic differential expression, all features. Causal priors from MRCONTROL (32 genes), TACO (411 DCE genes, 6 gene signature), SPLIT consensus drivers, alone or intersected or unioned with the invariant screen. On simulated data the invariant screen selected 15 of 20 true causal genes and zero of 40 study specific shortcut genes, while differential expression selected 12 causal and 8 shortcut genes.

Stage 4, mosaic multi omics integration (`chorus.integrate`). Primary route: late fusion of per layer leave one study out probabilities with layer presence indicators (missingness aware stacking). Comparators: diagonal integration of RNA and protein cohorts at gene product level through Ensembl to UniProt mapping with ranks; MOFA+ with cohorts as groups and layers as views, held out samples projected onto training loadings. Bilayer bridge cohorts: Eaton-Fitch 2024 (RNA plus Olink), Su 2022 (scRNA plus Olink), Li 2024 (MS plus scRNA), Cervia 2024 (affinity plus MS plus scRNA), Ai 2024 (controlled).

Stage 5, evaluation (`chorus.eval`). Leave one study out with harmoniser, selector and classifier fitted inside the training fold; nested inner leave one study out to choose the harmoniser; classifiers elastic net, random forest, TACO ensemble (TabPFN plus four), cross study ensemble (Patil and Parmigiani); metrics AUROC with stratified bootstrap CI, AUPRC, balanced accuracy, Brier, expected calibration error; within study five fold CV reported alongside as the optimism reference. Note that AUROC within a held out study is invariant to an additive shift of the score, so location only study effects hurt calibration and thresholded metrics more than AUROC; report Brier and balanced accuracy, not AUROC alone.

Primary experiment (T1). LC versus COVID recovered detection from blood bulk RNA-seq, whole blood primary, PBMC secondary, pooled sensitivity. Candidate cohorts (from your catalogue plus the review): GSE275334 Eaton-Fitch (165/495/72), GSE267625 Blankestijn (77/18), GSE215865 Thompson (acute samples with later LC labels; keep as the TACO replication and the acute transfer set), GSE226260 Aid 2025 PBMC (124 PASC/63 recovered, new), GSE270045 Shahbaz (19/17), GSE224615 Yin bulk (27/16), GSE251849 Greene PBMC (11/5/7), PRJNA1184005 Missailidis PBMC (7/7), GSE169687 Ryan (69 LC, never infected controls only; sensitivity contrast). Expect six to eight usable cohorts, about 1,100 to 1,400 participants.

Secondary experiment (T2). Plasma proteomics: Olink family (Zoodsma GitHub 186/350/61; Buggert Zenodo 70/70; Su Mendeley 309/209/457; Vijayakumar Mendeley 41/38/29; Eaton-Fitch Olink via GSE275334; Liew EGA controlled 424/233) and MS family (Chen PXD026442 118/94/109; Pushalkar PXD036597 29/29; Chowdhury PXD047623 25/0/9; Song PXD052898 50/0/10; Captur PXD036590 targeted 54/0/102; Li 2024 PXD044441 60/20). Compare directly with Maestri 2025 (AUROC 0.76 in validation cohorts).

Tertiary (T3). Mosaic integration on bilayer cohorts; transfer of pooled post acute models to acute phase GSE215865 against TACO within cohort.

---

## 5. Feasibility, data inventory and compute plan

Data access. Of the 101 catalogued studies, 52 have public or controlled access data. The transcriptomic arm needs only GEO or SRA public data (verify per series whether FASTQ are on SRA or restricted to dbGaP; the Mount Sinai series may provide counts only, which is acceptable with provenance flagged). Controlled access sets (Liew EGA, Balnis dbGaP, Ai GSA-Human) take weeks to months: apply now, treat as a revision or follow up addition.

Compute (you have HPC). Bulk RNA re-quantification: about 1,000 to 1,400 whole blood or PBMC libraries at 20 to 60 million reads, Salmon pseudo alignment via nf-core/rnaseq, roughly 0.5 to 1.5 CPU hours per library, 1,000 to 2,000 CPU hours total, transient FASTQ storage 5 to 15 TB (delete after quantification, keep quant files). Download from ENA with aria2c is the bottleneck; expect two to five days of transfers. DIA-NN re-search of five to six PRIDE plasma datasets (roughly 300 to 600 runs): one to three days on a 32 core node, or GPU. Benchmark itself: seconds to minutes per configuration; the full grid with TabPFN runs in hours on one GPU node.

Risks and mitigations. (i) Few cohorts have recovered controls at the same time window: pool across windows with days since infection as a covariate and report window stratified results. (ii) Whole blood versus PBMC: analyse within matrix first. (iii) Label heterogeneity: WHO strict sensitivity analysis. (iv) Participant overlap: de-duplicate; note that INCOV and Mount Sinai feed several papers. (v) Small cohorts: minimum five cases and five controls to be a held out set; small sets still contribute to training. (vi) If invariance does not help, report it; the harmonisation benchmark and resource carry the paper.

---

## 6. Ten week plan to 14 January 2027

Week 1 (4 to 10 Nov): confirm cohort list and accessions on the HPC (prompt CC1); submit controlled access applications; create the GitHub repository and push the package; agree the preregistered primary comparison with Thuc.
Weeks 2 and 3: download and re-quantify transcriptomic cohorts (prompt CC2); curate metadata with LLM drafts and verify every field (prompt CC3); export TACO, MRCONTROL and SPLIT gene lists into `priors/`.
Week 4: run T1 benchmark on the HPC (prompt CC4); first look at results in Claude Science (prompt CS1); draw Figure 2 (simulation) and Figure 3 (T1).
Weeks 5 and 6: proteomics arm (prompt CC5), MS re-search in parallel; T2 benchmark; comparison with Maestri 2025.
Week 7: mosaic integration and acute transfer (T3); stability and enrichment analyses of the invariant signature (prompt CS2).
Week 8: write Results and Discussion in the MRCONTROL register; freeze the repository with a tagged release and Zenodo DOI; make the BioRender overview figure.
Week 9: internal review by Thuc; run the reviewer audit checklist (Section 10); fix.
Week 10 (to 14 Jan): final compile in the OUP template, nine page check, alt text, EasyChair submission. Christmas falls in week 8; plan the writing before it.

---

## 7. Prompts for Claude Code and Claude Science

Use Claude Code (Fable, high effort for the pipelines; Opus, medium for small edits) inside the CHORUS repository clone on the HPC login node or a machine that can submit Slurm jobs. Each prompt is self contained.

### CC1. Verify accessions and raw data availability

```
You are working in the CHORUS repository (github.com/SindyPin/CHORUS). Read configs/lc_transcriptomics_blood.yaml and configs/lc_proteomics_plasma.yaml. For every accession listed, query NCBI E-utilities (GEO, SRA, BioProject), ENA (filereport API) and PRIDE (ws/archive/v2) and produce docs/data_inventory.csv with columns: study_id, accession, repository, n_samples, library_strategy, platform, raw_reads_public (yes, no, dbGaP or controlled), supplementary_count_matrix (file name or none), sample_characteristics_fields (comma separated), notes. Also check these additional candidate series and add them if they contain Long COVID blood transcriptomes with COVID recovered controls: GSE226260, GSE224615, GSE152418, GSE222253, GSE265753, GSE235938, GSE286325. Do not download sequencing data in this step. Flag any series whose raw reads are not public so that the author count matrix is used instead with provenance author_matrix. Summarise the inventory in docs/data_inventory.md with a table and a short list of blockers.
```

### CC2. Uniform re-quantification on the HPC

```
In the CHORUS repository, adapt scripts/hpc/01_fetch_sra.sh and scripts/hpc/02_nfcore_rnaseq.slurm to this cluster: load the right modules (Nextflow, Singularity or Apptainer), set the scratch root, and set the Slurm partition and account. Then, for every transcriptomic cohort in docs/data_inventory.csv with raw_reads_public = yes, fetch FASTQ from ENA with aria2c, verify md5, build the nf-core/rnaseq samplesheet, and submit nf-core/rnaseq 3.19.0 with Salmon pseudo alignment against GRCh38 and GENCODE v47 (one shared reference built once with --save_reference). Pin the pipeline version and reference in docs/provenance.md. After each run, copy salmon.merged.gene_counts.tsv and the MultiQC report to data/quant/<GSE>/ and write a QC table (mapping rate, rRNA fraction, globin fraction, read count) to data/quant/<GSE>/qc.csv. Fail loudly if a cohort has fewer than 80 percent of samples passing a 60 percent pseudo alignment rate. Do not modify the package code in this task.
```

### CC3. Metadata curation into the CHORUS schema

```
In the CHORUS repository, for every transcriptomic cohort in docs/data_inventory.csv, download the GEO family SOFT file and run scripts/hpc/00_curate_geo_metadata.py to produce data/meta/<GSE>_curation.csv. Then run it again with --draft using the Anthropic API (model claude-opus-5-5) to fill draft_* columns with supporting quotes. Produce for me a verification worksheet data/meta/<GSE>_verify.csv showing raw characteristics, the draft value and the quote side by side for every field, sorted so that disagreements between the draft and a rule based parse of the characteristics come first. I will verify by hand. After I return the verified file, write data/meta/<GSE>_meta.csv with exactly the chorus.schema.ClinicalSchema columns, validate it with chorus.schema.validate_metadata(strict=True), map SRA run accessions to GSM sample ids so that sample_id matches the count matrix columns, and record curation_provenance. Also build data/meta/subject_overlap.csv listing subjects that appear in more than one series (shared GSM titles, consortium names, identical identifiers) for de-duplication.
```

### CC4. Run the primary benchmark

```
In the CHORUS repository on the HPC, install the package with pip install -e ".[all]" in a fresh conda environment (python 3.11, torch with CUDA for TabPFN), run pytest, then run chorus run configs/lc_transcriptomics_blood.yaml --within as a Slurm job on a GPU node with 16 CPUs and 64 GB. Before running, export the TACO 411 DCE genes (github.com/SindyPin/TACO, DCE/DCE_COVID_LC_Results_annotated.csv) and the MRCONTROL 32 genes (github.com/SindyPin/Causal-Multiomics-Method) into priors/taco_411_dce_genes.txt and priors/mrcontrol_32_genes.txt as Ensembl gene ids without version, using the GENCODE v47 GTF for symbol to Ensembl mapping where needed. After the run, produce results/lc_transcriptomics_blood/report.md with: the cohort summary table; the LOSO summary (mean AUROC per harmoniser x selector x classifier); the per held out cohort AUROC with CIs for the preregistered comparison (nested harmoniser with invariant selection versus identity with variance selection, elastic net); the within study CV optimism gap; the harmoniser chosen by the nested loop per held out cohort; the number of strictly invariant genes per fold and their Jaccard overlap across folds; the overlap of invariant genes with the TACO, MRCONTROL and SPLIT lists; and Brier and balanced accuracy next to AUROC. Then rerun with preprocessing.matrices [whole_blood, PBMC] and covariates [age, sex, sample_matrix] as the pooled sensitivity analysis, and with the strict WHO mask (chorus.schema.strict_who_subset) as a second sensitivity analysis. Save all configs used under results/ for provenance.
```

### CC5. Proteomics arm

```
In the CHORUS repository, assemble the plasma proteomics cohorts listed in configs/lc_proteomics_plasma.yaml. For Olink cohorts, download the processed NPX tables from the listed repositories (Zoodsma GitHub, Buggert Zenodo, Su and Vijayakumar Mendeley, Eaton-Fitch GSE275334 supplementary), map OlinkID to UniProt with the Olink assay lists for Target 96, Explore 384, 1536 and 3072, write data/prot/<study>_npx_uniprot.tsv (samples x UniProt) and the matching schema metadata, and compute the set of assays shared across all Olink cohorts with the fraction below LOD per cohort. For MS cohorts, adapt scripts/hpc/03_diann_plasma.slurm to this cluster and re-search PXD026442, PXD036597, PXD047623, PXD052898 and PXD044441 with DIA-NN 2.x (or FragPipe for DDA sets) against one UniProt human release; write report.pg_matrix.tsv per dataset and map runs to samples from the PRIDE SDRF files. Then run chorus run configs/lc_proteomics_plasma.yaml for the Olink family, a second config for the MS family, and a third config that unions both at UniProt level with rank based operators only. Produce results/lc_proteomics_plasma/report.md as in CC4, and add a comparison table against Maestri and Kwon 2025 (summary AUROC 0.78 discovery, 0.76 validation, 60 protein signature) including how many of their 60 proteins are in our invariant set.
```

### CC6. Package hardening before submission

```
In the CHORUS repository: add type hints and docstrings where missing, make the test suite cover every harmoniser, selector and classifier name in HARMONISERS, SELECTORS and CLASSIFIERS including optional dependencies skipped when absent, add a GitHub Actions workflow running pytest on Python 3.10 to 3.12, add a CONTRIBUTING.md, write docs/usage.md with the config schema, pin a requirements lock, tag v0.1.0, create a GitHub release and mint a Zenodo DOI, and build and upload the wheel to PyPI as choruslc. Do not change any default parameter used in the paper results; if a bug fix changes results, record it in CHANGELOG.md and rerun the affected config.
```

### CS1. Claude Science: first analysis of T1 results

```
Attach results/lc_transcriptomics_blood/loso_results.csv, loso_predictions.csv, within_study_cv.csv and cohorts.csv. (1) For the preregistered comparison (elastic net; nested harmoniser with invariant selection versus identity with variance selection), compute per held out cohort AUROC differences with bootstrap CIs and a paired Wilcoxon signed rank test across cohorts; report Brier and balanced accuracy differences too. (2) Build the harmoniser x selector heatmap of mean held out AUROC for each classifier, and a panel showing within study CV AUROC against held out AUROC per cohort (optimism gap). (3) Test whether held out AUROC under identity depends on the held out cohort's case to control ratio and sample size. (4) Summarise which harmoniser the nested loop chose per cohort. (5) Flag any configuration where the transductive combat_pooled comparator exceeds combat_addon by more than 0.03 AUROC as evidence of leakage inflation. Write a results narrative in the register of my PLoS Comput Biol paper (we sentences, connected sentences, hedged interpretation, no em dashes, no hyphens in compound modifiers) that I can paste into the Results section, with figure legends.
```

### CS2. Claude Science: the invariant signature

```
Attach results/lc_transcriptomics_blood/selected_features.json and the invariant tables per fold, plus priors/taco_411_dce_genes.txt, priors/mrcontrol_32_genes.txt and priors/split_consensus_drivers.txt. (1) Compute the consensus invariant gene set (selected in all held out folds) and the Jaccard overlap between folds. (2) Test enrichment of the consensus set for the TACO, MRCONTROL and SPLIT lists against the filtered gene universe (hypergeometric), and run GO BP, KEGG and Reactome over representation. (3) Map the consensus genes to the convergent biology of our 101 study review (complement and coagulation, interferon, neutrophil programmes, mitochondrial and energy metabolism) and report which pathways are supported, with hedged language. (4) For the proteomic arm, check whether the nine complement proteins (C4A, C4BPB, C5, CFHR1, CFHR5, MBL2, FCN2, MASP2, C1Q) pass the invariance screen. Produce a table and a short Results subsection in the MRCONTROL register.
```

---

## 8. Repository, package and GitHub

The package is written and tested (8 passing tests): `src/chorus` with modules schema, io, harmonise (operators, transcriptomics, proteomics), features (invariant, priors), integrate (mosaic), models (classifiers), eval (loso, metrics), simulate, cli; `scripts/hpc` with the SRA fetch, nf-core and DIA-NN Slurm templates and the metadata curation script; `configs` for T1 and T2; `priors`; `paper` with the LaTeX skeleton and references; `docs` with this brief and the literature outputs.

GitHub. This session cannot create a new repository on your account (the token is scoped to existing repositories). Create an empty public repository named CHORUS under SindyPin (no README), then tell me and I will push the package from here; or run in Claude Code: `cd CHORUS && git init && git add . && git commit -m "CHORUS v0.1.0 scaffold" && gh repo create SindyPin/CHORUS --public --source=. --push`. The conference is not double blind, so no anonymisation is needed; keep the repository public and documented at submission because reviewers inspect it.

---

## 9. Paper plan and LaTeX template

Template. OUP Bioinformatics via the official Overleaf "OUP General Template" (overleaf.com/latex/templates/oup-general-template/ybpypwncdxyb) with `\documentclass[unnumsec,webpdf,contemporary,large]{oup-authoring-template}` and `\journaltitle{Bioinformatics}`. `paper/main.tex` and `paper/references.bib` in the repository are written for that class: open the Overleaf template, replace the main file with `main.tex`, add `references.bib` and compile. I can push them into your Overleaf through Chrome when you have the browser connected; tell me the project name you want.

Structure (nine pages). Structured abstract; Introduction (four paragraphs, drafted); Methods (overview figure, cohorts and schema, preprocessing, harmonisers, feature selection, classifiers, evaluation, mosaic, acute transfer, implementation; drafted); Results (simulation; harmonised cohorts table; pooling without harmonisation; harmonisation inside the fold; invariant versus variance features; mosaic and acute transfer; placeholders); Discussion (four paragraphs outlined); Conclusion; declarations. Figures: Fig 1 overview (BioRender, MRCONTROL style), Fig 2 simulation (already produced by `scripts/simulation_benchmark.py`), Fig 3 T1 results (two panels), Fig 4 invariant signature and proteomic arm. Table 1 cohorts.

Register. Written in the MRCONTROL register per your preference: "We" sentences, connectives (Specifically, Moreover, In contrast), hedged interpretation, no rhetorical openers, no "not X but Y", formal vocabulary, no em dashes, no hyphens in compound modifiers in prose.

---

## 10. Reviewer audit checklist (run before submission)

Novelty and framing
- Does the Introduction state precisely what exists (Maison, Varga, Maestri, Lammi) and what does not, and claim only the combination plus the validated resource?
- Is Nastl and Hardt 2024 cited and answered with the all features and differential baselines under identical harmonisation?
- Are MRCONTROL, TACO and SPLIT cited as the single cohort predecessors and is CHORUS positioned as their external validation step?

Leakage and evaluation
- Is every harmoniser fitted on training cohorts only, with addon transform of the held out cohort documented in Methods and enforced in code (tests)?
- Is the transductive ComBat comparator clearly flagged and is its inflation relative to addon ComBat reported?
- Was the harmoniser chosen by nested inner leave one study out, never with the held out cohort?
- Are subjects de-duplicated across series and is the de-duplication count reported?
- Is the preregistered primary comparison stated before the results, with a paired test across held out cohorts?
- Are within study CV results shown next to held out results (optimism gap)?
- Are Brier, calibration and balanced accuracy reported, given that AUROC is invariant to score shifts within a held out cohort?
- Are confidence intervals on every held out AUROC and is the number of cases and controls in each held out cohort shown?

Data and harmonisation
- Table 1 lists accession, platform, matrix, counts per group, window, case definition, provenance (re-quantified or author matrix) for every cohort.
- Is the gene filter applied per study and intersected, and the number of genes reported?
- Are whole blood and PBMC analysed within matrix first, with the pooled analysis labelled as sensitivity?
- Is NPX never pooled on one scale and is the rank level interpretation stated for Olink?
- Is MS missingness handled without imputation before correction and are missingness indicators described?
- Are case definition heterogeneity and the WHO strict sensitivity analysis reported?
- Is time since infection modelled or stratified, with anchors recorded?

Causal and invariance claims
- Is the invariance screen described as a single feature stability screen, not full ICP, and is invariance described as necessary rather than sufficient for causality?
- Are the numbers of strictly invariant features per fold and their cross fold overlap reported?
- Are results reported for the cases where invariant selection did not help?

Reproducibility (ISCB requires public code at submission)
- Public repository with README, install, configs, provenance, tagged release and Zenodo DOI; tests pass in CI; versions of nf-core/rnaseq, GENCODE, DIA-NN, UniProt and TabPFN pinned.
- All curated metadata and all result tables deposited; controlled access data described with the access route.

Format
- Nine pages including references in the OUP template; structured abstract at most 250 words with the five headings; alt text on figures; .bib references; no font or spacing changes; corresponding author named; LLM use disclosed per ISCB policy (writing assistance only, as in TACO).

Style (your preferences)
- No em dashes; no hyphens in compound modifiers in prose; S, E and AE abbreviations where safety, efficacy and adverse events are discussed; MRCONTROL register throughout; no AI sounding vocabulary.

---

## 11. Reusable master prompt for future papers

Save this as your standing prompt for "paper from scratch" sessions and fill the bracketed fields.

```
You are my research co-author for a computational biology paper. Work in this chat as the coordinator, use Claude Science and Claude Code only through prompts you write for me to run, and tell me in one line at the start if a different model or effort level suits the task better.

Context about me: [name, role, group, prior methods and their repositories, writing register to follow, formatting preferences]. My shared materials are in [folder]. Read every file before proposing anything and summarise what each contributes.

Topic: [one paragraph: the scientific question, the data I have or can access, the models I ran before and their performance, the improvement I want].

Target venue: [conference or journal and URL]. Check the call: tracks or areas, COSIs, deadlines, page limit, template, anonymisation policy, code and data requirements, review process. Pick the best fitting track and say why; if nothing fits, propose two or three feasible alternatives.

Deliver, in this order, and stop to ask me only when a decision cannot be undone:
1. Literature review: run parallel searches on (a) the disease or domain specific prior art, (b) each methodological component, (c) adjacent fields that solved the same problem, (d) negative results that a reviewer would raise. Return citable references with DOIs, a precise statement of what exists and what does not, and a one paragraph defensible novelty claim that claims the combination and the evaluation rather than the components.
2. Framework: stages, inputs, operators, baselines, evaluation design (no leakage, held out units, nested selection, preregistered primary comparison, optimism reference), sensitivity analyses, and the biological interpretation plan.
3. Feasibility: data inventory with access tiers, compute estimate for my HPC, risks and mitigations, and a week by week plan to the deadline.
4. Package: scaffold an installable package in [language] with tests, configs, HPC scripts and a README; run the tests and a simulation benchmark here; name it after checking PyPI or CRAN availability.
5. Repository: create or prepare the GitHub repository, respecting the venue's anonymisation policy.
6. Prompts: self contained prompts for Claude Code (data retrieval, pipelines, benchmark runs, package hardening) and Claude Science (analysis of results, figures, results narrative in my register).
7. Paper: pick the venue's template, write main.tex and references.bib with the Introduction and Methods drafted in my register and Results placeholders tied to the planned figures and tables; prepare it for Overleaf.
8. Audit: a reviewer checklist covering novelty, leakage, statistics, data, reproducibility, format and my style rules; run it against the draft and list the fixes.

Rules: cite only what you verified; mark unverified DOIs; no em dashes; no hyphens in compound modifiers; formal academic register modelled on [my reference paper]; keep deliverables as files in [folder] and give me each file as it is ready rather than at the end.
```
