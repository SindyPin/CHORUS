# Literature stream 6: clinical metadata and phenotype harmonisation for Long COVID cohorts (searched 4 Nov 2026)

## Papers, standards and tools
Case definitions
1. WHO 2021 clinical case definition (Soriano JB et al., Lancet Infect Dis 2022;22:e102, doi 10.1016/S1473-3099(21)00703-9; unverified this session). Symptoms at 3 months lasting at least 2 months.
2. NICE NG188 (ongoing symptomatic 4 to 12 weeks; post COVID syndrome over 12 weeks); CDC 4 weeks (unverified this session).
3. NASEM 2024 definition; Chu L et al., J Gen Intern Med 2025;40:2024, doi 10.1007/s11606-025-09415-8. Infection associated chronic condition present at least 3 months; no positive test or specific symptom required.
4. Thaweethai T et al., JAMA 2023 (PASC research index) and Geng LN et al. 2024, JAMA, doi 10.1001/jama.2024.24184 (11 symptoms, threshold 11, likely and possible LC tiers, five clusters).
5. Lovaglio PG et al., Front Epidemiol 2025;5:1597799, doi 10.3389/fepid.2025.1597799. Prevalence 67 to 80 percent across definitions in one centre; risk factors change with definition.
6. Wisk LE et al., JAMA Netw Open 2025, doi 10.1001/jamanetworkopen.2025.26506. Five definitions on INSPIRE (n about 4,700): prevalence 15 to 42 percent; up to a third of between study variance from definition.
7. Foulkes AS, Thaweethai T, Reeder HT, arXiv 2608.04919 (2026). Statistical considerations: infection as pseudo label, attenuation from low specificity definitions, left censoring of onset, confounding of time since infection with variant era; three level LCRI categorisation.

Common data elements and cross cohort harmonisation
8. Pike Welch HL et al., JAMIA Open 2025;8(3):ooaf046, doi 10.1093/jamiaopen/ooaf046. RADx-UP Long COVID CDE task force; RECOVER CDEs reusable.
9. Sigfrid L et al., BMJ Open 2021;11:e043887. ISARIC harmonised follow up protocol (discharge anchored).
10. Pfaff ER et al., Lancet Digit Health 2022, doi 10.1016/S2589-7500(22)00048-6 (N3C ML phenotype); Reese JT et al., eBioMedicine 2023, doi 10.1016/j.ebiom.2022.104413 (subtypes across N3C and RECOVER).
11. Rinaldi E et al., npj Digit Med 2022;5:75, doi 10.1038/s41746-022-00620-x (ORCHESTRA, 743 core elements, SNOMED CT, LOINC, ATC, ICD-10); Rinaldi E et al. 2024, Applied Medical Informatics 46(Suppl 2) (long COVID variables mapped to FHIR and OMOP).
12. Munblit D et al., Lancet Respir Med 2022;10:715, doi 10.1016/S2213-2600(22)00169-2. PC-COS core outcomes.
13. FDA COVID-19 RWD data elements harmonisation (no PASC elements); CDISC COVID-19 TAUG (content not verified).
14. PHOSP-COVID cohort profile, medRxiv 2023.05.08.23289442 (n 7,935; discharge anchored visits).
15. Sun J et al. (RECOVER), Front Syst Biol 2025;4:1422384, doi 10.3389/fsysb.2024.1422384. Multi omics strategy: matrices and timepoints (enrolment, 90, 180 days from index).

Metadata curation tools
16. Bernstein MN et al., MetaSRA, Bioinformatics 2017 (PMID 28535296).
17. Sundaram SS, Musen MA et al., arXiv 2404.05893 (BioSample GPT-4 plus CEDAR); Sundaram, Goncalves, Musen, arXiv 2504.05307 (2025); Hardi J et al., arXiv 2604.08552 (2026, ARMS).
18. Multi agent AI system for metadata curation, bioRxiv 10.1101/2025.06.10.658658 (GEO scRNA-seq, precision above 0.98 on most fields).
19. Pavlidis P et al., bioRxiv 10.64898/2026.07.30.741874. Gemma curation agent; factor recall 0.98; about USD 0.5 per dataset.
20. Verbitsky A, Boutet P, Eslami M, Bioinform Adv 2025;5:vbaf241, doi 10.1093/bioadv/vbaf241.
21. MetaHQ (Krishnan lab), arXiv 2602.07805, github krishnanlab/meta-hq.
22. LLM-GEO-Label-Extractor (github SciSpectator).
23. EMBERS, bioRxiv 10.1101/2024.10.26.620145.
24. doppelgangR (Waldron, Bioconductor): duplicate sample detection across expression datasets.

Variant era and vaccination
25. Carazo S et al., Clin Infect Dis 2026;82:e436, doi 10.1093/cid/ciaf549; Babicki M et al., Pathogens 2024;13:267, doi 10.3390/pathogens13030267.

## Proposed minimal harmonised schema (implemented in chorus.schema)
case_definition_used {WHO2021, NICE_NG188, CDC_4wk, NASEM2024, RECOVER_LCRI2023, RECOVER_LCRI2024, ICD10_U09.9, self_report, study_specific}; lc_status {LC_likely, LC_possible, LC_negative, unknown}; lcri_score 0 to 50; control_type {infected_recovered, never_infected_healthy, never_infected_disease_comparator, acute_COVID}; infection_confirmation {PCR, antigen, serology, clinical_only, none}; days_since_infection; time_anchor {symptom_onset, first_positive_test, hospital_admission, hospital_discharge, enrolment}; acute_severity_who 0 to 10; hospitalised {yes, no, ICU}; age_years; sex; bmi; vaccination_doses_before_infection and last_dose_days_before_infection; variant_era {pre_Alpha, Alpha, Delta, Omicron_BA1_2, Omicron_BA4_5, Omicron_XBB_later, unknown} plus infection_calendar_month; reinfection; comorbidities (SNOMED or ICD-10 codes plus flat flags); symptom_domains (PC-COS); sample_matrix (RECOVER list plus UBERON); pre_analytics; source_accession and participant_id_hash; curation_provenance (manual or LLM_assisted model and date, supporting quote).

## Gaps
No published reusable harmonised metadata schema for Long COVID omics reuse; RECOVER CDEs are for prospective collection, ORCHESTRA mapping has no public value sets, PC-COS covers outcomes. No LLM curation effort models disease specific fields (case definition, time since infection, comparator type) or reports a Long COVID use case. No weak supervision or latent class approach applied to pooled Long COVID classification. No participant level de-duplication method for overlapping published cohorts.
