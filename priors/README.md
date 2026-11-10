# Causal priors

Place the published gene lists here (one identifier per line, `#` comments allowed):

- `taco_411_dce_genes.txt`: the 411 DCE genes from TACO (repository SindyPin/TACO, the 411 gene columns of
  Data/411_DCEG_Paths_average.csv), Ensembl IDs without version. DCE/DCE_COVID_LC_Results_annotated.csv is the edge table
  (458 gene pairs, 416 genes, 399 of them in the 411) and is not the source. Written by `scripts/export_priors.py`.
- `taco_6_gene_signature.txt`: the six gene nested CV signature (provided, HGNC symbols).
- `mrcontrol_32_genes.txt`: the 32 putative causal genes from MRCONTROL (repository SindyPin/Causal-Multiomics-Method):
  16 at alpha = 1 (Mendelian randomisation, SM/SM1) and 16 at alpha = 0 (control theory, README), mapped to Ensembl through
  the HGNC complete set by `scripts/export_priors.py`. Each line is an Ensembl ID followed by `# SYMBOL`.
- `split_consensus_drivers.txt`: SPLIT driver genes supported by two or more of TWMR, CT and DCE.
- `review_complement_9.txt`: the nine complement proteins from the 101 study review (provided, UniProt accessions).

Identifier space must match the matrix. `chorus.features.priors.CausalPriors.to_feature_space` converts with a
two column mapping table (for example Ensembl gene to HGNC symbol from the GENCODE GTF used in quantification).
