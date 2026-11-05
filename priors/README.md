# Causal priors

Place the published gene lists here (one identifier per line, `#` comments allowed):

- `taco_411_dce_genes.txt`: the 411 DCE genes from TACO (repository SindyPin/TACO, file DCE/DCE_COVID_LC_Results_annotated.csv; export the gene column as Ensembl IDs without version).
- `taco_6_gene_signature.txt`: the six gene nested CV signature (provided, HGNC symbols).
- `mrcontrol_32_genes.txt`: the 32 putative causal genes from MRCONTROL (repository SindyPin/Causal-Multiomics-Method).
- `split_consensus_drivers.txt`: SPLIT driver genes supported by two or more of TWMR, CT and DCE.
- `review_complement_9.txt`: the nine complement proteins from the 101 study review (provided, UniProt accessions).

Identifier space must match the matrix. `chorus.features.priors.CausalPriors.to_feature_space` converts with a
two column mapping table (for example Ensembl gene to HGNC symbol from the GENCODE GTF used in quantification).
