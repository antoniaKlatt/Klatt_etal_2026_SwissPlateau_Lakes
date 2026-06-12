## Microbiota Analysis

**00_data:** Used data are archived in ENA (PRJEB114126) 

**01_rm_primers:** cutadapt was used to remove primers 

**02_dadapipe:** dada2 was used to filter and trim sequences, learn and infer the error model, remove chimeras and count ASVs 

**03_NCBI_taxonomy:** taxonomy was assigned with the NCBI database 

**04_filtering:** data were filtered by sample sequencing depth (390'000 sequences), normalised by TSS and outliers removed. Then ASVs were filtered by phytoplankton classes.

**05_OUTPUT:** phytoplankton data use for downstream analysis

**0_tree** phylogenetic tree of phytoplankton
