# Phytoplankton community responses to 20th century eutrophication of Swiss Plateau lakes: a proxy comparison from short sediment cores
## Antonia Klatt, Daniel B. Nelson, Jan Waelchli, Theresa Wietelmann, Nathalie Dubois & S. Nemiah Ladd

[Link]

In this study, we analyzed algal proxies in short sediment cores from Swiss Plateau lakes (Greifensee (GRE18, GRE 24), Murtensee (MUR23) and Zugersee (ZUG_N, ZUG_S)) to trace changes in phytoplankton community composition in response to the strong eutrophication during the 20th century. We used the abundance of phytol and phytosterols (phytol:sterol index, PSI) and compound-specific hydrogen isotope ratios of C16:0 fatty acid and phytol (δ2HC16:0 Acid/Phytol) to reconstruct relative proportions of cyanobacteria and eukaryotic algae. We then compared the lipid-based proxies with plastid and cyanobacterial 23S rRNA gene amplicon sequence variants (ASVs) and cross-validated these sedimentary proxies with published observational data from Buergi et al. (2003) (Aquat. Ecosyst. Health Manag. 6(2),147-158) and Merz et al. (2023) (https://doi.org/10.1038/s41558-023-01615-6) and other historical observations. 

#### This repository contains all datasets from this study and respective R scripts for analysis and visualization (except from raw DNA sequencing data).
#### Note: Raw DNA sequencing reads are uploaded at the European Nucleotide Archive (ENA) repository [Link].

### 01_Data
This folder contains lipid and sedDNA datasets from Greifensee, Murtensee and Zugersee. 

The sedDNA data included here have been processed by the dada2 pipeline incl. taxonomic assignments by the Basic Local Alignment Search Tool (BLAST) against the National Center for Biotechnology Information (NCBI) database (accessed on 20th November 2025) and filtered for phytoplankton groups while reads from non-algal groups have been excluded. Absolute algal ASV counts from each sample were normalized by total sum scaling. 

Raw sequencing reads from sedDNA samples can be found elsewhere (ENA repository [Link]).

R codes for running the dada2 pipeline and filtering for phytoplankton groups will be uploaded delayed in this repository or elsewhere...

Note that for Greifensee, lipid data are based on sediment cores collected in 2018, while sedDNA data are from a newly collected core from 2024. For Zugersee, data from two different sediment cores are included - one core was taken at the northern basin (ZUG_N), one core was taken at the southern basin (ZUG_S).

### 02_R_scripts
This folder contains R codes for data analysis and plotting.
