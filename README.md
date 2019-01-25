
This is an R package that implments a new method named BaSiC to discern intratumor heterogeneity (ITH) using genomic data (e.g., exome-seq or whole genome sequencing data) from both bulk tumor samples as well as single cells. BaSiC stands for **B**ulk tumor **a**nd **Si**ngle **C**ell. 

The main goal is to cluster somatic mutations based on their Variant Allele Frequencies (VAFs) from bulk tumor samples as well as somatic mutation calls from a set of single cells. 

An important novelty of our method is to based on observations that the percentage of missing values of somatic mutation calls from single cells are highly depend on read-depth information. Therefore we incorporate read-depth information to recover some missed somatic mutation calls. 

Two examples are using our package are provided in two R markddown files ex_Wang_2014.Rmd and ex_Gawad_2014.Rmd. 
