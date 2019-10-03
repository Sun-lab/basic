
## BaSiC

This is an R package that implments a new method named BaSiC to discern intratumor heterogeneity (ITH) using genomic data (e.g., exome-seq or whole genome sequencing data) from both bulk tumor samples as well as single cells. BaSiC stands for **B**ulk tumor **a**nd **Si**ngle **C**ell. 

The main goal is to cluster somatic mutations based on their Mutant Allele Frequencies (MAFs, which is often referred to as Variant Allele Frequencies (VAF)) from bulk tumor samples as well as somatic mutation calls from a set of single cells. 

An important novelty of our method is to based on observations that the percentage of missing values of somatic mutation calls from single cells are highly depend on read-depth information. Therefore we incorporate read-depth information to recover some missed somatic mutation calls. 

Here is the file structure of this repository
```
.
├── basic
├── basic_0.99.002.tar.gz
├── data
│   ├── B-SCITE_SCFile_n_1000_m_20.txt
│   ├── B-SCITE_SCFile_n_100_m_20.txt
│   ├── B-SCITE_SCFile_n_200_m_20.txt
│   ├── B-SCITE_bulkFile_n_100.txt
│   ├── B-SCITE_bulkFile_n_1000.txt
│   ├── B-SCITE_bulkFile_n_200.txt
│   ├── Final_Mutation_Calls_with_counts.tsv
│   ├── PhyloWGS_cnv_data.txt
│   ├── PhyloWGS_ssm_data.txt
│   ├── basic_simulation_100.matrices
│   ├── basic_simulation_100.samples
│   ├── basic_simulation_100.scores
│   ├── basic_simulation_100_ml0.gv
│   ├── basic_simulation_100_ml0.newick
│   ├── pat_1_cluster_dat
│   ├── pat_2_cluster_dat
│   ├── pat_3_cluster_dat
│   ├── pat_4_cluster_dat
│   ├── pat_5_cluster_dat
│   └── pat_6_cluster_dat
├── ex_Gawad_2014.Rmd
├── ex_Gawad_2014.html
├── ex_Wang_2014.Rmd
├── ex_Wang_2014.html
├── ex_figures
│   ├── ex_simu_PhyloWGS_n_subclone.png
│   └── ex_simu_PhyloWGS_tree.png
├── ex_simulation.Rmd
├── ex_simulation.html
├── ex_simulation_comparison.Rmd
├── ex_simulation_comparison.html
├── ex_simulation_comparison.pdf
└── ex_simulation_sensitivity
    ├── ex_simulation_sensitivity1.Rmd
    ├── ex_simulation_sensitivity1.html
    ├── ex_simulation_sensitivity2.Rmd
    ├── ex_simulation_sensitivity2.html
    ├── ex_simulation_sensitivity3.Rmd
    ├── ex_simulation_sensitivity3.html
    ├── ex_simulation_sensitivity4.Rmd
    ├── ex_simulation_sensitivity4.html
    ├── ex_simulation_sensitivity_cn1.Rmd
    ├── ex_simulation_sensitivity_cn1.html
    ├── ex_simulation_sensitivity_cn2.Rmd
    └── ex_simulation_sensitivity_cn2.html
```

* basic: the R package

* data: data used or generated for various analysis. 

  + Final_Mutation_Calls_with_counts.tsv is mutation calls from bulk tumor for the anlaysis of ex_Gawad_2014. 
  
  + pat_*_cluster_dat is the single cell mutation data for ex_Gawad_2014. 
  
* ex_simulation
  
  + use simulated data to evaluate Basic

* ex_simulation_comparison

    + use simulated data to evaluate two alternative methods: B-SCITE and PhyloWGS. 

* ex_simulation_sensitivity

  + ex_simulation_sensitivity1-4.Rmd: sensitivity analysis by chaning the read-depth cutoff to call muations. The default is a muation is called if the read-depth equal to or larger than 20. Here we evaluated four other cutoffs of 16, 18, 22, and 24
  
  + ex_simulation_sensitivity_cn1-2.Rmd: sensitivity analysis after adding noise in the copy number calls. 
  
