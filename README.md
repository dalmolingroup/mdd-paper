# Repository for scripts used in *Gene- and transcript-level analyses reveal sex-specific trancriptional alterations in prefrontal cortex in Major Depressive Disorder*. 

## Requisites
    
We use Miniconda (or Anaconda) for environment management (See `environment.yml` file). To create the environment:
`conda env create -f environment.yml`

We use R version 4.1.2 and Bioconductor version 3.14.

All R packages are installed through conda. 

Clone this repository to your local machine and open the project through `mdd.Rproj`.

## Overall project structure

This project is organized in the following directories:

 - `data/`: Holds the RNA-seq processed data, study metadata, and genome references;  
    - `genome/`: auxiliary files for kallisto quantification.  
        - `Homo_sapiens.GRCh38.97.gtf.gz`: GTF file from human gencode/Ensembl version 97.  
        - `Homo_sapiens.GRCh38.cdna.all.fa.gz`: Human transcriptome gencode/Ensembl version 97.  
    - `kallisto/`: output from kallisto quantification mode (`abundance.tsv`). Quantification files are separated by brain region.  
        - `aINS/`:  
        - `Cg25/`:  
        - `dlPFC/`:  
        - `OFC/`:  
        - `Nac/`:  
        - `Sub/`:  
    - `meta/`:  metadata for BioProject SRP115956.   
        - `SraRunTable.txt` : study metadata;
        - `Metadata.csv` :  Information on individual samples. This information was previously available here: http://neuroscience.mssm.edu/nestler/contecenter/listchromatingenedatabase.html. 
    - `vcfs`: VCF files filtered by the variants showed in Table 1 from the paper.

 - `scripts/`: Holds all the scripts used in the analyses;  
 - `results/`: Holds results from each step analysis.  

## Description of scripts used 

- **Quantification**

Script for download, processing and quantification by kallisto;  

- **Exploratory analysis**: 
    - Metadata:  
        - `metadata`: Organizes metadata for further steps.  
    - Transcript and gene estimates:  
        - `tx_gene`: Uses `tximport` to summarise gene counts and prepare data for further steps.  
        - `tx_tx`: Uses `tximport` to prepare data to further steps.  
    - Outlier identification:    
        - `robust_pca`: Performs robust PCA analysis of samples using `rrcov`;   
        - `remove_outliers_samples`: Removes the outliers samples chosen by Robust PCA analysis.   
    - Covariates selection:  
        - `impute_meta`: Imputes data for missing values found in some metadata covariates.  
        - `rank_variables`: Performs covariate analysis.  


- **TAG**:  
    - Feature-wise outlier detection:  
        - `outliers_edge_ppcseq_gene`: Performs identification of outlier genes by `ppcseq`. Outlier genes were removed from DGE analysis.  
        - `outliers_edge_ppcseq_tx`:  Performs identification of outlier transcripts by `ppcseq`. Outlier transcripts were removed from DTE and DTU analyses.  
    - Differential gene expression:  
        - `edger_diff_gene`:  Differential gene expression with `edgeR`.
    - Differential transcript expression:    
        - `edger_diff_tx`: Differential transcript expression with `edgeR`.  
        - `diff_tx_correct`: Performs multiple hypothesis correction with `stageR`.  
    - Differential transcript usage:  
        - `ISA/`: scripts for differential transcript usage using `IsoformSwitchAnalyzeR` are stored in this directory.   
    - Gather results from three methods:  
        - `organize_dge_dte_after_filtering`: Filters the outlier genes and transcripts identified by `ppcseq` from DGE and DTE results.  
        - `summarise_results_dge_dte_dtu`: Removes outlier transcripts from DTU analysis and gathers results from three methods.  

- **Functional analyses**:  
    - `network`: Network inference using stringDB. Visualization by `RedeR` and `ggraph`.  
    - `enrichment`: Enrichment analysis of transcriptionally altered genes using `clusterProfiler`.  
    - `gwas_intersections`: Get genes with genomic regions related to depression using `gwasrapidd`.  
    - `intersection_analysis`: intersection analysis by sex, brain region, and method used.   

- **Additional scripts**:  
    - `plots` and `plot_dtu`: Description of figures produced to the paper.
    - `supp_fig_variants_by_donors`: Description of Supplementary Figures 9 and 10, which represent the presence of depression-associated SNPs on the samples considered. 








    







