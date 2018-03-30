# peer_factor_analysis

This R script is used to run PEER factor analysis on gene expression 

## Prerequisite
- [R 3.0+](https://www.r-project.org)
- [data.table](https://github.com/Rdatatable/data.table): `install.packages('data.table')`
- [dplyr](https://github.com/tidyverse/dplyr): `install.packages('dplyr')`
- [tidyr](http://tidyr.tidyverse.org): `install.packages('tidyr)`
- [peer](https://github.com/PMBio/peer/wiki/Installation-instructions): `R CMD INSTALL R_peer_source_1.3.tgz`

## Input
- The data is assumed to have N rows and G columns, where N is the number of samples, and G is the number of genes.

## Run
```
	Rscript run_peer_args.R <input file path> <project name> <number of peer factors> <output directory> 

	for example: Rscript run_peer_args.R /home/jiamaoz/cancer_project/data/mRNA/mRNA_updated/TCGA-OV/TCGA-OV_mRNA.csv.gz TCGA-OV 5 /home/jiamaoz/cancer_project/data/mRNA/mRNA_updated/
```

## outputs 
- {project_name}_expression_PEER_residuals.csv.gz
- {project_name}_expression_PEER_alpha.csv.gz
- {project_name}_expression_PEER_covariates.csv.gz
