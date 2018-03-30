# Author: Jiamao Zheng (jiamaoz@yahoo.com) 03/30/2018

# Adapted from https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl and https://github.com/PMBio/peer/wiki/Tutorial

# usages: Rscript run_peer_args.R /home/jiamaoz/cancer_project/data/mRNA/mRNA_updated/TCGA-OV/TCGA-OV_mRNA.csv.gz TCGA-OV 5 /home/jiamaoz/cancer_project/data/mRNA/mRNA_updated/src

library(data.table)
library(dplyr)
library(tidyr)
library(peer)

argv <- commandArgs(trailingOnly = TRUE) 
input_file_path <- argv[1]
project_name <- argv[2]
PEER_factors_number <- as.numeric(argv[3])
output_dir <- argv[4]

cat('Project: ', project_name, '\n')

# load data  
# The data matrix is assumed to have N rows and G columns, where N is the number of samples, and G is the number of genes.
cat('PEER: loading mRNA data ... \n')
expression <- fread(paste('zcat ', input_file_path, sep=''))  
expression_colnames <- colnames(expression)[2:ncol(expression)]
expression_rownames <- expression[, 1]$aliquot_barcode

expression <- expression[, -1]
expression <- as.matrix(expression)
dimnames(expression) = NULL 
cat('done. \n\n')


# PEER Factors - PEER Factors were calculated using the GTEx pipeline docker
#     container.  The number of PEER factors was chosen according to GTEx v7
#     protocol with the number of samples being the determining factor. If the
#     number of samples was greater than or equal to 350, we used 60 PEER factors.
#     If the number of samples was between 250 and 350, we used 45. Between 150
#     and 250, we used 30, and less than 150 we used 15.

# we choose 60 as number of BRCA individuals is larger than 350 

# https://github.com/PMBio/peer/wiki/Tutorial

# Run PEER 
cat('PEER: runing PEER ... \n')
model = PEER()
invisible(PEER_setPhenoMean(model,expression))
invisible(PEER_setNk(model,PEER_factors_number))   # we choose 60 for BRCA and 45 for OV   
invisible(PEER_setPriorAlpha(model, 0.001, 0.01))
invisible(PEER_setPriorEps(model, 0.1, 10))
invisible(PEER_setNmax_iterations(model, 10000))
time <- system.time(PEER_update(model))


# GET DATA 
cat('PEER: getting factors ...\n')
factors = PEER_getX(model)
c <- paste0("InferredCov",1:ncol(factors))
rownames(factors) = expression_rownames
colnames(factors) = c
cat('factors: samples', ncol(factors), 'PEER factors\n\n')

cat('PEER: getting precision ... \n')
precision = PEER_getAlpha(model)
rownames(precision) = c
colnames(precision) = 'Alpha'
precision <- as.data.frame(precision)
precision$Relevance <- 1.0 / precision$Alpha 
cat('precision: ', dim(precision), '\n\n') 

cat('PEER: getting residuals ... \n')
residuals = t(PEER_getResiduals(model))
rownames(residuals) = expression_colnames
colnames(residuals) = expression_rownames
cat('residuals:', dim(residuals), 'genes x samples ', '\n\n')

# write results
cat("PEER: writing results ... \n")
output_dir <- paste(output_dir, '/', project_name, sep='')
if (!dir.exists(output_dir)){
        dir.create(output_dir)
        cat('\n')
} else {
        cat(output_dir, " Dir already exists!\n\n")
}
write.csv(t(factors), file = gzfile(paste(output_dir, '/', project_name, '_expression_PEER_covariates.csv.gz', sep = '')),  quote = FALSE, row.names = T)
write.csv(precision, file = gzfile(paste(output_dir, '/', project_name, '_expression_PEER_alpha.csv.gz', sep = '')),  quote = FALSE, row.names = T)
write.csv(residuals, file = gzfile(paste(output_dir, '/', project_name, '_expression_PEER_residues.csv.gz', sep = '')),  quote = FALSE, row.names = T)

cat('done!\n\n')





