# AdaTiSS
AdaTiSS is a R package to calculate tissue specificity scores, analyzed in the paper ["AdaTiSS: A Novel Data-Adaptive Robust Method for Quantifying Tissue Specificity Scores"] (https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab460/6306407?login=true), published in Bioinformatics, 2021.

Please contact Meng Wang by email <mengw1@stanford.edu> for questions. 

## Dependence
* [R](https://www.r-project.org/) (version >= 3.3.0)

## Usage
To illustrate the usage of our package, we took a small set of gene expression in raw TPM from [GTEx project] (https://gtexportal.org/home/datasets) in version 7 as an example.

The data folder includes raw TPM expression for 2000 genes from 181 samples across 32 tissues and a pheonotype info of the tissue type for each sample.

`source('./R/AdaTiSS_fn.R')`

To load expression data

`dat.rna = readRDS(file='./data/raw_tpm_expression.rds')`

To preprocess raw data and obtain filtered gene expression in log scale (for more options, see the file 'AdaTiSS_fn.R')
`X = preproc.filter.fn (dat.rna, dat.type = "TPM or RPKM", proc.zero = 'ceiled to 1', filter.col.prp = 1, exp.thres=1)`

To load phenotype data
`p.dat = read.csv('../data/sample_phenotype.csv')`

To obtain gene expression in tissue level
`tiss.abd = tiss.abd.fn(X, p.dat)`

To apply AdaTiSS (for more options, see the file 'AdaTiSS_fn.R')
`result = AdaTiSS(X, tiss.abd=tiss.abd, adjust=TRUE, adjust.opt=0)`

Output: 

sample normalized scores
`head(result$ada.s)`

tissue specificity scores
`head(result$ada.z)`

population fitting info
`head(result$pop.fit.mx)`

input: 

X --- expression matrix with rows of proteins and columns of samples (prefered in log scale, after preprocessing steps)

tiss.abd --- tissue level abundance (default: NULL)
       
output: 

ada.s --- score for X

ada.z --- score for tiss.abd

pop.fit.mx --- population fitted info (Note: to take another care on the genes with 'pi0.hat' <= 0.5)

`source('./R/AdaTiSS_fn.R')`

`adatiss.result = AdaTiSS(X, tiss.abd)`
