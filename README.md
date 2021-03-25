# AdaTiSS
AdaTiSS is a R package to calculate tissue specificity scores, analyzed in the paper "AdaTiSS: A Novel Data-Adaptive Robust Method for Quantifying Tissue Specificity Scores".

Please contact Meng Wang by email <mengw1@stanford.edu> for questions. 

## Dependence
* [R](https://www.r-project.org/) (version >= 3.3.0)

## Usage
input: 

X --- expression matrix with rows of proteins and columns of samples (prefered in log scale, after preprocessing steps)

tiss.abd --- tissue level abundance (default: NULL)
       
output: 

ada.s --- score for X

ada.z --- score for tiss.abd

pop.fit.mx --- population fitted info (Note: to take another care on the genes with 'pi0.hat' <= 0.5)

`source('./R/AdaTiSS_fn.R')`

`adatiss.result = AdaTiSS(X, tiss.abd)`
