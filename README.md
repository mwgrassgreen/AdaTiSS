# AdaTiSS
input: X --- expression matrix (prefered in log scale, after preprocessing steps)
       tiss.abd --- tissue level abundance
output: ada.s --- score for X
        ada.z --- score for tiss.abd
        pop.fit.mx --- population fitted info (Note: to take another care on the genes with 'pi0.hat' <= 0.5)

source('AdaTiSS_mw_20210128.R')
adatiss.result = AdaTiSS(X, tiss.abd=NULL)
