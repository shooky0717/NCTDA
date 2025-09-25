# NCTDA
NCTDA is a statistical analysis framework for detecting both SVGs and ctSVGs, which is built on the statistical advances enabling computationally scalable (linear-time and storage) parameter estimation in the spatial covariance functions of the Gaussian Process through the nearest neighbor Gaussian process (NNGP) model.
```R
library(Matrix)
library(ggplot2)
library(reshape2)
library(dplyr)
library(here)
library(BiocParallel)
library(BRISC)
```
## Load in data
```R
data(mOB)
```
## Create NCTDA object
```R
NCTDA <- creatobject(counts = mOB$counts,
                        pos = mOB$pos,
                        prop = mOB$prop,
                        covariates = NULL)
```
## Preprocess the input data
```R
nnCTSVG <- data_preprocess(object = nnCTSVG,
                            gene.threshold = 0.05, spot.threshold = 10,
                            normalized = FALSE)
```
## Build kernel matrix and perform NCTDA test
```R
nnCTSVG<- Test(object =nnCTSVG, X = NULL,
                Cell_types_to_test = NULL, 
                correction = F, 
                pv.adjust = "BY",
                cov.model = "exponential",
                n_threads = 1)
```
## Obtain the ctSVG list
```R
ct_genes_list <- lapply(NCTDA@cell_types, function(i) {
  df <- NCTDA@Test[[i]]
  rows_use <- df[df$padj < 0.05 , ]
  rownames(rows_use)
})
names(ct_genes_list) <- NCTDA@cell_types
ctSVG_list <- unique(unlist(ct_genes_list))
```
