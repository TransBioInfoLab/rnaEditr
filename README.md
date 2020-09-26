## rnaEditr: an R package for RNA-seq data analysis

## Description
`rnaEditr` is an R package that identifies genomic sites and genomic regions 
that are differentially edited in RNA-seq datasets. `rnaEditr` can analyze 
studies with continuous, binary, or survival phenotypes, along with multiple 
covariates and/or interaction effects. To identify hyper-edited regions, 
`rnaEditr` first determines co-edited sub-regions without using any phenotype 
information. Next, `rnaEditr` tests association between RNA editing levels 
within the co-edited regions with binary, continuous or survival phenotypes.

## Installation

The latest version can be installed by

```
devtools::install_github("TransBioInfoLab/rnaEditr")
```
After installation, the rnaEditr package can be loaded into R using:

```{r}
library(rnaEditr)
```

