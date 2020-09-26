# Example script to test function AddMetaData
# Lanyu Zhang
# created on 2020-08-12
# updated on 2020-08-12

data(rnaedit_df)

# get GRanges for genes
genes_gr <- TransformToGR(
  genes_char = c("PHACTR4", "CCR5", "METTL7A"),
  type = "symbol"
)

# find close-by regions within the genes
closeby_gr <- AllCloseByRegions(
  regions_gr = genes_gr,
  rnaEditMatrix = rnaedit_df
)

# identify co-edited regions within the genes
coedited_gr <- AllCoeditedRegions(
  regions_gr = genes_gr,
  rnaEditMatrix = rnaedit_df,
  output = "GRanges",
  method = "spearman"
)

# summarize editing levels within each gene by maximum
summarizedRegions_df <- SummarizeAllRegions(
  regions_gr = coedited_gr,
  rnaEditMatrix = rnaedit_df,
  selectMethod = MaxSites
)  

pheno_df<- readRDS(
  system.file(
    "extdata",
    "pheno_df.RDS",
    package = 'rnaEditr',
    mustWork = TRUE
  )
)  

# test summarized editing levels against survival outcome
results_df <- TestAssociations(
  rnaEdit_df = summarizedRegions_df,
  pheno_df = pheno_df,
  responses_char = "sample_type",
  covariates_char = NULL,
  respType = "binary"
)  

AnnotateResults(
  results_df = results_df,
  closeByRegions_gr = closeby_gr,
  inputRegions_gr = genes_gr
) 
