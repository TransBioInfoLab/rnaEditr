context("AnnotateResults")

data(rnaedit_df)

# get GRanges for genes
genes_gr <- TransformToGR(
  genes_char = c("PHACTR4", "CCR5", "METTL7A"),
  type = "symbol",
  genome = "hg19"
)

# find close-by regions within the genes
closebyRegions_gr <- AllCloseByRegions(
  regions_gr = genes_gr,
  rnaEditMatrix = rnaedit_df
)

# identify co-edited regions within the genes
coedited_gr <- AllCoeditedRegions(
  regions_gr = closebyRegions_gr,
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

exm_pheno <- readRDS(
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
  pheno_df = exm_pheno,
  responses_char = "sample_type",
  covariates_char = NULL,
  respType = "binary"
)



test_that("AnnotateResults creates correct column names", {
  
  out <- AnnotateResults(
    results_df = results_df,
    closeByRegions_gr = NULL,
    inputRegions_gr = NULL,
    genome = "hg19"
  )
  
  testthat::expect_equal(
    colnames(out),
    c(
      "seqnames", "start", "end", "width", "symbol", "estimate", "stdErr",
      "pValue", "fdr"
    )
  )

})

