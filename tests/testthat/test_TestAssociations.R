context("TestAssociations")

data(rnaedit_df)

genes_gr <- TransformToGR(
  genes_char = c("PHACTR4", "CCR5", "METTL7A"),
  type = "symbol",
  genome = "hg19"
)

exm_regions <- AllCoeditedRegions(
  regions_gr = genes_gr,
  rnaEditMatrix = rnaedit_df,
  output = "GRanges",
  method = "spearman"
)

sum_regions <- SummarizeAllRegions(
  regions_gr = exm_regions,
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



test_that("TestAssociations returns a dataframe", {
  
  dat_df <- TestAssociations(
    rnaEdit_df = sum_regions,
    pheno_df = exm_pheno,
    responses_char = "sample_type",
    covariates_char = NULL,
    respType = "binary"
  )
  
  expect_output(
    str(dat_df),
    "'data.frame'"
  )
  
})
