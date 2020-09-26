context("SummarizeAllRegions")

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

test_that("SummarizeAllRegions returns a dataframe with class 'rnaEdit_df'", {
  
  dat_df <- SummarizeAllRegions(
    regions_gr = exm_regions,
    rnaEditMatrix = rnaedit_df
  )
  
  expect_output(
    print(class(dat_df)),
    '"rnaEdit_df" "data.frame"'
  )
  
})
