context("AllCloseByRegions")

data(rnaedit_df)

genes_gr <- TransformToGR(
  genes_char = c("PHACTR4", "CCR5", "METTL7A"),
  type = "symbol",
  genome = "hg19"
)


test_that("AllCloseByRegions returns a GRanges object", {

  out_gr <- AllCloseByRegions(
    regions_gr = genes_gr,
    rnaEditMatrix = rnaedit_df
  )

  expect_output(str(out_gr), "Formal class 'GRanges'")

})
