context("CreateEditingTable")

data(rnaedit_df)

test_that("CreateEditingTable returns a dataframe with class 'rnaEdit_df'", {
  
  dat_df <- CreateEditingTable(rnaEditMatrix = rnaedit_df)
  
  expect_output(
    print(class(dat_df)),
    '"rnaEdit_df" "data.frame"'
  )
  
})
