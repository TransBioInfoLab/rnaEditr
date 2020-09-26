context("TransformToGR")


test_that("TransformToGR gives correct errors", {
  expect_error(
    TransformToGR(genes_char = c("Hello"), type = "symbol"),
    "No gene found."
  )

  expect_error(
    TransformToGR(genes_char = c("GENE1", "GENE2"), type = "symbol"),
    "No gene found."
  )

})

# test_that("TransformToGR gives correct warnings", {
#   expect_message(
#     TransformToGR(
#       genes_char = c("ABCA4", "GENE1", "PLD6"),
#       type = "symbol"
#     ),
#     "1 gene(s) not found. These gene(s) are:"
#   )
# 
#   expect_message(
#     TransformToGR(
#       genes_char = c("GENE2", "PLD6", "GENE1", "ABCA4"),
#       type = "symbol"
#     ),
#     " 2 gene(s) not found. These gene(s) are:"
#   )
# 
# })
