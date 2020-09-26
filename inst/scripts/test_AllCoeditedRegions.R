# Example script to test function AllCoeditedRegions
# Lanyu Zhang
# 2020-07-27

data(rnaedit_df)

input_gr <- TransformToGR(
  genes_char = c("CCR5", "PHACTR4"),
  type = "symbol"
)

AllCoeditedRegions(
  regions_gr = input_gr,
  rnaEditMatrix = rnaedit_df,
  output = "GRanges",
  method = "spearman"
)

# Warning message:
# In .Seqinfo.mergexy(x, y) :
#   The 2 combined objects have no sequence levels in common. (Use
#   suppressWarnings() to suppress this warning.)

# will need to fix function AllCloseByRegions in the same way, because it used
#   the same code.

