# Example script to test function AddMetaData
# Lanyu Zhang
# created on 2020-08-26
# updated on 2020-08-26

##### Grab sites based on region_df ############################################

# note: range of sites in chr1 in example dataset rnaedit_df:
#   chr1:28661572-chr1:28826284.
data(rnaedit_df)
rnaEditMatrix <- rnaedit_df
region_df <- data.frame(
  seqnames = "chr1",
  start =  28661573,
  end = 28826200,
  stringsAsFactors = FALSE
)

# original code
rnaEditMatrix_sites_df <- OrderSitesByLocation(
  row.names(rnaEditMatrix),
  output = "dataframe"
)

sitesOrdered_df <- rnaEditMatrix_sites_df[
  (rnaEditMatrix_sites_df$chr == region_df$seqnames) &
    (rnaEditMatrix_sites_df$pos >= region_df$start) &
    (rnaEditMatrix_sites_df$pos <= region_df$end),
]

# Tiago's code using granges and subsetByOverlaps
rnaEditMatrix_sites_df <- OrderSitesByLocation(
  row.names(rnaEditMatrix),
  output = "dataframe"
)

sites_gr <- makeGRangesFromDataFrame(
  df = rnaEditMatrix_sites_df,
  start.field = "pos",
  end.field = "pos"
)   

sitesOrdered_gr <- subsetByOverlaps(
  x = sites_gr,
  ranges = makeGRangesFromDataFrame(df = region_df)
)

sitesOrdered2_df <- data.frame(
  "chr" = as.character(seqnames(sitesOrdered_gr)),
  "pos" = start(sitesOrdered_gr)
)

identical(sitesOrdered_df$pos, sitesOrdered2_df$pos)

