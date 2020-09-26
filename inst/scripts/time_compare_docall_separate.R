# Time comparison of tidyr::separate() with do.call()
# Lanyu Zhang
# 2019-08-24

library(rnaEditr)

##### Call in datasets
inputRegions_gr <- readRDS(
  system.file(
    "extdata",
    "annoGene_gr.RDS",
    package = 'rnaEditr',
    mustWork = TRUE
  )
)

rnaEdit_ge5_df <- readRDS(
  "~/Dropbox (BBSR)/Rpackage-rnaEditr/BRCA_analysis/data_created/rnaEdit_ge5_df.RDS"
)
dim(rnaEdit_ge5_df)



##### Run AllCloseByRegions using tidyr::separate() in inner function 
#####   OrderSitesByLocation
a <- Sys.time()
closeByRegions_gr <- AllCloseByRegions(
  regions_gr = inputRegions_gr, #length:23349
  rnaEditMatrix = rnaEdit_ge5_df, #dim:33806 221
  maxGap = 50,
  minSites = 3,
  cores = 1
)
Sys.time() - a
# tidyr::separate: 28.20043 mins for 33806 x 221


##### Run AllCloseByRegions using do.call() in inner function 
#####   OrderSitesByLocation
OrderSitesByLocation <- function(sites_char, output = c("dataframe","vector")){
  
  output <- match.arg(output)
  
  ### Separate site into chromosome and position. Eg.
  ###   "chr22:41327462" --> c("chr22", "41327462")
  
  # Original script
  sites_mat <- do.call(rbind, strsplit(sites_char, split = ":"))
  sites_df <- data.frame(
    site = sites_char,
    chr = sites_mat[, 1],
    pos = as.integer(sites_mat[, 2]),
    stringsAsFactors = FALSE
  )
  
  # sites_df <- data.frame("sites" = sites_char) %>% 
  #   # tidyr::separate()
  #   separate(
  #     col = "sites",
  #     into = c("chr","pos"),
  #     convert = TRUE,
  #     remove = FALSE
  #   )
  
  ### Order data
  sites_order_df <- sites_df[order(sites_df$chr, sites_df$pos), ]
  
  ### Select and return output ###
  if (output == "dataframe") {
    out <- sites_order_df
  } else {
    out <- sites_order_df$site
  }
  
  out
  
}

b <- Sys.time()
closeByRegions2_gr <- AllCloseByRegions(
  regions_gr = inputRegions_gr, #length: 23349
  rnaEditMatrix = rnaEdit_ge5_df, #dim: 33806 221
  maxGap = 50,
  minSites = 3,
  cores = 1
)
Sys.time() - b
# do.call(): 26.8293 mins for 33806 x 221


# Note from Gabriel: it appears that the time difference between tidyr::
#   separate() and base:: do.call() is negligible. I then recommend that any
#   uses of the Reduce() function  be replaced with do.call, except in the two
#   cases that Reduce() is used on a list of GRanges objects in conjunction with
#   the c() function.

