#' Order RNA editing sites by their genomic locations.
#'
#' @description Split RNA editing sites locations into chromosomes and 
#'   positions, and then order the sites by their genomic locations.
#' 
#' @param sites_char A character vector of RNA editing sites. site IDs should 
#'   be in the form of "chrAA:XXXXXXXX".
#' @param output Type of output data. Defaults to \code{"dataframe"}.
#'
#' @return Depends on the output type. When \code{output} is set as
#'   \code{"vector"}, a character vector of ordered input RNA editing sites
#'   will be returned. When \code{output} is set as \code{"dataframe"}, a
#'   data frame of ordered RNA editing sites with following columns will be
#'   returned:
#'   \itemize{
#'     \item{\code{site} : }{site ID.}
#'     \item{\code{chr} : }{chromosome number.}
#'     \item{\code{pos} : }{genomic location.}
#'   }
#' 
#' @export
#' @keywords internal
#'
#' @examples
#'   exm_sites <- c(
#'     "chr22:41327462", "chr22:24969087",
#'     "chr22:29538891", "chr22:45736763"
#'   )
#'                   
#'   OrderSitesByLocation(
#'     sites_char = exm_sites,
#'     output = "dataframe"
#'   )
#'
OrderSitesByLocation <- function(sites_char, output = c("dataframe","vector")){
  
  output <- match.arg(output)
  
  # Separate site into chromosome and position. Eg.
  #   "chr22:41327462" --> c("chr22", "41327462")

  sites_mat <- do.call(rbind, strsplit(sites_char, split = ":"))
  sites_df <- data.frame(
    site = sites_char,
    chr = sites_mat[, 1],
    pos = as.integer(sites_mat[, 2]),
    stringsAsFactors = FALSE
  )
  
  # Order data
  sites_order_df <- sites_df[order(sites_df$chr, sites_df$pos), ]
  
  # Select and return output ###
  if (output == "dataframe") {
    out <- sites_order_df
  } else {
    out <- sites_order_df$site
  }
  
  out
  
}
