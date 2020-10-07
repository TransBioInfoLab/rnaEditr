#' Extract RNA editing sites located in a genomic region.
#' 
#' @description Extract and order RNA editing sites located within an input 
#'   genomic region.
#' 
#' @param region_df A data frame with the input genomic region. Please make 
#'   sure columns \code{seqnames}, \code{start}, and \code{end} are included in 
#'   the data frame.
#' @param rnaEditMatrix A matrix (or data frame) of RNA editing level values on
#'   individual sites, with row names as site IDs in the form of
#'   "chrAA:XXXXXXXX", and column names as sample IDs. Please make sure to
#'   follow the format of example dataset (\code{data(rnaedit_df)}).
#' @param output Type of output data. Defaults to \code{"locationsOnly"}.
#'
#' @return When \code{output} is set to \code{"locationsOnly"}, a data frame of 
#'   extracted and ordered RNA editing sites with columns \code{chr} and 
#'   \code{pos} will be returned. 
#'   
#'   When \code{output} is set to \code{"locationsAndValues"}, a data frame of
#'   RNA editing level values from the extracted and ordered sites will be
#'   returned. Please note that site IDs will be in row names of the output
#'   data frame.
#' 
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start
#' 
#' @export
#' @keywords internal
#'
#' @examples
#'   data (rnaedit_df)
#'   
#'   exm_region <- data.frame(
#'      seqnames = "chr1",
#'      start =  28000000,
#'      end = 28826881, 
#'      stringsAsFactors = FALSE
#'   )
#'   
#'   GetSitesLocations(
#'      region_df = exm_region,
#'      rnaEditMatrix = rnaedit_df,
#'      output = "locationsOnly"
#'   )[1:3, ]
#'    
GetSitesLocations <- function(region_df,
                              rnaEditMatrix,
                              output = c("locationsOnly", "locationsAndValues")
){
  
    output <- match.arg(output)
    
    rnaEditMatrix_sites_df <- OrderSitesByLocation(
      sites_char = row.names(rnaEditMatrix),
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
    
    sitesOrdered_df <- data.frame(
      "chr" = as.character(seqnames(sitesOrdered_gr)),
      "pos" = start(sitesOrdered_gr)
    )
    
    if (nrow(sitesOrdered_df) == 0) {
      return(NULL)
    }
    
    if (output == "locationsOnly"){
      out <- sitesOrdered_df
    } else {
      
      # Extract rna editing cluster by ordered input sites
      out <- rnaEditMatrix[
        match(
          paste0(sitesOrdered_df$chr,":",sitesOrdered_df$pos),
          row.names(rnaEditMatrix)
        ),
      ]
      
    }
    
    out
  
}
