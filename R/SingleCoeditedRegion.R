#' Extracts contiguous co-edited genomic regions from a single 
#'   genomic region.
#' 
#' @description Extracts contiguous co-edited genomic regions from an input
#'   genomic region.
#' 
#' @param region_df A data frame with the input genomic region. Please make 
#'   sure columns \code{seqnames}, \code{start}, and \code{end} are included in
#'   the data frame.
#' @param rnaEditMatrix A matrix (or data frame) of RNA editing level values on
#'   individual sites, with row names as site IDs in the form of
#'   "chrAA:XXXXXXXX", and column names as sample IDs. Please make sure to
#'   follow the format of example dataset (\code{data(rnaedit_df)}).
#' @param output Type of output data, can be "GRanges" or "dataframe". Defaults
#'   to "GRanges".
#' @param rDropThresh_num Threshold for minimum correlation between RNA editing
#'   levels of one site and the mean RNA editing levels of the rest of the 
#'   sites. Please set a number between 0 and 1. Defaults to 0.4.
#' @param minPairCorr Minimum pairwise correlation coefficient of a cluster
#'   is used as a filter to select clusters for output. Only clusters with all 
#'   pairwise correlations between sites more than \code{minPairCorr} will be 
#'   selected for output. 
#'   To use this filter, set this argument to a number between -1 and 1
#'   (defaults to 0.1). To turn it off, please set the argument to -1.
#' @param minSites Minimum number of sites to be considered a region. Only
#'   regions with more than \code{minSites} number of sites will be returned.
#' @param method Method for computing correlations. Defaults to 
#'   \code{"spearman"}.
#' @param minEditFreq Threshold for minimum percentage of samples for a given 
#'   site. The \code{r_drop} value of the sites with frequency lower than 
#'   \code{minEditFreq} will be set as NA. Please set a number between 0 and 1. 
#'   Defaults to 0.05.
#' @param returnAllSites When no co-edited region is found in
#'   an input genomic region, \code{returnAllSites = TRUE} indicates
#'   outputting all the sites from the input region, while
#'   \code{returnAllSites = FALSE} indicates not returning any site from the
#'   input region. Defaults to FALSE.
#' @param verbose Should messages and warnings be displayed? Defaults to TRUE,
#'   but is set to FALSE when called from within \code{AllCoeditedRegions()}.
#'
#' @return When \code{output} is set to
#'   \code{"GRanges"}, a GRanges object with \code{seqnames}, \code{ranges} and
#'   \code{strand} of the contiguous co-edited regions will be returned. 
#'   
#'   When \code{output} is set to \code{"dataframe"}, a data frame with
#'   following columns will be returned:
#'   \itemize{
#'     \item{\code{site} : }{site ID.}
#'     \item{\code{chr} : }{chromosome.}
#'     \item{\code{pos} : }{genomic location.}
#'     \item{\code{r_drop} : }{the correlation between RNA editing levels of 
#'     one site and the mean RNA editing levels of the rest of the sites.}
#'     \item{\code{keep} : }{indicator for co-edited sites, The sites with
#'     \code{keep = 1} belong to the contiguous and co-edited region.}
#'     \item{\code{keep_contiguous} : }{contiguous co-edited region number.}
#'     \item{\code{regionMinPairwiseCor} : }{the minimum pairwise correlation 
#'     of a co-edited region.}
#'     \item{\code{keep_regionMinPairwiseCor} : }{equals 1 for contiguous
#'     co-edited subregions. The regions with \code{keepminPairwiseCor = 1}
#'     are the ones that passed the \code{regionMinPairwiseCor} filter and will 
#'     be returned as a co-edited sub-region.}
#'   }
#' 
#' @export
#' @keywords internal
#'
#' @examples
#'   data(rnaedit_df)
#'   
#'   exm_region <- data.frame(
#'     seqnames = "chr1",
#'     start =  28691093,
#'     end = 28826881, 
#'     stringsAsFactors = FALSE
#'   )
#'   
#'   SingleCoeditedRegion(
#'     region_df = exm_region,
#'     rnaEditMatrix = rnaedit_df,
#'     minPairCorr = 0.25,
#'     output = "dataframe",
#'     method = "spearman"
#'   )
#'    
SingleCoeditedRegion <- function(region_df,
                                 rnaEditMatrix,
                                 output = c("GRanges", "dataframe"),
                                 rDropThresh_num = 0.4,
                                 minPairCorr = 0.1,
                                 minSites = 3,
                                 method = c("spearman", "pearson"),
                                 minEditFreq = 0.05,
                                 returnAllSites = FALSE,
                                 verbose = TRUE){
  
    output <- match.arg(output)
    method <- match.arg(method)
    
    
    # Get the rnaEditing Values for the Region
    rnaEditCluster_df <- GetSitesLocations(
      region_df = region_df,
      rnaEditMatrix = rnaEditMatrix,
      output = "locationsAndValues"
    )
    
    if (anyNA(rnaEditCluster_df)) {
      stop("Complete data is required for rnaEditMatrix.")
    }
    
    if (is.null(rnaEditCluster_df)) {
      return(NULL)
    } else if (nrow(rnaEditCluster_df) < minSites) {
      return(NULL)
    }
    
    
    # Transpose RNA editing level matrix
    rnaEditClusterT_mat <- t(rnaEditCluster_df)
    
    # Mark co-edited sites
    keepSites_df <- MarkCoeditedSites(
      rnaEditCluster_mat = rnaEditClusterT_mat,
      rDropThresh_num = rDropThresh_num,
      method = method,
      minEditFreq = minEditFreq,
      verbose = verbose
    )
    
    # Find contiguous and co-edited subregions
    keepContiguousSites_df <- FindCorrelatedRegions(
      sites_df = keepSites_df,
      featureType = "site",
      minSites_int = 3
    )
    
    # Split sites by subregion
    keepContiguousSites_ls <- split(
      x = keepContiguousSites_df$site,
      f = keepContiguousSites_df$subregion
    )
    
    # if no co-edited region is found, then not going over the next step with
    #   function GetMinPairwiseCor()
    if (max(keepContiguousSites_df$subregion) == 0) {
      minPairCorr <- -1
    }
    
    # Calculate min pairwise correlation and subset for each subregion and add
    #   information to df
    pairwiseResults <- GetMinPairwiseCor(
      rnaEditCluster_mat = rnaEditClusterT_mat,
      minPairCorr = minPairCorr,
      probes_ls = keepContiguousSites_ls,
      method = method
    )
    
    # Create output 
    if (output == "GRanges") {
      
      # Create output GRanges
      out <- SitesToRegion(
        sitesSubregion_df = keepContiguousSites_df,
        sitesAreOrdered = TRUE,
        keepminPairwiseCor_df = pairwiseResults$keepminPairwiseCor_df,
        returnAllSites = returnAllSites,
        verbose = verbose
      )
      
    } else {
      
      # Create output data frame
      out <- CreateOutputDF(
        keepSites_df = keepSites_df,
        keepContiguousSites_df = keepContiguousSites_df,
        keepminPairwiseCor_df = pairwiseResults$keepminPairwiseCor_df,
        returnAllSites = returnAllSites,
        verbose = verbose
      )
      
    }
    
    # Return
    out
  
}

