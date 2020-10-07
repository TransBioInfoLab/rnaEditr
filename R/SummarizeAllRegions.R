#' Summarize RNA editing levels from multiple sites in regions.
#'
#' @description A wrapper function to summarize RNA editing levels from 
#'   multiple sites in regions.
#'
#' @param regions_gr A GRanges object of input genomic regions.
#' @param rnaEditMatrix A matrix (or data frame) of RNA editing level values 
#'   for individual sites, with row names as site IDs in the form of
#'   "chrAA:XXXXXXXX", and column names as sample IDs. Please make sure to
#'   follow the format of example dataset (\code{data(rnaedit_df)}).
#' @param selectMethod Method for summarizing regions. Available options are 
#'   \code{"MaxSites", "MeanSites", "MedianSites", "PC1Sites"}. Please see
#'   \code{\link{RegionSummaryMethod}} for more details.
#' @param ... Dots for additional internal arguments (currently unused).
#' @param progressBar Name of the progress bar to use. There are currently five
#'   types of progress bars: \code{"time"}, \code{"none"}, \code{"text"},
#'   \code{"tk"}, and \code{"win"}. Defaults to \code{"time"}. See
#'   \code{\link[plyr]{create_progress_bar}} for more details.
#'
#' @return A data frame of the class \code{rnaEdit_df}, includes
#'   variables \code{seqnames, start, end, width} and summarized RNA editing
#'   levels in each sample.
#'   
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom plyr adply
#' 
#' @export
#'
#' @seealso \code{\link{TransformToGR}}, \code{\link{AllCloseByRegions}}, 
#'   \code{\link{AllCoeditedRegions}}, \code{\link{CreateEditingTable}}, 
#'   \code{\link{TestAssociations}}, \code{\link{AnnotateResults}}
#'   
#' @examples
#'   data(rnaedit_df)
#'   
#'   genes_gr <- TransformToGR(
#'     genes_char = c("PHACTR4", "CCR5", "METTL7A"),
#'     type = "symbol",
#'     genome = "hg19"
#'   )
#'   
#'   exm_regions <- AllCoeditedRegions(
#'     regions_gr = genes_gr,
#'     rnaEditMatrix = rnaedit_df,
#'     output = "GRanges",
#'     method = "spearman"
#'   )
#'   
#'   SummarizeAllRegions(
#'     regions_gr = exm_regions,
#'     rnaEditMatrix = rnaedit_df
#'   )[1:3, 1:6]
#'
SummarizeAllRegions <- function(regions_gr,
                                rnaEditMatrix,
                                selectMethod = MedianSites,
                                progressBar = "time", ...){

    # parallel <- register_cores(cores)
    
    sites_mat <- do.call(
      rbind, strsplit(row.names(rnaEditMatrix), split = ":")
    )
    sites_df <- data.frame(
      site = row.names(rnaEditMatrix),
      chr = sites_mat[, 1],
      pos = as.integer(sites_mat[, 2]),
      stringsAsFactors = FALSE
    )
    
    sites_gr <- makeGRangesFromDataFrame(
      df = sites_df,
      start.field = "pos",
      end.field = "pos"
    )
    
    hits <- data.frame(
      findOverlaps(
        query = regions_gr,
        subject = sites_gr
      )
    )
    
    regions_df <- data.frame(regions_gr)
    regions_df$seqnames <- as.character(regions_df$seqnames)
    
    summarizedRegions_df <- adply(
      .data = regions_df,
      .margins = 1,
      .fun = function(row){
        SummarizeSingleRegion(
          region_df = row,
          rnaEditMatrix = rnaEditMatrix[unique(hits$subjectHits),],
          selectMethod = selectMethod,
          ...
        )
      },
      # .parallel = parallel,
      .progress = progressBar
    )
  
    keepCol_lgl <- colnames(summarizedRegions_df) != "strand"
    sum_df <- summarizedRegions_df[, keepCol_lgl]
    
    class(sum_df) <- c("rnaEdit_df", "data.frame")
    
    sum_df

}
