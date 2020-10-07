#' Extract clusters of RNA editing sites located closely in genomic 
#'   regions.
#'
#' @description A wrapper function to extract clusters of RNA editing sites
#'   that are located closely in genomic regions.
#'
#' @param regions_gr A GRanges object of input genomic regions.
#' @param rnaEditMatrix A matrix (or data frame) of RNA editing level values on
#'   individual sites, with row names as site IDs in the form of
#'   "chrAA:XXXXXXXX", and column names as sample IDs. Please make sure to
#'   follow the format of example dataset (\code{data(rnaedit_df)}).
#' @param maxGap An integer, genomic locations within \code{maxGap} from each
#'   other are placed into the same cluster. Defaults to 50.
#' @param minSites An integer, minimum number of RNA editing sites within each
#'   resulting cluster. Defaults to 3.
#' @param progressBar Name of the progress bar to use. There are currently five
#'   types of progress bars: \code{"time"}, \code{"none"}, \code{"text"},
#'   \code{"tk"}, and \code{"win"}. Defaults to \code{"time"}. See
#'   \code{\link[plyr]{create_progress_bar}} for more details.
#'
#' @return A GRanges object containing genomic regions of RNA editing sites
#'   located closely within each input pre-defined genomic region.
#'
#' @details The algorithm of this function is based on the
#'   \code{\link[bumphunter]{clusterMaker}} function in the \code{bumphunter}
#'   R package. Each cluster is essentially a group of site locations such that
#'   two consecutive locations in the cluster are separated by less than
#'   \code{maxGap}.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom plyr alply
#'
#' @export
#' 
#' @seealso \code{\link{TransformToGR}}, \code{\link{AllCoeditedRegions}}, 
#'   \code{\link{CreateEditingTable}}, \code{\link{SummarizeAllRegions}}, 
#'   \code{\link{TestAssociations}}, \code{\link{AnnotateResults}}
#'
#' @examples
#'   data(rnaedit_df)
#'   
#'   exm_regions <- TransformToGR(
#'     genes_char = c("PHACTR4", "CCR5", "METTL7A"),
#'     type = "symbol",
#'     genome = "hg19"
#'   )
#'   
#'   AllCloseByRegions(
#'     regions_gr = exm_regions,
#'     rnaEditMatrix = rnaedit_df,
#'     maxGap = 50,
#'     minSites = 3,
#'     progressBar = "time"
#'   )  
#'    
AllCloseByRegions <- function(regions_gr,
                              rnaEditMatrix,
                              maxGap = 50,
                              minSites = 3,
                              progressBar = "time"){
  
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
    
    closeByRegions_ls <- alply(
      .data = unique(hits$queryHits),
      .margins = 1,
      .fun = function(idx){
        SingleCloseByRegion(
          region_df = data.frame(regions_gr[idx]),
          rnaEditMatrix = rnaEditMatrix[unique(hits$subjectHits),],
          maxGap = maxGap,
          minSites = minSites
        )
      },
      # .parallel = parallel,
      .progress = progressBar
    )
    
    # Delete null elements and duplicated elements in the list
    closeByRegions_ls <- unique(
      closeByRegions_ls[lengths(closeByRegions_ls) > 0]
    )
    
    # Turn the list of GRanges to a single GRanges.
    # We suppress warnings because the .Seqinfo.mergexy() function expects that
    #   there will be an overlap between the combined ranges. This will never
    #   be the case for us. The "c" method called here eventually dispatches to
    #   this merge function internally. See below for more information:
    # https://rdrr.io/bioc/GenomeInfoDb/src/R/Seqinfo-class.R
    suppressWarnings(Reduce("c", closeByRegions_ls))
  
}
