#' Add Annotations to site-specific or region-based analysis results.
#'
#' @description Add annotations to site-specific or region-based analysis 
#'   results from function \code{\link{TestAssociations}}.
#'
#' @param results_df An output data frame from function 
#'   \code{TestAssociations}, which includes variables for locations and result
#'   of statistical tests for the genomic sites or regions. 
#' @param closeByRegions_gr An output GRanges object from function
#'   \code{AllCloseByRegions}, defaults to \code{NULL}.
#' @param inputRegions_gr A GRanges object for input genomic
#'   regions, defaults to \code{NULL}.
#' @param genome Use \code{"hg19"} or \code{"hg38"} gene reference. Defaults 
#'   to \code{"hg38"}.
#' @param analysis Results type. Defaults to \code{"region-based"}. When it's
#'   set to \code{"site-specific"}, arguments \code{closeByRegions_gr} and 
#'   \code{inputRegions_gr} will not be used and set to NULL automatically.
#'
#' @return A data frame with locations of the genomic sites or regions
#'   (\code{seqnames, start, end, width}), annotations for locations
#'   (\code{inputRegion, closeByRegion, symbol}), test statistics
#'   (\code{estimate, stdErr} or \code{coef, exp_coef, se_coef}), \code{pValue}
#'   and false discovery rate (\code{fdr}).
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' 
#' @export
#'
#' @seealso \code{\link{TransformToGR}}, \code{\link{AllCloseByRegions}}, 
#'   \code{\link{AllCoeditedRegions}}, \code{\link{CreateEditingTable}},
#'   \code{\link{SummarizeAllRegions}}, \code{\link{TestAssociations}}
#'   
#' @examples
#'   data(rnaedit_df)
#'   
#'   # get GRanges for genes
#'   genes_gr <- TransformToGR(
#'     genes_char = c("PHACTR4", "CCR5", "METTL7A"),
#'     type = "symbol",
#'     genome = "hg19"
#'   )
#'   
#'   # find close-by regions within the genes
#'   closebyRegions_gr <- AllCloseByRegions(
#'     regions_gr = genes_gr,
#'     rnaEditMatrix = rnaedit_df
#'   )
#'   
#'   # identify co-edited regions within the genes 
#'   coedited_gr <- AllCoeditedRegions(
#'     regions_gr = closebyRegions_gr,
#'     rnaEditMatrix = rnaedit_df,
#'     output = "GRanges",
#'     method = "spearman"
#'   )
#'   
#'   # summarize editing levels within each gene by maximum
#'   summarizedRegions_df <- SummarizeAllRegions(
#'     regions_gr = coedited_gr,
#'     rnaEditMatrix = rnaedit_df,
#'     selectMethod = MaxSites
#'   )
#'   
#'   exm_pheno <- readRDS(
#'     system.file(
#'     "extdata",
#'     "pheno_df.RDS",
#'     package = 'rnaEditr',
#'     mustWork = TRUE
#'     )
#'   )
#'   
#'   # test summarized editing levels against survival outcome
#'   results_df <- TestAssociations(
#'     rnaEdit_df = summarizedRegions_df,
#'     pheno_df = exm_pheno,
#'     responses_char = "sample_type",
#'     covariates_char = NULL,
#'     respType = "binary"
#'   )
#'   
#'   AnnotateResults(
#'     results_df = results_df,
#'     closeByRegions_gr = closebyRegions_gr,
#'     inputRegions_gr = genes_gr,
#'     genome = "hg19"
#'   )
#'   
AnnotateResults <- function(results_df,
                            closeByRegions_gr = NULL,
                            inputRegions_gr = NULL,
                            genome = c("hg38", "hg19"),
                            analysis = c("region-based", "site-specific")){

    genome <- match.arg(genome)
    analysis <- match.arg(analysis)
    
    results_gr <- makeGRangesFromDataFrame(
      df = results_df,
      keep.extra.columns = TRUE
    )
    
    analysisSite_logi <- analysis == "site-specific"
    nullInput_logi <- is.null(inputRegions_gr)
    nullClose_logi <- is.null(closeByRegions_gr)
    
    
    ###  Add Annotations  ###
    # Add inputRegion annotation if available
    if (analysisSite_logi | nullInput_logi) {
      addInput_gr <- results_gr
    } else {
      
      addInput_gr <- AddMetaData(
        target_gr = results_gr,
        annot_gr = inputRegions_gr,
        annotType_char = "region",
        genome = genome,
        annotLabel_char = "inputRegion"
      )
      
    }
    
    # Add closeByRegion annotation if available
    if (analysisSite_logi | nullClose_logi) {
      addMetaData_gr <- addInput_gr
    } else {
      
      addMetaData_gr <- AddMetaData(
        target_gr = addInput_gr,
        annot_gr = closeByRegions_gr,
        annotType_char = "region",
        genome = genome,
        annotLabel_char = "closeByRegion"
      )
      
    }
    
    # Add symbol annotation
    coeditedAnno_gr <- AddMetaData(
      target_gr = addMetaData_gr,
      annot_gr = NULL,
      annotType_char = "geneSymbol",
      genome = genome,
      annotLabel_char = "symbol"
    )
    
    
    ###  Wrangle Column Names  ###
    # Order column names of final dataset
    colPre_char <- c("estimate", "stdErr")
    
    # maybe turn this statement into a switch function once there are more
    #   different types of outcome estimates later
    if (!(all(colPre_char %in% colnames(results_df)))) {
      colPre_char <- c("coef", "exp_coef", "se_coef")
    }
    
    # Add on symbols
    colPre_char <- c("symbol", colPre_char, "pValue", "fdr")
    # Because the column order, if applicable, is "inputRegion" then
    #   "closeByRegion", we want to left-concatenate them in reverse order
    if (!analysisSite_logi) {
      
      if (!nullClose_logi) {
        colPre_char <- c("closeByRegion", colPre_char)
      }
      if (!nullInput_logi) {
        colPre_char <- c("inputRegion", colPre_char)
      }
      
    }
    colOrder_char <- c("seqnames", "start", "end", "width", colPre_char)
    
    # Organize final annotation dataframe
    dat <- data.frame(coeditedAnno_gr)
    dat$seqnames <- as.character(dat$seqnames)
    
    dat[ ,colOrder_char]
  
}
