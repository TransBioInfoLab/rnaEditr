#' Test associations between phenotype and RNA editing levels.
#'
#' @description A wrapper function to test associations between phenotype and
#'   RNA editing levels in single-site analysis or summarized RNA editing 
#'   levels in region-based analysis.
#' 
#' @param rnaEdit_df A data frame with class \code{rnaEdit_df}, which is a
#'   output from function \code{\link{CreateEditingTable}()} or function
#'   \code{\link{SummarizeAllRegions}()}. This data frame should include RNA
#'   editing level values, with row names as site IDs or region IDs, and column
#'   names as sample IDs.
#' @param pheno_df A data frame with phenotype and covariates, which should
#'   include all the samples in \code{rnaEdit_df}. Please make sure the input
#'   \code{pheno_df} has the variable named \code{"sample"} to indicate sample
#'   IDs.  
#' @param responses_char A character vector of names of response variables in
#'   \code{pheno_df}. When respType is set as \code{"survival"},
#'   \code{responses_char} should have length 2. The first element must be the
#'   name of the variable with following up time, and the second element must 
#'   be status indicator. Status indicator should be coded as 0/1(1=death),
#'   TRUE/FALSE(TRUE=death), or 1/2(death). Please make sure variable names are 
#'   in this order. We have not tested this code on interval-censored data; use 
#'   at your own risk. See \code{\link[survival]{Surv}} for more details.
#' @param covariates_char A character vector of names of covariate variables in
#'   \code{pheno_df}.
#' @param respType Type of outcome. Defaults to \code{"binary"}.
#' @param progressBar Name of the progress bar to use. There are currently five
#'   types of progress bars: \code{"time"}, \code{"none"}, \code{"text"},
#'   \code{"tk"}, and \code{"win"}. Defaults to \code{"time"}. See
#'   \code{\link[plyr]{create_progress_bar}} for more details.
#' @param orderByPval Sort co-edited regions by model p-value or not? Defaults
#'   to TRUE.
#' 
#' @return A data frame with locations of the genomic regions or sites
#'   (\code{seqnames, start, end, width}), test statistics
#'   (\code{estimate, stdErr} or \code{coef, exp_coef, se_coef}), \code{pValue}
#'   and false discovery rate (\code{fdr}).
#' 
#' @importFrom plyr alply
#' 
#' @export
#'
#' @seealso \code{\link{TransformToGR}}, \code{\link{AllCloseByRegions}}, 
#'   \code{\link{AllCoeditedRegions}}, \code{\link{CreateEditingTable}},
#'   \code{\link{SummarizeAllRegions}}, \code{\link{AnnotateResults}}
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
#'   sum_regions <- SummarizeAllRegions(
#'     regions_gr = exm_regions,
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
#'   TestAssociations(
#'     rnaEdit_df = sum_regions,
#'     pheno_df = exm_pheno,
#'     responses_char = "sample_type",
#'     covariates_char = NULL,
#'     respType = "binary"
#'   )
#'
TestAssociations <- function(rnaEdit_df,
                             pheno_df,
                             responses_char,
                             covariates_char = NULL,
                             respType = c("binary", "continuous", "survival"),
                             progressBar = "time",
                             orderByPval = TRUE){

  respType <- match.arg(respType)
  
  if (!("rnaEdit_df" %in% class(rnaEdit_df))) {
    stop("Please use function CreateEditingTable() or SummarizeAllRegions() to 
         create rnaEdit_df.")
  }
  
  metaDataCol_char <- c("seqnames", "start", "end", "width")
  dat <- rnaEdit_df[, !(colnames(rnaEdit_df) %in% metaDataCol_char)]
  
  # Make model formula
  formula_char <- MakeModelFormula(
    responses_char = responses_char,
    covariates_char = covariates_char,
    respType = respType
  )

  # Make model preparation list
  model_ls <- switch(
    
    respType,
    "binary" = {
      
      # To fit GLM with categorical outcomes, must factor character response.
      pheno_df[, responses_char] <- as.factor(
        pheno_df[, responses_char]
      )
      
      # Fit model
      minSize <- CountSamplesPerGroup(
        pheno_df = pheno_df,
        responses_char = responses_char,
        covariates_char = covariates_char
      )
      
      list(
        modelFormula_char = formula_char,
        pheno_df = pheno_df,
        minSize = minSize
      )
      
    },
    
    "continuous" = {
      
      list(
        modelFormula_char = formula_char,
        pheno_df = pheno_df,
        minSize = NULL
      )
      
    },
    
    "survival" = {
      
      # update: add a check to make sure second ele of responses_char only
      #   have two unique values/options/levels. (adter dropping missingness)
      #   make a warning, add link to coxph documentation
      
      # The `[` operator works slightly differently on tibbles, so be aware.
      levs <- sum(
        !is.na(unique(pheno_df[, responses_char[2]]))
      )
      
      if (levs != 2) {
        warning(
          "Please make sure status indicator has only two values or levels!
         For more details, please see documentation for argument 
         'responses_char' by using ?TestAssociations.",
          immediate. = TRUE
        )
      }
      
      
      # Create a survival object, used as a response variable in Cox model
      pheno_df$surv_object <- Surv(
        time = pheno_df[, responses_char[1]],
        event = pheno_df[, responses_char[2]]
      )
      
      list(
        modelFormula_char = formula_char,
        pheno_df = pheno_df,
        minSize = NULL
      )
      
    }
    
  )
  
  # Fit model
  results_ls <- alply(
    .data = dat,
    .margins = 1,
    .fun = function(row){
      TestSingleRegion(
        # we use unlist() instead of as.numeric() here to preserve colnames
        rnaEdit_num = unlist(row), 
        modelPrep_ls = model_ls,
        respType = respType
      )
    },
    .progress = progressBar
  )

  results_df <- do.call(rbind, results_ls)
  results_df$fdr <- p.adjust(results_df$pValue, method = "fdr")
  
  results_meta_df <- cbind(rnaEdit_df[, metaDataCol_char], results_df)
  
  # Order final results
  if (orderByPval) {
    results_meta_df[order(results_meta_df$pValue), ]
  } else {
    results_meta_df
  }
 
}
