#' Find minimum sample Size per group.
#'
#' @description Find minimum sample size for each group which is decided by the
#'   combination of variables with class character or factor.
#' 
#' @param pheno_df A data frame with phenotype and covariates, which should
#'   include all the samples in \code{rnaEdit_df}. Please make sure the input
#'   \code{pheno_df} has the variable named \code{"sample"} to indicate sample
#'   IDs.
#' @param responses_char A character vector of names of response variables in
#'   \code{pheno_df}. When respType is set as \code{"survival"},
#'   \code{responses_char} should have length 2. The first element must be the
#'   name of the variable with follow up time, and the second element must be
#'   the status indicator. Status indicator should be coded as 0/1(1=death),
#'   TRUE/FALSE(TRUE=death), or 1/2(2=death). Please make sure variable names 
#'   are in this order. This code has not been tested for interval-censored 
#'   data yet.
#' @param covariates_char A character vector of names of covariate variables in
#'   \code{pheno_df}.
#' 
#' @return An integer.
#'   
#' @importFrom stats complete.cases
#' 
#' @export
#' @keywords internal
#'
#' @examples
#'   exm_pheno <- readRDS(
#'     system.file(
#'     "extdata",
#'     "pheno_df.RDS",
#'     package = 'rnaEditr',
#'     mustWork = TRUE
#'     )
#'   )
#'   
#'   CountSamplesPerGroup(
#'     pheno_df = exm_pheno,
#'     responses_char = "sample_type",
#'     covariates_char = "race"
#'   )
#'
CountSamplesPerGroup <- function(pheno_df,
                                 responses_char,
                                 covariates_char){
  
    # Limit df to variables in responses_char  or covariates with char or 
    #   factor class.
    cov_df <- pheno_df[
      , colnames(pheno_df) %in% covariates_char, drop = FALSE
    ]
    
    select_char <- vapply(cov_df, class, character(1))
    select_names <- names(
      select_char[select_char %in% c("character", "factor", "logical")]
    )
    
    # Limit final table to rows without missingness: because if there's
    #   missingness then model would ignore incomplete data.
    table_df <- pheno_df[, c(responses_char, select_names), drop = FALSE]
    table_df <- table_df[complete.cases(table_df), , drop = FALSE]
    
    # allCombos_char <- apply(
    #   X = table_df,
    #   MARGIN = 1,
    #   FUN = function(row) { paste(row, collapse = "_") }
    # )
    
    # Find the minimum sample size in each group.
    # min(table(allCombos_char))
    min(table(table_df))
  
}
