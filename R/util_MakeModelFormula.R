#' Make model formula.
#'
#' @description Make model formula for different types of phenotype responses.
#' 
#' @param responses_char A character vector of the response variable.
#' @param covariates_char A character vector of the covariate variables.
#' @param respType Type of outcome. Defaults to \code{"binary"}.
#' 
#' @return A character vector of the model formula.
#' 
#' @details When \code{respType} is set as \code{"survival"},
#'   \code{"surv_object"} is only a placeholder here , which will be defined 
#'   later in \code{TestSingleRegion()}.
#' 
#' @export
#' @keywords internal
#'
#' @examples
#' 
#' MakeModelFormula(
#'   responses_char = "age",
#'   covariates_char = c("sex", "tumor_type"),
#'   respType = "continuous"
#' )
#' 
#' MakeModelFormula(
#'   responses_char = "sample_type",
#'   covariates_char = c("sex", "tumor_type"),
#'   respType = "binary"
#' )
#' 
#' MakeModelFormula(
#'   responses_char = c("OS.time", "OS"),
#'   covariates_char = c("sex", "tumor_type"),
#'   respType = "survival"
#' )
#'
MakeModelFormula <- function(responses_char,
                             covariates_char = NULL,
                             respType = c("binary", "continuous", "survival")){
  
  respType <- match.arg(respType)
  
  if (respType == "survival") {
    responses_char <- "surv_object"
  }
  
  mainMod_char <- paste(responses_char, "rnaEditSummary", sep = " ~ ")
  
  if (is.null(covariates_char)) {
    mainMod_char
  } else {
    
    covMod_char <- paste(covariates_char, collapse = " + ")
    paste(mainMod_char, covMod_char, sep = " + ")
    
  }
  
}
