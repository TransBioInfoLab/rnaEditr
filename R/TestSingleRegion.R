#' Test associations between phenotype and RNA editing levels.
#'
#' @description Test associations between phenotype and
#'   RNA editing levels in a single site or summarized RNA editing levels
#'   in a single region.
#'
#' @param rnaEdit_num A named numeric vector of (summarized) RNA editing level
#'   values with sample IDs as names.
#' @param modelPrep_ls A list includes \code{modelFormula_char} which is 
#'   created by function \code{\link{MakeModelFormula}}, \code{pheno_df} which 
#'   is the input phenotype data frame in \code{\link{TestAssociations}}, and
#'   \code{minSize} (minimum sample size per group to use regular logistic 
#'   regression) which is created by function
#'   \code{\link{CountSamplesPerGroup}} when \code{respType} is \code{"binary"}.
#' @param respType Type of outcome. Defaults to \code{"binary"}. 
#' 
#' @return a dataframe with test statistics (\code{estimate, stdErr, pValue} or
#'   \code{coef, exp_coef, se_coef, pValue}).
#' 
#' @details \code{minSize} is used by function \code{TestSingleRegion} to 
#'   decide on whether to use regular logistic regression or Firth corrected 
#'   logistic regression (\url{"https://www.jstor.org/stable/2336755"}).  
#'
#' @importFrom stats as.formula binomial coef glm lm p.adjust drop1 vcov
#' @importFrom stats gaussian 
#' @importFrom survival Surv coxph
#' @importFrom logistf logistf 
#' 
#' @export
#' @keywords internal
#'
#' @examples
#'   data(rnaedit_df)
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
#'   exm_model <- list(
#'     modelFormula_char = "age_at_diagnosis ~ rnaEditSummary",
#'     pheno_df = exm_pheno,
#'     minSize = NULL
#'   )
#'   
#'   TestSingleRegion(
#'     rnaEdit_num = unlist(rnaedit_df[2,]),
#'     modelPrep_ls = exm_model,
#'     respType = "continuous"
#'   )
#'
TestSingleRegion <- function(rnaEdit_num,
                             modelPrep_ls,
                             respType = c("binary", "continuous", "survival")){
  
  respType <- match.arg(respType)
  
  # Convert rnaEdit_num to rnaEditOne_df with sample names to be merged
  #   into final df later.
  rnaEditOne_df <- data.frame(
    sample = names(rnaEdit_num),
    rnaEditSummary = rnaEdit_num,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # Merge pheno_df and rnaEditOne_df
  rnaEditOnePheno_df <- merge(
    x = modelPrep_ls$pheno_df,
    y = rnaEditOne_df,
    by = "sample"
  )
  
  # Fit model
  switch(
    respType,
    "binary" = {
      
      if (is.null(modelPrep_ls$minSize)) {
        stop(
"minSize can't be NULL. Please use function CountSamplesPerGroup to find it!")
      } else if (modelPrep_ls$minSize < 10) {
        
        .LOGISTFTest(
          modelFormula_char = modelPrep_ls$modelFormula_char,
          rnaEditOnePheno_df = rnaEditOnePheno_df
        )
        
      } else {
        
        .GLMTest(
          modelFormula_char = modelPrep_ls$modelFormula_char,
          rnaEditOnePheno_df = rnaEditOnePheno_df,
          family = binomial(link = "logit")
        )
        
      }
      
    },
    
    "continuous" = {
      
      .GLMTest(
        modelFormula_char = modelPrep_ls$modelFormula_char,
        rnaEditOnePheno_df = rnaEditOnePheno_df,
        family = gaussian(link = "identity")
      )
    
    },
    
    "survival" = {
      
      .CoxPHTest(
        modelFormula_char = modelPrep_ls$modelFormula_char,
        rnaEditOnePheno_df = rnaEditOnePheno_df
      )
      
    }
  )
}


# Firth Bias-Corrected Logistic
.LOGISTFTest <- function(modelFormula_char, rnaEditOnePheno_df){
  
  f <- tryCatch({
    logistf(
      formula = as.formula(modelFormula_char),
      data = rnaEditOnePheno_df
    )
  }, warning = function(w){
    NULL
  }, error = function(e){
    NULL
  })
  
  if (is.null(f)) {
    
    result_df <- data.frame(
      estimate = NA_real_,
      stdErr = NA_real_,
      pValue = 1
    )
    
  } else {
    
    result_df <- data.frame(
      estimate = coef(f)["rnaEditSummary"],
      stdErr = sqrt(diag(vcov(f)))["rnaEditSummary"],
      pValue = drop1(f)["rnaEditSummary", "P-value"],
      stringsAsFactors = FALSE
    )
    rownames(result_df) <- NULL
    
  }
  
  result_df
  
}


# Ordinary Least Squares / Logistic
.GLMTest <- function(modelFormula_char,
                     rnaEditOnePheno_df,
                     family){
  
  f <- tryCatch({
    glm(
      formula = as.formula(modelFormula_char),
      family = family,
      data = rnaEditOnePheno_df
    )
  }, warning = function(w){
    NULL
  }, error = function(e){
    NULL
  })
  
  if (is.null(f)) {
    
    result_df <- data.frame(
      estimate = NA_real_,
      stdErr = NA_real_,
      pValue = 1
    )
    
  } else {
    
    result_df <- data.frame(
      estimate = coef(summary(f))["rnaEditSummary", "Estimate"],
      stdErr = coef(summary(f))["rnaEditSummary", "Std. Error"],
      # Regardless of the statistical test, the p-value is in column four
      pValue = coef(summary(f))["rnaEditSummary", 4]
    )
    rownames(result_df) <- NULL
    
  }
  
  result_df
  
}


# CoxPH
.CoxPHTest <- function(modelFormula_char, rnaEditOnePheno_df){
  
  f <- tryCatch({
    coxph(
      formula = as.formula(modelFormula_char),
      data = rnaEditOnePheno_df
    )
  }, warning = function(w){
    NULL
  }, error = function(e){
    NULL
  })
  
  if (is.null(f)) {
    
    result_df <- data.frame(
      coef = NA_real_,
      exp_coef = NA_real_,
      se_coef = NA_real_,
      pValue = 1
    )
    
  } else {
    
    result_df <- data.frame(
      coef = coef(summary(f))["rnaEditSummary", "coef"],
      exp_coef = coef(summary(f))["rnaEditSummary", "exp(coef)"],
      se_coef = coef(summary(f))["rnaEditSummary", "se(coef)"],
      pValue = coef(summary(f))["rnaEditSummary", "Pr(>|z|)"]
    )
    rownames(result_df) <- NULL
    
  }
  
  result_df
  
}

