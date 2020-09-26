#' Calculates R-drop values for RNA editing sites.
#' 
#' @description Calculates the correlation coefficient between RNA editing 
#'   levels of one site and the mean RNA editing levels of the rest of the 
#'   sites in a region.
#'
#' @param data A data frame of RNA editing level values on individual sites,
#'   with row names as sample IDs and column names as site IDs in the form of
#'   "chrAA:XXXXXXXX".
#' @param method Method for computing correlation. Defaults to 
#'   \code{"spearman"}.
#' @param minEditFreq Threshold for minimum percentage of edited samples 
#'  for a given 
#'   site. The \code{r_drop} value of the sites with frequency lower than 
#'   \code{minEditFreq} will be set as NA. Please set a number between 0 and 1. 
#'   Defaults to 0.05.
#' @param verbose Should messages and warnings be displayed? Defaults to TRUE.
#' 
#' 
#' @importFrom stats cor
#' 
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item{\code{site} : }{site ID.}
#'     \item{\code{r_drop} : }{the correlation between RNA editing levels of 
#'     one site and the mean RNA editing levels of the rest of the sites.}
#'   }
#' 
#' @export
#' @keywords internal
#'
#' @examples
#'   data(t_rnaedit_df)
#'   
#'   ordered_cols <- OrderSitesByLocation(
#'     sites_char = colnames(t_rnaedit_df),
#'     output = "vector"
#'   )
#'   exm_data <- t_rnaedit_df[, ordered_cols]
#'   
#'   CreateRdrop(
#'     data = exm_data,
#'     method = "spearman"
#'   )
#'    
CreateRdrop <- function(data,
                        method = c("spearman", "pearson"),
                        minEditFreq = 0.05,
                        verbose = TRUE){
  
  method <- match.arg(method)
  col_data <- colnames(data)
  
  out_num <- vapply(
    X = seq_along(col_data),
    FUN = function(column){
      
      # Remove site i and then compute row mean
      data_no_i <- data[, -column]
      data_i <- data[, column]
      
      # Calculate row mean for data_no_i
      data_no_i_mean <- rowMeans(data_no_i, na.rm = TRUE)
      
      # Check the frequence of zero rna editing levels in data_i
      data_i_freq <- mean(data_i != 0, na.rm = TRUE)
      
      if (data_i_freq < minEditFreq) {
        r_drop <- NA_real_
      } else {
        
        # Correlate dat_i and dat_no_i_mean
        r_drop <- cor(
          x = data_i,
          y = data_no_i_mean,
          method = method,
          use = "pairwise.complete.obs"
        )
        
      }
      
      r_drop
      
    },
    FUN.VALUE = numeric(1)
  )
  
  if (verbose & anyNA(out_num)) {
    outMessage <-
      "Some of the sites have fewer than 5% samples with RNA editing levels.
       The r_drop values for these sites have been set to NA."
    warning(outMessage)
  }
  
  data.frame(
    site = col_data,
    r_drop = out_num,
    stringsAsFactors = FALSE
  )
  
}


