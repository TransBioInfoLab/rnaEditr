#' Transposed breast cancer example dataset.
#'
#' @description A subset of the TCGA breast cancer RNA editing dataset for 20
#'   randomly selected RNA editing sites and 50 randomly selected subjects from
#'   example dataset \code{rnaedit_df}. Please note that this is only a
#'   computational testing dataset for inner functions of this package. To test
#'   main functions, please use dataset \code{rnaedit_df} instead.
#'
#' @format A data frame containing RNA editing levels for 50 subjects (in the
#'   rows) at 20 edited sites (in the columns). Row names are sample IDs and 
#'   column names are site IDs.
#'
#' @source Synapse database ID: syn2374375.
#'
"t_rnaedit_df"
