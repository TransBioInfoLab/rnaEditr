# Example script to test function TestSingleRegion
# Lanyu Zhang
# 2020-07-21

library(rnaEditr)

data(rnaedit_df)

pheno_df <- readRDS(
  system.file(
    "extdata",
    "pheno_df.RDS",
    package = 'rnaEditr',
    mustWork = TRUE
  )
)
  
######  1. For continuous outcome  ############################################
model_ls <- list(
  modelFormula_char = "age_at_diagnosis ~ rnaEditSummary",
  pheno_df = pheno_df,
  minSize = NULL
) 

TestSingleRegion(
  rnaEdit_num = unlist(rnaedit_df[1,]),
  modelPrep_ls = model_ls,
  respType = "continuous"
)

###  Output  ###
#    estimate   stdErr     pValue
# 1 -95.59037 57.11043 0.09558196



######  2. For binary outcome  ################################################
pheno_df$isHispanic <- dplyr::case_when(
  pheno_df$ethnicity == "HISPANIC OR LATINO" ~ 1L,
  pheno_df$ethnicity == "NOT HISPANIC OR LATINO" ~ 0L,
  TRUE ~ NA_integer_
)
pheno_df$isWhite <- dplyr::case_when(
  pheno_df$race == "WHITE" ~ 1L,
  pheno_df$race == "[Not Available]" ~ NA_integer_,
  TRUE ~ 0L
)


responses_char <- "isWhite"
pheno_df[, responses_char] <- as.factor(
  pheno_df[, responses_char]
)
minSize_num <- CountSamplesPerGroup(
  pheno_df = pheno_df,
  responses_char = "isWhite",
  covariates_char = NULL
)
model_ls <- list(
  modelFormula_char = "isWhite ~ rnaEditSummary",
  pheno_df = pheno_df,
  minSize = minSize_num
) 

TestSingleRegion(
  rnaEdit_num = unlist(rnaedit_df[1,]),
  modelPrep_ls = model_ls,
  respType = "binary"
)

###  Output  ###
#   estimate   stdErr    pValue
# 1 27.71487 17.05264 0.1041084



######  3. For survival outcome  ##############################################
responses_char <- c("OS.time", "OS")

levs <- sum(
  !is.na(unique(pheno_df[, responses_char[2]]))
)

if (levs != 2) {
  warning(
    "Please make sure status indicator has only two values or levels!
           For more details, please see documentation for function Surv().",
    immediate. = TRUE
  )
}

pheno_df$surv_object <- survival::Surv(
  time = pheno_df[, responses_char[1]],
  event = pheno_df[, responses_char[2]]
)

model_ls <- list(
  modelFormula_char = "surv_object ~ rnaEditSummary",
  pheno_df = pheno_df,
  minSize = NULL
) 

TestSingleRegion(
  rnaEdit_num = unlist(rnaedit_df[1,]),
  modelPrep_ls = model_ls,
  respType = "survival"
)

###  Output  ###
#       coef     exp_coef  se_coef    pValue
# 1 -18.2723 1.159949e-08 11.99725 0.1277483
