# Creating example dataset
# Lanyu Zhang
# Created: April 2019
# Last update 09/04/2020

###### 1. Call in Data  ########################################################
dat <- readRDS("inst/extdata/rnaEdit_ge5_df.RDS")

###### 2. Get chr and position info from each row  #############################
library(stringr)

annotations_chr <- data.frame(
  str_split_fixed(row.names(dat), "\\:", 2),
  stringsAsFactors = FALSE
)

colnames (annotations_chr) <- c("chr", "position")

dat_new <- cbind(annotations_chr, dat)

###### 3. Order sites by chr and position  #####################################
data_df <- dat_new[order(dat_new[, "chr"], dat_new[, "position"] ), ]

###### 4. Create example datasets  #############################################

### Example 1
## limit sites that in Genes "PHACTR4", "CCR5", and "METTL7A" (hg19), and then 
##   add few more random sites:
 #   "PHACTR4"    -- chr1:28691093-28826881
 #   "CCR5"       -- chr3:46406633-46417697
 #   "METTL7A"    -- chr12:51313534-51326300
##

PHACTR4_df <- data_df[
  data_df$chr == "chr1" &
    as.numeric(data_df$position) >= 28600000 &
    as.numeric(data_df$position) <= 28827000,
]

CCR5_df <- data_df[
  data_df$chr == "chr3" &
    as.numeric(data_df$position) >= 46000000 &
    as.numeric(data_df$position) <= 48000000,
]

METTL7A_df <- data_df[
  data_df$chr == "chr12" &
    as.numeric(data_df$position) >= 51000000 &
    as.numeric(data_df$position) <= 51400000,
]

rnaedit_df <- rbind(
  PHACTR4_df, CCR5_df, METTL7A_df
)

rnaedit_df <- rnaedit_df[, 3:ncol(rnaedit_df)]

## limit columns to randomly selected 50 samples
set.seed(50)
ranSample <- sample(
  colnames(rnaedit_df),
  size=50
)

rnaedit2_df <- rnaedit_df[, ranSample] #dim: 272 50

### Example 2

## transposed rnaedit_df
library(pathwayPCA)
t_rnaedit_df <- TransposeAssay(
  rnaedit2_df,
  omeNames = "rowNames"
) #dim:50 272

## randomly select 20 sites
set.seed(20)
ranSite <- sample(
  colnames(t_rnaedit_df),
  size=20
)

t_rnaedit_df <- t_rnaedit_df[, ranSite] #dim: 50 20


library(usethis)
usethis::use_data(rnaedit_df, t_rnaedit_df, overwrite = TRUE)

