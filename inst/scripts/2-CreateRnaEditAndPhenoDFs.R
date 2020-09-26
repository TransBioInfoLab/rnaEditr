# Create RNA Editing level data set and Phenotype information data set
# & Creating phenotype dataset
# Lanyu Zhang
# Created: Feb 2019
# Last update 8/13/2020

###### 1. Data management for rna editing level matrix #########################
library(tidyverse)

### 1.1 Call in dataset

# I'm using read_delim() instead of read.delim() here because it's faster.
brca <- read_delim(
  file = "inst/extdata/BRCA_AG_annovar_nomutation.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Delete the last empty column and convert tibble into dataframe because we
#   need to add row names in the next step.
brca_df <- brca %>% select(-X944) %>% as.data.frame()


### 1.2 Create row names

# Row names will be the combination of the first two parts in column "Sample",
#   e.g.: in "chr1|109795329|Synonymous|CELSR2|+|no-repetitive|no-conserve",
#   we only need "chr1" and "109795329", and the row name will be
#   "chr1:109795329".
brca_final_df <- brca_df %>%
  separate(
    col = Samples,
    into = c("chr", "position"),
    sep = "([|])",
    extra = "drop"
  )

rownames(brca_final_df) <- paste(
  brca_final_df$chr, brca_final_df$position, sep = ":"
)


### 1.3 Exclude sites on Chr Y then drop off "chr" and "position" columns
brca_final_df <- brca_final_df %>%
  rownames_to_column() %>%
  filter(chr != "chrY") %>%
  select(-chr, -position) %>%
  column_to_rownames()


### 1.4 Drop off samples with over 50% missing sites, then drop off sites with
###   over 1/3 missing samples.

# why 1/3? here is the link to the reference paper
#   https://www.nature.com/articles/nature24041

missingSites <- apply(
  X = brca_final_df,
  MARGIN = 2,
  FUN = function(x) mean(is.na(x))
)
missingSites50 <- missingSites[missingSites < 0.5]

brca_final_df <- brca_final_df[
  ,match(names(missingSites50), colnames(brca_final_df))
]

missingSamples <- apply(
  X = brca_final_df,
  MARGIN = 1,
  FUN = function(x) mean(is.na(x))
)

missingSamples33 <- missingSamples[missingSamples < 1/3]

brca_final_df <- brca_final_df[
  match(names(missingSamples33), row.names(brca_final_df)),
]




###### 2. Impute rna editing level matrix ######################################

### 2.1 Transpose datasets for imputation purpose

# we are using pathwayPCA::TransposeAssay() here instead of t() to preserve row
#   and column names when transposing the data frame.

# BiocManager::install("pathwayPCA")
library(pathwayPCA)

brca_final_t_df <- TransposeAssay(
  assay_df = brca_final_df,
  omeNames = "rowNames"
)


### 2.2 Register a parallel backend and impute the dataset
library(parallel)
library(doParallel)

cl <- makeCluster(detectCores()/2)
registerDoParallel(cl)


# Default setting for number of trees in missForest is 100. As the OOB errors
#   were similar from using 100, 500, 1000 trees, we decided to use default 
#   number of trees in the package
library(missForest)

brca_final_t_imp_df <- missForest(
  xmis = brca_final_t_df,
  verbose = TRUE,
  parallelize = "variables"
)

stopCluster(cl)

brca_final_imp_df <- TransposeAssay(
  assay_df = brca_final_t_imp_df$ximp,
  omeNames = "rowNames"
)




###### 3. Create phenotype data ################################################

library(TCGAbiolinks)
library(dplyr)

### 3.1 Call in datasets
brca_final_imp_df <- readRDS(paste0(data, "/rnaEdit_df.RDS"))

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
pheno_df <- clinical.BCRtab.all$clinical_patient_brca

survival <- read_tsv(
  file = paste0(data.raw, "/TCGA-BRCA.survival.tsv")
)
survival_spec <- spec(survival)
survival_df <- read_tsv(
  file = paste0(data.raw, "/TCGA-BRCA.survival.tsv"),
  col_types = survival_spec
)

### 3.2 Get submitter id and sample type info from column names of
###   brca_final_imp_df
sample <- data.frame(
  sample = colnames(brca_final_imp_df),
  stringsAsFactors = FALSE
)

info_df <- tidyr::separate(
  data = sample,
  col = sample,
  into = c("BRCA", "sample_type", "TCGA", "tissue_source_site", "patient_id")
)

sample_df <- data.frame(
  sample = sample$sample,
  submitter_id = paste(
    info_df$TCGA, info_df$tissue_source_site, info_df$patient_id, sep = "-"),
  sample_type = ifelse(info_df$sample_type == "Normal",
                       "Solid Tissue Normal", "Primary Tumor"),
  stringsAsFactors = FALSE
)

### 3.3 Select needed variables and merge datasets to get final phenotype 
###   dataset

### select needed variables in pheno_df
pheno_df <- pheno_df %>%
  select(
    bcr_patient_barcode, gender, race, ethnicity,
    age_at_diagnosis, er_status_by_ihc
  )

### merge sample_df with pheno_df
sample_pheno_df <- merge(
  sample_df, pheno_df,
  by.x = "submitter_id",
  by.y = "bcr_patient_barcode",
  sort = FALSE
)

### select needed variables in survival_df
survival_df <- survival_df[,c("_PATIENT", "OS", "OS.time")]
survival_df <- unique(survival_df)

### merge sample_pheno_df with survival_df
sample_pheno_surv_df <- merge(
  sample_pheno_df, survival_df,
  by.x = "submitter_id",
  by.y = "_PATIENT",
  sort = FALSE,
  all.x = TRUE
)

sample_pheno_surv_df$age_at_diagnosis <- as.numeric(
  sample_pheno_surv_df$age_at_diagnosis
)




###### 4. Data management for rna editing matrix and phenotype dataset #########

### 4.1 Delete males from final rna editing matrix and phenotype dataset
sample_pheno_surv_df <- sample_pheno_surv_df %>%
  filter(gender == "FEMALE")
brca_final_imp_df <- brca_final_imp_df[, sample_pheno_surv_df$sample]

### 4.2 Re-categorize variables in sample_pheno_surv_df
sample_pheno_surv_df$race <- ifelse(
  sample_pheno_surv_df$race == "[Not Available]",
  NA,
  sample_pheno_surv_df$race
)

sample_pheno_surv_df$ethnicity <- ifelse(
  sample_pheno_surv_df$ethnicity == "[Not Available]",
  NA,
  sample_pheno_surv_df$ethnicity
)

sample_pheno_surv_df$er_status_by_ihc <- ifelse(
  sample_pheno_surv_df$er_status_by_ihc == "[Not Evaluated]",
  NA,
  sample_pheno_surv_df$er_status_by_ihc
)

### 4.3. Limit brca_final_imp_df to sites with >= 5% rna editing levels
rnaEditFreq <- data.frame(
  site = row.names(brca_final_imp_df),
  freq = apply(brca_final_imp_df, 1, function(x) mean(x != 0)),
  stringsAsFactors = FALSE
)

rnaEditFreq <- rnaEditFreq[rnaEditFreq$freq >= 0.05,]

rnaEdit_ge5_df <- brca_final_imp_df[
  row.names(brca_final_imp_df) %in% rnaEditFreq$site,
]

# saveRDS(rnaEdit_ge5_df, "inst/extdata/rnaEdit_ge5_df.RDS")
# saveRDS(sample_pheno_surv_df, paste0("inst/extdata/pheno_df.RDS"))
