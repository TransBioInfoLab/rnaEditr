# Creating gene annotation datasets
# Lanyu Zhang
# Created: Sept 2020
# Last update 09/14/2020

###### 1. Create a function ####################################################
# BiocManager::install("biomaRt")
library(biomaRt)

get_gene_information_biomart <- function(genome = c("hg38","hg19")){
  # check_package("biomaRt")
  genome <- match.arg(genome)
  tries <- 0L
  msg <- character()
  while (tries < 3L) {
    gene.location <- tryCatch({
      host <- ifelse(
        genome == "hg19",
        "grch37.ensembl.org",
        "www.ensembl.org"
      )
      mirror <- list(NULL, "useast", "uswest", "asia")[[tries + 1]]
      ensembl <- tryCatch({
        message(ifelse(is.null(mirror),
                       paste0("Accessing ",
                              host, " to get gene information"),
                       paste0("Accessing ",
                              host, " (mirror ", mirror, ")")))
        biomaRt::useEnsembl(
          "ensembl",
          dataset = "hsapiens_gene_ensembl",
          host = host, mirror = mirror)
      }, error = function(e) {
        message(e)
        return(NULL)
      })
      attributes <- c(
        "ensembl_gene_id",
        "external_gene_name",
        "chromosome_name",
        "strand",
        "end_position",
        "start_position",
        "gene_biotype"
      )
      db.datasets <- biomaRt::listDatasets(ensembl)
      description <- db.datasets[
        db.datasets$dataset == "hsapiens_gene_ensembl", 
      ]$description
      message(
        paste0(
          "Downloading genome information (try:", 
          tries, ") Using: ", description)
      )
      gene.location <- biomaRt::getBM(
        attributes = attributes,
        filters = "chromosome_name",
        values = c(seq_len(22),"X","Y"),
        mart = ensembl
      )
      gene.location
    }, error = function(e) {
      msg <<- conditionMessage(e)
      tries <<- tries + 1L
      NULL
    })
    if (!is.null(gene.location)) break
    if (tries == 3L) stop(
      "failed to get URL after 3 tries:", "\n  error: ", msg
    )
  }
}


###### 2. Create final datasets for hg19 #######################################
get_gene_information_biomart(genome = "hg19")

# Data management
library(dplyr)
hg19_df <- gene.location %>% 
  filter(gene_biotype == "protein_coding") %>%
  dplyr::rename(seqnames = chromosome_name) %>%
  dplyr::rename(start = start_position) %>%
  dplyr::rename(end = end_position) %>%
  dplyr::rename(symbol = external_gene_name) %>%
  dplyr::mutate(seqnames = paste0("chr", seqnames)) %>%
  dplyr::select(-strand, -gene_biotype)
  
# Save as GRanges
# BiocManager::install("GenomicRanges")
library(GenomicRanges)
hg19_gr <- makeGRangesFromDataFrame(
  df = hg19_df,
  keep.extra.columns = TRUE,
  ignore.strand = TRUE
)
saveRDS(hg19_gr, "inst/extdata/hg19_annoGene_gr.RDS")


###### 2. Create final datasets for hg38 #######################################
get_gene_information_biomart(genome = "hg38")

# Data management
library(dplyr)
hg38_df <- gene.location %>% 
  filter(gene_biotype == "protein_coding") %>%
  dplyr::rename(seqnames = chromosome_name) %>%
  dplyr::rename(start = start_position) %>%
  dplyr::rename(end = end_position) %>%
  dplyr::rename(symbol = external_gene_name) %>%
  dplyr::mutate(seqnames = paste0("chr", seqnames)) %>%
  dplyr::select(-strand, -gene_biotype)

# Save as GRanges
# BiocManager::install("GenomicRanges")
library(GenomicRanges)
hg38_gr <- makeGRangesFromDataFrame(
  df = hg38_df,
  keep.extra.columns = TRUE,
  ignore.strand = TRUE
)
saveRDS(hg38_gr, "inst/extdata/hg38_annoGene_gr.RDS")
