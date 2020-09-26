#' Transform gene symbols or region ranges into GRanges object.
#' 
#' @description Transform a character vector of gene symbols or region ranges
#'   into a GRanges object.
#' 
#' @param genes_char A character vector of gene symbols or region ranges. If 
#'   you select \code{type} to be \code{"symbol"}, then please make sure your 
#'   input of \code{genes_char} is in the format of c("ABCB10", "PEX26"). If 
#'   you select \code{type} to be \code{"region"}, then please make sure your 
#'   input of \code{genes_char} is in the format of
#'   c("chr1:33772367-33791699", "chr22:18555686-18573797").
#' @param type What is the type of \code{genes_char}. Can be \code{"symbol"} 
#'   (default) or \code{"region"}.
#' @param genome Use \code{"hg19"} or \code{"hg38"} gene reference. Defaults 
#'   to \code{"hg38"}. It's only used when \code{type} is set to 
#'   \code{"symbol"}
#' 
#' @return A GRanges object with \code{seqnames}, \code{ranges} and
#'   \code{strand}.
#'  
#' @details \code{TransformToGR()} uses the hg19/hg38 genes to associate gene 
#'   symbols with their genomic region ranges. The pre-processed dataset is 
#'   saved in inst/extdata in this package.
#'   
#'   Users who wish to add gene symbols to the GRanges created using
#'   function \code{TransformToGR()} can use function \code{AddMetaData()}. 
#'   Please see \code{\link{AddMetaData}} for details.
#'  
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' 
#' @export
#' 
#' @seealso \code{\link{AllCloseByRegions}}, \code{\link{AllCoeditedRegions}}, 
#'   \code{\link{CreateEditingTable}}, \code{\link{SummarizeAllRegions}}, 
#'   \code{\link{TestAssociations}}, \code{\link{AnnotateResults}}
#'
#' @examples
#'   TransformToGR(
#'     genes_char = c("PHACTR4", "CCR5", "METTL7A"),
#'     type = "symbol",
#'     genome = "hg19"
#'   )
#'   
#'   TransformToGR(
#'     genes_char = c("chr22:18555686-18573797", "chr22:36883233-36908148"),
#'     type = "region",
#'     genome = "hg19"
#'   )
#'  
TransformToGR <- function(genes_char,
                          type = c("symbol", "region"),
                          genome = c("hg38", "hg19")){
  
  genome <- match.arg(genome)
  type <- match.arg(type)

  if (type == "symbol") {
    
    if (genome == "hg38") {
      genes_gr <- readRDS(
        system.file(
          "extdata",
          "hg38_annoGene_gr.RDS",
          package = 'rnaEditr',
          mustWork = TRUE
        )
      )
    } else {
      genes_gr <- readRDS(
        system.file(
          "extdata",
          "hg19_annoGene_gr.RDS",
          package = 'rnaEditr',
          mustWork = TRUE
        )
      )
    }
    
    genes_df <- data.frame(genes_gr)
    genes_df$seqnames <- as.character(genes_df$seqnames)
    
    gene_sub_df <- genes_df[genes_df$symbol %in% genes_char, ]
    
    if (nrow(gene_sub_df) == 0){
      
      stop("No gene found.")
      
    } else {
      
      geneNoMiss <- genes_char[genes_char %in% gene_sub_df$symbol]
      geneMiss <- genes_char[!(genes_char %in% gene_sub_df$symbol)]
      
      if (length(geneMiss) > 0) {
        
        message(
          sprintf("%i gene(s) not found. These gene(s) are:", length(geneMiss))
        )
        print(geneMiss)
        
      }
      
      gene_sub_df$symbol <- factor(
        x = gene_sub_df$symbol,
        levels = geneNoMiss,
        ordered = TRUE
      )
      gene_sub_df <- gene_sub_df[
        order(gene_sub_df$symbol), 
      ]
      gene_sub_df$symbol <- as.character(
        gene_sub_df$symbol
      )
      
    }
    
    # Then turn the data frame of ranges (seqnames,start,end) into GRanges
    dat_gr <- GRanges(
      seqnames = gene_sub_df$seqnames,
      ranges = IRanges(
        gene_sub_df$start,
        gene_sub_df$end
      ),
      symbols = gene_sub_df$symbol
    )
    
  } else {
    
    # Split genes_char into seqname, start and end. For example,
    #   "chr22:18555686-18573797" --> "chr22", "18555686", "18573797"
    site_mat <- do.call(rbind, strsplit(genes_char, split = c("\\:|\\-")))
    
    gene_sub_df <- data.frame(
      seqnames = site_mat[, 1],
      start = as.integer(site_mat[, 2]),
      end = as.integer(site_mat[, 3]),
      stringsAsFactors = FALSE
    )
    
    # Then turn the data frame of ranges (seqnames,start,end) into GRanges
    dat_gr <- GRanges(
      seqnames = gene_sub_df$seqnames,
      ranges = IRanges(
        gene_sub_df$start,
        gene_sub_df$end
      )
    )
    
  }
  
  dat_gr
  
}
