#' Add metadata columns to GRanges object.
#' 
#' @description Add metadata information to GRanges object.
#' 
#' @param target_gr A GRanges object that will be annotated with metadata
#'   
#' @param annot_gr A GRanges object that includes the metadata
#'   information. When \code{annotType_char} = \code{"geneSymbol"}, this 
#'   argument can be left as NULL, and the gene annotation file saved in the 
#'   package will be used to annotate \code{target_gr}. When 
#'   \code{annotType_char} = \code{"region"}, this argument must be specified, 
#'   each row in \code{target_gr} will be annotated with rows in 
#'   \code{annot_gr} that overlap with it.  
#' @param annotType_char Type of the metadata column, defaults to
#'   \code{"geneSymbol"}.
#' @param genome Use \code{"hg19"} or \code{"hg38"} gene reference. Defaults 
#'   to \code{"hg38"}.
#' @param annotLabel_char Name of the metadata column, defaults to 
#'   \code{"symbol"} which corresponds to default setting \code{"geneSymbol"} 
#'   for argument \code{annotType_char}.
#'  
#' @return A GRanges object with \code{seqnames, ranges, region}, and supplied
#'   metadata information.
#' 
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits mcols "mcols<-"
#' 
#' @export
#' @keywords internal
#'   
#' @examples
#'   data(rnaedit_df)
#'   
#'   input_gr <- TransformToGR(
#'     genes_char = "PHACTR4",
#'     type = "symbol",
#'     genome = "hg19"
#'   )
#'  
#'   # identifies co-edited region within input_gr 
#'   coedited_gr <- AllCoeditedRegions(
#'     regions_gr = input_gr,
#'     rnaEditMatrix = rnaedit_df,
#'     output = "GRanges",
#'     method = "spearman"
#'   )
#'   
#'   # identify input regions for co-edited regions
#'   AddMetaData(
#'     target_gr = coedited_gr,
#'     annot_gr = input_gr,
#'     annotType_char = "region",
#'     annotLabel_char = "inputRegion",
#'     genome = "hg19"
#'   )
#'   
AddMetaData <- function(target_gr,
                        annot_gr = NULL,
                        annotType_char = c("geneSymbol", "region"),
                        annotLabel_char = "symbol",
                        genome = c("hg38", "hg19")){
  # browser()
  
  genome <- match.arg(genome)
  annotType_char <- match.arg(annotType_char)
  
  if (annotType_char == "region" & is.null(annot_gr)) {
    stop(
      "annot_gr can't be NULL when annotType_char is 'region'. Please specify!"
    )
  }
  
  if (annotType_char == "geneSymbol" & is.null(annot_gr)) {
    
    if (genome == "hg38") {
      annot_gr <- readRDS(
        system.file(
          "extdata",
          "hg38_annoGene_gr.RDS",
          package = 'rnaEditr',
          mustWork = TRUE
        )
      )
    } else {
      annot_gr <- readRDS(
        system.file(
          "extdata",
          "hg19_annoGene_gr.RDS",
          package = 'rnaEditr',
          mustWork = TRUE
        )
      )
    }
    
  }
  
  # Find overlaps of target_gr in annot_gr
  #   in OverlapIndex, "subjectHits" (row names of annot_gr) tells us where 
  #   each row of target_gr locations in in annot_gr. We use "subjectHits" to 
  #   extract the correct metadata info for target_gr in the next steps.
  OverlapIndex <- findOverlaps(
    query = target_gr,
    subject = annot_gr,
    ignore.strand = TRUE
  )
  
  # Convert OverlapIndex into dataframe for convenience
  OverlapIndex_df <- data.frame(OverlapIndex)
  
  if (annotType_char == "region") {
    
    # Create region name based on seqnames, start and end in annot_df.
    annot_df <- data.frame(annot_gr, stringsAsFactors = FALSE)
    annot_df$region <- paste0(
      annot_df$seqnames, ":",
      annot_df$start, "-",
      annot_df$end
    )
    
    # Match each region ranges of annot_gr into target_gr based on OverlapIndex
    #   In GRanges, anything except for seqnames and ranges will go into meta
    #   data, so here we add meta data columns with $ into target_gr.
    OverlapIndex_df$region <- annot_df$region[subjectHits(OverlapIndex)]
    
  } else {
    # Find symbol names for each query in OverlapIndex_df
    OverlapIndex_df$region <- annot_gr$symbol[subjectHits(OverlapIndex)]
  }
  
  # In GRanges, anything except for seqnames and ranges will go into meta data,
  # so here we add meta data columns with $ into target_gr

  target_gr$region <- vapply(
    X = seq_len(length(target_gr)),
    FUN = function(query){
      OverlapIndexRow_df <- OverlapIndex_df[
        OverlapIndex_df$queryHits == query,
      ]
      
      if (nrow(OverlapIndexRow_df) == 0) {
        ""
      } else {
        paste(OverlapIndexRow_df$region, collapse = ";")
      }
      
    },
    FUN.VALUE = character(1)
  )
  
  # Name the newly added annotaion metadata column
  if (!is.null(annotLabel_char)) {
    
    # Find the position of the newly added annotation metadata column. As it's
    # always added the last. I just need to find the length of meta columns
    annotNameLen <- length(mcols(target_gr))
    
    # Only change the name of the last added meta column from "region" to input
    # annotLabel_char
    names(mcols(target_gr))[annotNameLen] <- annotLabel_char 
    
  }
  
  target_gr

}


