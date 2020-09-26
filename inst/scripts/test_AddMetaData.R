# Example script to test function AddMetaData
# Lanyu Zhang
# created on 2020-07-27
# updated on 2020-08-12

##### 1. when annotType_char is set as "region" ################################
data(rnaedit_df)

input_gr <- TransformToGR(
  genes_char = c("CCR5", "PHACTR4"),
  type = "symbol"
)

coedited_gr <- AllCoeditedRegions(
  regions_gr = input_gr,
  rnaEditMatrix = rnaedit_df,
  output = "GRanges",
  method = "spearman"
)

# test if function will successfully give error when annotType_char is set
#   set as "region" and annot_gr is NULL.
AddMetaData(
  target_gr = coedited_gr,
  annot_gr = NULL,
  annotType_char = "region",
  annotLabel_char = "inputRegion"
)

# test if function will leave meta data column name as region when
#   annotType_char is set as "region" and annotLabel_char is NULL.
AddMetaData(
  target_gr = coedited_gr,
  annot_gr = input_gr,
  annotType_char = "region",
  annotLabel_char = NULL
)

# test if function works well when all arguments are set correctly.
AddMetaData(
  target_gr = coedited_gr,
  annot_gr = input_gr,
  annotType_char = "region",
  annotLabel_char = "inputRegion"
)

# save results which will be used in section 3.
coedited_annot_gr <- AddMetaData(
  target_gr = coedited_gr,
  annot_gr = input_gr,
  annotType_char = "region",
  annotLabel_char = "inputRegion"
)




##### 2. when annotType_char is set as "geneSymbol" ###########################
regions1_gr <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(761586:767902)
)
regions2_gr <- TransformToGR(
  genes_char = "PHACTR4",
  type = "symbol"
)
regions_gr <- c(regions1_gr, regions2_gr)

# test if function works well when all arguments are set correctly.
AddMetaData(
  target_gr = regions_gr,
  annot_gr = NULL,
  annotType_char = "geneSymbol",
  annotLabel_char = "symbol"
)

# test if function will change meta data column name successfully when
#   annotLabel_char is not set as default "symbol" but something else.
AddMetaData(
  target_gr = regions_gr,
  annot_gr = NULL,
  annotType_char = "geneSymbol",
  annotLabel_char = "geneSymbol"
)




##### 3. add more than one meta data column ####################################
AddMetaData(
  target_gr = coedited_annot_gr,
  annot_gr = NULL,
  annotType_char = "geneSymbol",
  annotLabel_char = "geneSymbol"
)
