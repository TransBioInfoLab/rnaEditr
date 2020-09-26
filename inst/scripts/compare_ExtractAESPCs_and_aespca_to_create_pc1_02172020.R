library(rnaEditr)
data(rnaedit_df)

a1 <- Sys.time()

dat <- pathwayPCA::TransposeAssay(rnaedit_df, omeNames = "row")

### Create a pathwayCollection list
clust_ls <- list(region = colnames(dat))
terms <- "cluster1"
aclustpath <- pathwayPCA::CreatePathwayCollection(
  clust_ls,
  terms,
  setType = "pathways"
)

### Create OmicsPath object
dat2 <- tibble::rownames_to_column(dat, "Sample")
methOmics <- pathwayPCA::CreateOmics(
  dat2,
  aclustpath,
  response = NULL,
  respType = "none",
  centerScale = c(FALSE, FALSE)
)

### Extract PC1 from each pathway
methPCs <- pathwayPCA::ExtractAESPCs(
  methOmics, numPCs = 1,
  # parallel = TRUE,
  # numCores = 2,
  standardPCA = TRUE
)

### Organize final results
ff <- unlist(methPCs$PCs,recursive = FALSE)
dat3 <- as.data.frame(t(do.call(rbind,ff)))
rownames(dat3) <- rownames(dat)

final <- tibble::deframe(dat3)
names(final) <- rownames(dat3)

final

Sys.time() - a1 # 0.707994 secs



######  pathwayPCA internal version  ######
a2 <- Sys.time()

methPCs2 <- pathwayPCA:::aespca(
  X = t(rnaedit_df)
)
final2 <- methPCs2$oldScore$V1
# standardPCA <- TRUE
# if(standardPCA) {
#   final2 <- methPCs2$oldScore$V1
# } else {
#   final2 <- methPCs2$aesScore$V1
# }

names(final2) <- colnames(rnaedit_df)
final2

Sys.time() - a2 # 0.1442709 secs

identical(final, final2)
