## dataLoadingFunctions.R

## Erich S. Huang
## erich@post.harvard.edu
## erich.huang@sagebase.org
## Sage Bionetworks
## Seattle, Washington

loadBreastCaRnaSeqData <- function(){
  require(synapseClient)
  cat('Loading TCGA Breast Cancer RPKM Data Object from Synapse metaGenomics Project\n')
  brcaEnt <- loadEntity('syn595321')
  brcaEset <- brcaEnt$objects$eset
  brcaRnaSeq <- exprs(brcaEset)
}

loadPik3caIndicator <- function(){
  require(synapseClient)
  cat('Loading TCGA Breast Cancer PIK3CA Mutation Indicator\n')
  pikEnt <- loadEntity('syn1589828')
  pik3caInd <- read.table(list.files(pikEnt$cacheDir, full.names = T))
}