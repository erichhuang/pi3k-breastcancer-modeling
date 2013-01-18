## pilotPIK3CA-Model.R

## Erich S. Huang
## erich@post.harvard.edu
## erich.huang@sagebase.org
## Sage Bionetworks
## Seattle, Washington

## REQUIRED LIBRARIES
require(synapseClient)
require(rGithubClient)
require(Biobase)
require(randomForest)

## SOURCE IN THE GITHUB CODE REPOSITORY FOR THIS PROJECT
projectRepo <- getRepo(repository = '/erichhuang/pi3k-breastcancer-modeling')
sourceRepoFile(projectRepo, 'Functions/dataLoadingFunctions.R')

## LOAD DATA
cat('Loading the TCGA breast cancer RNA Seq and PIK3CA binary mutation data\n')
# tcga breast cancer data: rpkm normalized RNA seq data
breastCaRpkm <- loadBreastCaRnaSeqData()

# simple binary indicator of pik3ca mutations as provided by cbio
pik3caInd <- loadPik3caIndicator()

## INTERSECT THE DATA
cat('Intersecting the unique sample IDs\n')
# first make rownames for the pik3ca indicator
rownames(pik3caInd) <- pik3caInd[ , 1]

# remove extraneous characters from the RNA seq colnames
colnames(breastCaRpkm) <- sapply(strsplit(colnames(breastCaRpkm), '-'), function(x){
  paste(x[1:3], collapse="-")
})

## actual intersect
idIntersect <- intersect(rownames(pik3caInd), colnames(breastCaRpkm))
pik3caSubset <- pik3caInd[idIntersect, ]
breastCaRnaSubset <- breastCaRpkm[ , idIntersect]

## INSPECT FOR OUTLIERS & REMOVE THEM
cat('Identifying and removing two outlier samples\n')
svdDat <- svd(breastCaRnaSubset)
plot(svdDat$v[ , 1], svdDat$v[ , 2])

## Looks like removing samples > 0.3 on Eigengene 2 is an easy win
outlierSamps <- grep('TRUE', svdDat$v[ , 2] > 0.3)

breastCaRnaSubset <- breastCaRnaSubset[ , -outlierSamps]
pik3caSubset <- pik3caSubset[-outlierSamps, ]

## RANDOMLY DIVIDE THE DATA INTO TRAINING AND VALIDATION COHORTS
cat('Randomly dividing the dataset into firewalled training and validation cohorts\n')
set.seed(130118)
cohortInd <- sample(1:dim(breastCaRnaSubset)[2], dim(breastCaRnaSubset)[2]/2)

trainClass <- pik3caSubset[cohortInd, 2]
trainExpress <- breastCaRnaSubset[ , cohortInd]

validClass <- pik3caSubset[-cohortInd, 2]
validExpress <- breastCaRnaSubset[ , -cohortInd]

## GENERATE THE NEW RANDOM FOREST MODEL PHASE 1
cat('Building a preliminary random forest model on PIK3CA mutant versus PIK3CA wildtype samples\n')
pi3kPrelimModel <- randomForest(t(trainExpress),
                                as.factor(trainClass),
                                ntree = 500,
                                do.trace = 2,
                                importance = TRUE,
                                proximity = TRUE)

## Obtain train hat predictions
# trainScoreHat <- predict(pi3kPrelimModel, t(trainExpress), type = 'prob')

## VALIDATE THE MODEL
validScoreHat <- predict(pi3kPrelimModel, t(validExpress), type = 'prob')

# # Nonparametric ranksum
# validNumeric <- as.numeric(validClass)
# rankSum <- wilcox.test(validScoreHat[validNumeric == 0, 2], 
#                        validScoreHat[validNumeric == 1, 2])
#   
