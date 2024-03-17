# Part 1: Setting up

# Part 1.1: Sets the location of the data to be used and where the packages should be put
datadir <- "E:/jmcco/Downloads/BNF 300.2 Data/GSE208338_RAW/"
probedir <- "E:/jmcco/Downloads/BNF 300.2 Data/Affymetrix_HuEx/"
setwd(datadir)
package_loc <- paste(datadir, "lib", sep = "")

# Part 1.2: Package and library installations
# Bioconductor Installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = package_loc)

BiocManager::install(
  c("tzdb", "vroom", "readr", "affy", "ggplot2", 
    "backports", "Hmisc", "ggcorrplot", "locfit", "oligo", 
    "pd.huex.1.0.st.v2", "BiocGenerics",
    
    # Part 3: Dimension Reduction
    "withr", "corrr", "idm", "irlba", "PCAtools", 
    "RMTstat", "biomaRt", "pROC",
    
    # Part 4: Analysis
    "caret", "e1071", "lattice",  "naivebayes",
    "class", "gmodels", 
    "rpart", "rpart.plot",
    "nnet",
    "randomForest"), 
  #force = TRUE, 
  dependencies = TRUE, 
  lib = package_loc
)

library(Biobase, lib.loc = package_loc); library(tzdb, lib.loc = package_loc); library(vroom, lib.loc = package_loc); library(readr, lib.loc = package_loc)
library(S4Vectors, lib.loc = package_loc); library(IRanges, lib.loc = package_loc); library(XVector, lib.loc = package_loc); library(GenomeInfoDb, lib.loc = package_loc); library(Biostrings, lib.loc = package_loc)
library(backports, lib.loc = package_loc); library(Hmisc, lib.loc = package_loc); library(ggcorrplot, lib.loc = package_loc)
library(BiocGenerics, lib.loc = package_loc)
library(oligoClasses, lib.loc = package_loc); library(memoise, lib.loc = package_loc); library(pd.huex.1.0.st.v2, lib.loc = package_loc); library(oligo, lib.loc = package_loc, attach.required = TRUE)

# Part 2: Exploration
library(psych, lib = package_loc)

# Part 3: Dimension Reduction
library(locfit, lib.loc = package_loc)
library(corrr, lib.loc = package_loc); library(idm, lib.loc = package_loc); library(irlba, lib.loc = package_loc) 
library(PCAtools, lib.loc = package_loc); library(RMTstat, lib.loc = package_loc); library(biomaRt, lib.loc = package_loc)
library(pROC, lib.loc = package_loc); library(withr, lib.loc = package_loc); 

# Part 4: Analysis
library(caret, lib.loc = package_loc); library(dplyr, lib.loc = package_loc); library(e1071, lib.loc = package_loc); library(naivebayes, lib.loc = package_loc)
library(class, lib.loc = package_loc); library(gmodels, lib.loc = package_loc)
library(parallel, lib.loc = package_loc); library(doParallel, lib.loc = package_loc)
library(rpart, lib.loc = package_loc); library(rpart.plot, lib.loc = package_loc)
library(nnet, lib.loc = package_loc)
library(randomForest, lib.loc = package_loc)

'%ni%' <- Negate('%in%')  # define 'not in' func

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
clusterEvalQ(cluster, .libPaths("E:/jmcco/Downloads/BNF 300.2 Data/GSE208338_RAW/lib"))

# Part 2: Preparing data
# Part 2.1: Get all Affymetrix CEL files from the directory
celFiles <- list.celfiles(paste(datadir, "CEL", sep = ""), full.names=TRUE)
# Get the data from all CEL files
# Returns ExonFeatureSet
rawData <- read.celfiles(celFiles)
# Part 2.2: Get information about the dataset
str(rawData)
head(rawData)
dim(pm(rawData))
# Platform used for gene information
annotation(rawData)

# Get the probe-level information of ExonFeatureSet
# mm - mismatch matrix
# pm - perfect match matrix
get_matrix <- function(type, data) {
  switch (type,
          'intensity' = intensity(data),
          'mm' = mm(data, subset = NULL, target='core'),
          'pm' = pm(data, subset = NULL, target='core'),
          'bg' = bg(data, subset = NULL),
          'all' = c(intensity(data), mm(data, subset = NULL, target='core'), pm(data, subset = NULL, target='core'), bg(data, subset = NULL))
  )
}

# Part 2.3: Get probe IDs
ensemblIDs <- read_tsv(paste(probedir, 'Affymetrix_HuEx_microarray_probeset_IDs_to_Ensemble_IDs.tsv', sep = ''))
# Part 2.4: Get unique gene IDs
transcriptIDs <- ensemblIDs[c('id_internal_huex', 'transcript_id')]
transcriptIDs <- transcriptIDs[!duplicated(transcriptIDs$id_internal_huex), ]
# Retain only the genes which probeset ID are found in the DLPFC dataset
transcriptIDs <- transcriptIDs[transcriptIDs$id_internal_huex %in% colnames(data), ]
# Combine the probeIDs from DLPFC dataset that have no gene IDs
transcriptIDs.missing <- colnames(data)[colnames(data) %ni% transcriptIDs$id_internal_huex]
transcriptIDs.missing <- data.frame(matrix(c(transcriptIDs.missing, transcriptIDs.missing), ncol = 2))
names(transcriptIDs.missing) <- c('id_internal_huex', 'transcript_id')
transcriptIDs <- rbind(transcriptIDs, transcriptIDs.missing)
transcriptIDs <- arrange(transcriptIDs, id_internal_huex)
remove(transcriptIDs.missing)

# Part 2.5: Create a metadata array to be used for PCA
data.pca.metadata <- matrix(c(transcriptIDs$id_internal_huex, transcriptIDs$transcript_id), nrow = 2, ncol = nrow(transcriptIDs), byrow = TRUE)
data.pca.metadata <- data.frame(data.pca.metadata)
names(data.pca.metadata) <- data.pca.metadata[1, ]
data.pca.metadata <- data.pca.metadata[-1, ]
