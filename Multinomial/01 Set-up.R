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
  c("tzdb", "vroom", "readr", "ggplot2", 
    "backports", "Hmisc", "ggcorrplot", "locfit", "oligo", 
    "pd.huex.1.0.st.v2", "BiocGenerics",
    
    # Part 2: Explore
    "ggpubr", "genefilter", "limma",
    "multtest", "annaffy", "hgu95av2.db",
    
    # Part 3: Dimension Reduction
    "withr", "corrr", "idm", "irlba", "PCAtools", 
    "RMTstat", "biomaRt", "pROC", "nFactors",
    "EFA.dimensions",
    
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
library(ggpubr, lib = package_loc); library(genefilter, lib = package_loc); library(limma, lib = package_loc)
library(multtest, lib = package_loc); library(annaffy, lib = package_loc); library(hgu95av2.db, lib = package_loc)

# Part 3: Dimension Reduction
library(locfit, lib.loc = package_loc)
library(corrr, lib.loc = package_loc); library(idm, lib.loc = package_loc); library(irlba, lib.loc = package_loc) 
library(PCAtools, lib.loc = package_loc); library(RMTstat, lib.loc = package_loc); library(biomaRt, lib.loc = package_loc)
library(pROC, lib.loc = package_loc); library(withr, lib.loc = package_loc); 
library(EFA.dimensions, lib.loc = package_loc)

# Part 4: Analysis
library(caret, lib.loc = package_loc); library(dplyr, lib.loc = package_loc); library(e1071, lib.loc = package_loc); library(naivebayes, lib.loc = package_loc)
library(class, lib.loc = package_loc); library(gmodels, lib.loc = package_loc)
library(parallel, lib.loc = package_loc); library(doParallel, lib.loc = package_loc)
library(rpart, lib.loc = package_loc); library(rpart.plot, lib.loc = package_loc)
library(nnet, lib.loc = package_loc)
library(randomForest, lib.loc = package_loc)

# Initialize cores for parallel processing
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
clusterEvalQ(cluster, .libPaths("E:/jmcco/Downloads/BNF 300.2 Data/GSE208338_RAW/lib"))
registerDoSEQ()
on.exit(stopCluster(cluster))

'%ni%' <- Negate('%in%')  # define 'not in' func

# Part 2: Preparing data
# Part 2.1: Get all Affymetrix CEL files from the directory
celFiles <- list.celfiles(paste(datadir, "CEL", sep = ""), full.names=TRUE)
# Get the data from all CEL files
# Returns ExonFeatureSet
rawData <- read.celfiles(celFiles)

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
ids.ensembl <- read_tsv(paste(probedir, 'Affymetrix_HuEx_microarray_probeset_IDs_to_Ensemble_IDs.tsv', sep = ''))
# Part 2.4: Get unique gene IDs
ids.transcript <- ids.ensembl[c('id_internal_huex', 'transcript_id')]
ids.transcript <- ids.transcript[!duplicated(ids.transcript$id_internal_huex), ]
# Retain only the genes which probeset ID are found in the DLPFC dataset
ids.transcript <- ids.transcript[ids.transcript$id_internal_huex %in% unique(probeNames(rawData)), ]
# Combine the probeIDs from DLPFC dataset that have no gene IDs
ids.transcript.missing <- unique(probeNames(rawData))[unique(probeNames(rawData)) %ni% ids.transcript$id_internal_huex]
ids.transcript.missing <- data.frame(matrix(c(ids.transcript.missing, ids.transcript.missing), ncol = 2))
names(ids.transcript.missing) <- c('id_internal_huex', 'transcript_id')
ids.transcript <- rbind(ids.transcript, ids.transcript.missing)
ids.transcript <- arrange(ids.transcript, id_internal_huex)
remove(ids.transcript.missing)