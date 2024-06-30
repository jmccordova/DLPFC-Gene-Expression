# Part 1: Setting up
# Part 1.1: Sets the location of the data to be used and where the packages should be put
datadir <- "E:/jmcco/Downloads/BNF 300.2 Data/GSE208338_RAW/"
#probedir <- "E:/jmcco/Downloads/BNF 300.2 Data/Affymetrix_HuEx/"
probedir <- "E:/jmcco/Downloads/BNF 300.2 Data/HuEx-1_0-st-v2-na36-hg19 Probeset/"
setwd(datadir)
package_loc <- paste(datadir, "lib", sep = "")

# Part 1.2: Package and library installations
# Bioconductor Installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = package_loc)
library(BiocManager, lib.loc = package_loc)

BiocManager::install(
  c("tzdb", "vroom", "readr", "ggplot2", 
    "backports", "ggcorrplot", "locfit", "oligo", "ggpubr", 
    "pd.huex.1.0.st.v2", "BiocGenerics","Biobase",
    
    # Part 2: Explore
    "genefilter", "limma", "ggvenn",
    "multtest", "annaffy",
    
    # Part 3: Dimension Reduction
    "withr", "corrr", "idm", "irlba", "PCAtools", 
    "RMTstat", "biomaRt", "pROC", "nFactors",
    "EFA.dimensions", 
    "corrplot", "factoextra", "car", 
    "jsonlite",
    
    # Part 4: Analysis
    "caret", "e1071", "lattice",  "naivebayes",
    "class", "gmodels", 
    "rpart", "rpart.plot", "Hmisc", 
    "nnet", "rminer",
    "randomForest",
    "MASS", "Metrics"), 
  force = TRUE, 
  dependencies = TRUE, 
  lib = package_loc
)

library(Biobase, lib.loc = package_loc); library(BiocSingular, lib.loc = package_loc); library(tzdb, lib.loc = package_loc); library(vroom, lib.loc = package_loc); library(readr, lib.loc = package_loc)
library(S4Vectors, lib.loc = package_loc); library(IRanges, lib.loc = package_loc); library(XVector, lib.loc = package_loc); library(GenomeInfoDb, lib.loc = package_loc); library(Biostrings, lib.loc = package_loc)
library(withr, lib = package_loc); library(backports, lib.loc = package_loc); library(ggcorrplot, lib.loc = package_loc); library(ggpubr, lib = package_loc); 
library(BiocGenerics, lib.loc = package_loc); library(dplyr, lib.loc = package_loc); 
library(oligoClasses, lib.loc = package_loc); library(memoise, lib.loc = package_loc); library(pd.huex.1.0.st.v2, lib.loc = package_loc); library(oligo, lib.loc = package_loc, attach.required = TRUE)

# Part 2: Exploration
library(withr, lib = package_loc); library(ggplot2, lib = package_loc); 
library(genefilter, lib = package_loc); library(memoise, lib = package_loc); library(limma, lib = package_loc); library(labeling, lib = package_loc); library(farver, lib = package_loc); library(ggvenn, lib = package_loc)
library(multtest, lib = package_loc);  library(pkgconfig, lib = package_loc);  library(GO.db, lib = package_loc); library(annaffy, lib = package_loc)

# Part 3: Dimension Reduction
library(locfit, lib.loc = package_loc)
library(corrr, lib.loc = package_loc); library(idm, lib.loc = package_loc); library(irlba, lib.loc = package_loc) 
library(PCAtools, lib.loc = package_loc); library(RMTstat, lib.loc = package_loc); library(rappdirs, lib.loc = package_loc); library(biomaRt, lib.loc = package_loc); library(cowplot, lib.loc = package_loc); library(ggplotify, lib.loc = package_loc)
library(pROC, lib.loc = package_loc); library(withr, lib.loc = package_loc); 
library(EFA.dimensions, lib.loc = package_loc)
library(corrplot, lib.loc = package_loc); library(factoextra, lib.loc = package_loc); library(car, lib.loc = package_loc)
library(jsonlite, lib.loc = package_loc); library(backports, lib.loc = package_loc); library(Hmisc, lib.loc = package_loc); 

# Part 4: Analysis
library(e1071, lib.loc = package_loc); library(naivebayes, lib.loc = package_loc)
library(class, lib.loc = package_loc); library(gmodels, lib.loc = package_loc)
library(parallel, lib.loc = package_loc); library(doParallel, lib.loc = package_loc)
library(rpart, lib.loc = package_loc); library(rpart.plot, lib.loc = package_loc)
library(nnet, lib.loc = package_loc); library(rminer, lib.loc = package_loc)
library(randomForest, lib.loc = package_loc)
library(MASS, lib.loc = package_loc); library(Metrics, lib.loc = package_loc); 
library(caret, lib.loc = package_loc); 

# Part 1.3: Initialize cores for parallel processing
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
clusterEvalQ(cluster, .libPaths("E:/jmcco/Downloads/BNF 300.2 Data/GSE208338_RAW/lib"))
registerDoSEQ()
on.exit(stopCluster(cluster))

# Part 1.4: Define 'not in' func
'%ni%' <- Negate('%in%')

# Part 1.5: Extracting data
# Part 1.5.1: Get all Affymetrix CEL files from the directory
celFiles <- list.celfiles(paste(datadir, "CEL", sep = ""), full.names=TRUE)
# Part 1.5.2: Get the data from all CEL files
# Returns ExonFeatureSet
rawData <- read.celfiles(celFiles)

# Part 1.5.3: Get the probe-level information of ExonFeatureSet
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

# Part 1.5.4: Extracting probeset information
huex.probes <- read.csv(paste(probedir, 'HuEx-1_0-st-v2.na36.hg19.probeset.csv', sep = ''), comment.char = "#", header = TRUE)

# A which for multidimensional arrays.
# Mark van der Loo 16.09.2011
#
# A Array of booleans
# returns a sum(A) x length(dim(A)) array of multi-indices where A == TRUE
#
multi.which <- function(A){
  if ( is.vector(A) ) return(which(A))
  d <- dim(A)
  T <- which(A) - 1
  nd <- length(d)
  t( sapply(T, function(t){
    I <- integer(nd)
    I[1] <- t %% d[1]
    sapply(2:nd, function(j){
      I[j] <<- (t %/% prod(d[1:(j-1)])) %% d[j]
    })
    I
  }) + 1 )
}

show_perfect_collinearity <- function(df) {
  df <- cor(t(df))
  print(multi.which(df == 1))
  remove(df)
}