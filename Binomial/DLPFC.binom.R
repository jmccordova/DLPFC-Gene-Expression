# Part 1: Setting up

# Part 1.1: Sets the location of the data to be used and where the packages should be put
datadir <- "E:/jmcco/Downloads/GSE208338_RAW/"
setwd(datadir)
package_loc <- paste(datadir, "lib", sep = "")

# Part 1.2: Package and library installations
# Bioconductor Installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = package_loc)

BiocManager::install(c("affy", "caret", "e1071", "genefilter", "ggplot2", "locfit", "nnet", "oligo", "pd.huex.1.0.st.v2", "BiocGenerics", "withr", "randomForest", "pROC"), force = TRUE, lib = package_loc)

library(Biobase, lib.loc = package_loc)

library(S4Vectors, lib.loc = package_loc)
library(IRanges, lib.loc = package_loc)
library(XVector, lib.loc = package_loc)
library(GenomeInfoDb, lib.loc = package_loc)
library(Biostrings, lib.loc = package_loc)

library(BiocGenerics, lib.loc = package_loc)

library(oligoClasses, lib.loc = package_loc)

library(memoise, lib.loc = package_loc)
library(pd.huex.1.0.st.v2, lib.loc = package_loc)
library(oligo, lib.loc = package_loc, attach.required = TRUE)

library(caret, lib.loc = package_loc)
library(locfit, lib.loc = package_loc)

library(randomForest, lib.loc = package_loc)
library(pROC, lib.loc = package_loc)

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

# Part 3: Analysis
# Part 3.1: Pre-processing
# Part 3.1.1: Turn probe-level information to gene-level data by normalization using RMA
# RMA - Robust Multi-Array Average (http://www.sthda.com/english/wiki/affymetrix-cel-files)
# Returns ExpressionSet
e <- rma(rawData)
# Remove rawData variable to save space
rm(rawData, celFiles)

# Show details about the preprocessed values
str(e)
head(e)
boxplot(e)

# Part 3.1.2: Create factors for each diagnosis
# Part 3.1.2.1: Replace string to integer classification
# Diagnosis from each sample in GEO
# CTL - 9, BPD - 1, MDD - 2, SCZ - 3
diagnosis <- c(
  # 299-309
  9, 3, 2, 3, 1, 1, 2, 3, 3, 9, 3,
  # 310-320
  2, 1, 3, 9, 3, 2, 2, 1, 1, 9, 1,
  # 321-331
  1, 9, 1, 9, 1, 9, 9, 9, 3, 2, 9,
  # 332-342
  2, 9, 3, 2, 9, 2, 3, 9, 3, 2, 2,
  # 343-353
  9, 2, 9, 3, 9, 3, 3, 3, 3, 9, 3,
  # 354-364
  9, 9, 3, 9, 3, 3, 3, 9, 3, 9, 3,
  # 365-375
  3, 9, 9, 3, 9, 9, 9, 3, 3, 9, 3,
  # 376-386
  9, 3, 3, 9, 3, 3, 3, 3, 3, 9, 3,
  # 387-397
  9, 3, 3, 9, 3, 9, 3, 9, 3, 3, 3,
  # 398-408
  9, 9, 3, 3, 1, 9, 9, 9, 9, 2, 1,
  # 409-419
  9, 9, 2, 1, 9, 2, 1, 1, 3, 2, 1,
  # 420-430
  9, 3, 9, 3, 2, 3, 3, 2, 3, 3, 2,
  # 431-441
  3, 2, 3, 2, 3, 2, 3, 2, 9, 3, 2,
  # 442-452
  9, 9, 3, 9, 9, 3, 3, 9, 3, 9, 3,
  # 453-463
  3, 3, 9, 3, 9, 3, 9, 9, 9, 9, 9,
  # 464-476
  9, 9, 3, 3
)
# Part 3.1.2.2: Create factors for each diagnosis
diagnosis.fact <- factor(diagnosis, ordered = FALSE)

# Part 3.1.2.3: Create factors for with disorder and control
withDisorder <- diagnosis
# Factor non-control to DIS which means disorder
# CTL - 9, DIS - 1
withDisorder[withDisorder != 9] <- 1
withDisorder.fact <- factor(withDisorder, ordered = FALSE)
withDisorder.fact <- relevel(withDisorder.fact, ref = "9")

# Part 3.1.3: Get gene expressions
gene_expressions <- exprs(e)
data <- gene_expressions
data <- t(data)
data <- as.data.frame(data)
print(paste("After RMA: ", paste(c("Samples: ", "Expressions: "), dim(data), collapse = " ")))
# Append diagnosis.fact factors to the dataset
data.multinomial <- data
data.multinomial[, 'diagnosis'] <- diagnosis.fact
# Append withDisorder.fact factors to the dataset
data.binomial <- data
data.binomial[, 'diagnosis'] <- withDisorder.fact

# Part 3.1.4: Alpha computation
# Perform Bonferroni correction
get_alpha <- function(alpha_orig, n) {
  return(alpha_orig/n)
}
alpha <- 0.05
alpha <- get_alpha(alpha, ncol(data))

# Part 3.1.5: Splitting dataset
## For binomial
'%ni%' <- Negate('%in%')  # define 'not in' func
options(scipen=999)  # prevents printing scientific notations.
set.seed(100)
index.binomial <- createDataPartition(data.binomial$diagnosis, p=0.75, list = F)
trainset.binomial <- data.binomial[index.binomial, ]
testset.binomial <- data.binomial[-index.binomial, ]
## For multinomial
'%ni%' <- Negate('%in%')  # define 'not in' func
options(scipen=999)  # prevents printing scientific notations.
set.seed(100)
index.multinomial <- createDataPartition(data.multinomial$diagnosis, p=0.75, list = F)
trainset.multinomial <- data.multinomial[index.multinomial, ]
testset.multinomial <- data.multinomial[-index.multinomial, ]

# Part 3.1.6: Dimension reduction using Random Forest
set.seed(100)
# Random Forest Function
#
# @param Dataframe train - training dataset
# @param Dataframe test - testing dataset
# @param Integer n - number of trees for random forest (ntree)
# @param Integer m - number of tries for random forest (mtry)
# @param Boolean getImportantVars
# @return model
randomForestFunction <- function(train, test, n, m, getImportantVars) {
  model <- randomForest(
    x = train[, colnames(train) != "diagnosis"],
    y = train$diagnosis, 
    ntree = n, 
    mtry = m,
    importance = getImportantVars
  )
  
  pred <- pred <- predict(model, newdata = test, type = "response")
  confMatrix <- table(test$diagnosis, pred)
  
  importantVars <- c()
  if(getImportantVars) {
    varImpPlot(model)
    importantVars <- importance(model)
  }
  
  return(list(
    model = model,
    confMatrix = confMatrix,
    importantVars = importantVars
  ))
}

# Metrics for model efficiency
# @param Integer[][] confusionMatrix
# @return Float[]
getMetrics <- function(confusionMatrix, printMetrics = TRUE, asPercent = TRUE) {
  tp <- confusionMatrix[1,1]
  tn <- confusionMatrix[2,2]
  fn <- confusionMatrix[1,2]
  fp <- confusionMatrix[2,1]
  
  accuracy <- (tp + tn)/(tp + tn + fp + fn)
  sensitivity <- (tp) / (tp + fn)
  specificity <- (tn) / (tn + fp)
  precision <- (tp) / (tp + fp)
  fscore <- 2 * ((precision * sensitivity) / (precision + sensitivity))
  
  if (asPercent) {
    accuracy <- paste(accuracy * 100, '%', sep = '')
    sensitivity <- paste(sensitivity * 100, '%', sep = '')
    specificity <- paste(specificity * 100, '%', sep = '')
    precision <- paste(precision * 100, '%', sep = '')
    fscore <- paste(fscore * 100, '%', sep = '')
  }
  
  if (printMetrics) {
    print(paste("Accuracy: ", accuracy))
    print(paste("Sensitivity: ", sensitivity))
    print(paste("Specificity: ", specificity))
    print(paste("Precision: ", precision))
    print(paste("F-score: ", fscore))
  }
  
  return(c(accuracy = accuracy, sensitivity = sensitivity, specificity = specificity, precision = precision, fscore = fscore))
}

# Parameter tuning
# Get 5 random mtries from 1 to 22011
print(paste("mtry: ", mtries <- sort.int(sample(ncol(data)-1, 5))))
# ntree default is 500 but we will do a trial and error
# The range will go less than 500 (http://lisa.ulb.ac.be/image/publifiles/100/Mcs2001_final.pdf) and greater than 500
# Feb 18 2024 update: It looks like ntree = 100 has the highest accuracy. Let's try lowering it.
# Feb 25 2024 update: ntree = 25 has the highest accuracy. ntree = 22 seems promising, too, with low OOB.
for(ntree in c(20, 50, 100, 200, 500, 1000)) {
  for(mtry in mtries) {
    print(paste("ntree = ", ntree, " mtry = ", mtry))
    rf.model <- randomForestFunction(train = trainset.binomial, test = testset.binomial, n = ntree, m = mtry, getImportantVars = FALSE)
    rf.metrics <- getMetrics(rf.model$confMatrix, TRUE, TRUE)
    print('-----------------------')
  }
}

# Choose the best parameters for the model
# Tuned
ntree <- 50; mtry <- 6860    # OOB Error = 36.72% Accuracy = 63.41%
#tuneRF(x = trainset.binomial[, colnames(trainset.binomial) != "diagnosis"],
#       y = trainset.binomial$diagnosis, mtryStart = mtry, ntreeTry = ntree, stepFactor=2, improve=0.05,
#       trace=TRUE, plot=TRUE)

rf.model <- randomForestFunction(trainset.binomial, testset.binomial, ntree, mtry, TRUE)
rf.metrics <- getMetrics(rf.model$confMatrix, TRUE, TRUE)
rf.model$model$importance
# Importance by Accuracy (get only MeanDecreaseAccuracy > 0)
features.acc <- rownames(rf.model$model$importance[which(rf.model$model$importance[, 3] > 0), ])
# Importance by GINI (get only MeanDecreaseGINI > 0)
features.gini <- rownames(rf.model$model$importance[which(rf.model$model$importance[, 4] > 0), ])
# Combine all features into one
features <- unique(c(features.acc, features.gini))


# I.2.2. Data split
library(withr, lib.loc = package_loc)
library(ggplot2, lib.loc = package_loc)
library(caret, lib.loc = package_loc)
'%ni%' <- Negate('%in%')  # define 'not in' func
# options(scipen=999)  # prevents printing scientific notations.
# set.seed(100)
# index <- createDataPartition(data_genefilter$withDisorder, p=0.75, list = F)
# trainset <- data_genefilter[index, ]
# testset <- data_genefilter[-index, ]


# II.1. Logistic Regression for t-test gene filter
library(e1071, lib.loc = package_loc, attach.required = TRUE)
modellogit <- glm(withDisorder ~ ., family = binomial(link = "logit"), data = data_genefilter.t_test)
summary(modellogit)
print(cbind("coeff" = modellogit$coefficients, "odds ratio" = (exp(modellogit$coefficients) - 1) * 100)) # Odds ratio
# Predict
pred.modellogit <- predict(modellogit, newdata = data_genefilter.t_test[, colnames(data_genefilter.t_test) %ni% "withDisorder"], type = "response") > 0.5
confMatrix.modellogit <- table(data_genefilter.t_test$withDisorder, pred.modellogit)
print(confMatrix.modellogit)
err_metric(confMatrix.modellogit[1,1], confMatrix.modellogit[2,2], confMatrix.modellogit[1,2], confMatrix.modellogit[2,1])

# I.2.2 Gene filtering through ANOVA
data_genefilter.anova <- c()
data_genefilter.anova.colnames <- c()
i <- 1
for(feature in data) {
  p_value <- summary(aov(feature ~ diagnosis.fact))[[1]][,5][1]
  if (p_value < alpha) {
    data_genefilter.anova <- cbind(data_genefilter.anova, feature)
    data_genefilter.anova.colnames <- cbind(data_genefilter.anova.colnames, colnames(data[i]))
    # print(paste(i, colnames(data[i]), p_value))
  }
  
  i <- i+1
}

data_genefilter.anova <- data.frame(data_genefilter.anova, row.names = row.names(data))
colnames(data_genefilter.anova) <- data_genefilter.anova.colnames
remove(data_genefilter.anova.colnames)
print(paste("After gene filter: ", paste(c("Expressions: ", "Samples: "), dim(data_genefilter.anova), collapse = " ")))
print("Significant Features: ")
print(colnames(data_genefilter.anova))
print(with(data_genefilter.anova, table(diagnosis)))
# Uncomment when using glm()
data_genefilter.anova <- cbind(data_genefilter.anova, diagnosis)
data_genefilter.anova[, 'diagnosis'] <- as.factor(data_genefilter.anova[, 'diagnosis'])

# II.2. Logistic Regression for ANOVA gene filter
data_genefilter.anova$diagnosis <- relevel(data_genefilter.anova$diagnosis, ref = "9")
library(nnet, lib.loc = package_loc, attach.required = TRUE)
modellogit <- multinom(diagnosis ~ ., data = data_genefilter.anova)
# Does not have convergence
# library(e1071, lib.loc = package_loc, attach.required = TRUE)
# modellogit <- glm(diagnosis ~ ., family = poisson(link = "log"), data = data_genefilter.anova)
summary(modellogit)
modellogit.z <- summary(modellogit)$coefficients/summary(modellogit)$standard.errors
modellogit.p <- (1 - pnorm(abs(modellogit.z), 0, 1)) * 2
print(modellogit.p)
print(cbind("coeff" = coef(modellogit), "odds ratio" = exp(coef(modellogit)))) # Odds ratio
# Predict
pred.modellogit <- predict(modellogit, newdata = data_genefilter.anova[, colnames(data_genefilter.anova) %ni% "diagnosis"], type = "probs") > 0.5
# Create similar table
data_genefilter.anova.diagnosis.logi <- pred.modellogit
i <- 1
for (sample in rownames(data_genefilter.anova)) {
  if(data_genefilter.anova.diagnosis.logi[sample, toString(diagnosis[i])] != TRUE) {
    data_genefilter.anova.diagnosis.logi[sample, toString(diagnosis[i])] <- TRUE
    data_genefilter.anova.diagnosis.logi[sample, colnames(data_genefilter.anova.diagnosis.logi) %ni% toString(diagnosis[i])] <- FALSE
  }
  i <- i + 1
}
confMatrix.modellogit <- table(data_genefilter.anova.diagnosis.logi, pred.modellogit)
print(confMatrix.modellogit)
err_metric(confMatrix.modellogit[1,1], confMatrix.modellogit[2,2], confMatrix.modellogit[1,2], confMatrix.modellogit[2,1])

