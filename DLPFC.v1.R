basedir <- "E:/jmcco/Downloads/GSE208338_RAW/"
setwd(basedir)
package_loc <- paste(basedir, "lib", sep = "")

# Bioconductor Installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = package_loc)
  BiocManager::install(version = "3.17")

BiocManager::install(c("affy", "caret", "e1071", "genefilter", "ggplot2", "locfit", "nnet", "oligo", "pd.huex.1.0.st.v2", "BiocGenerics", "withr"), force = TRUE, lib = package_loc)

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

library(locfit, lib.loc = package_loc, attach.required = TRUE)

# Get all Affymetrix CEL files
celFiles <- list.celfiles(paste(basedir, "CEL", sep = ""), full.names=TRUE)
# Get the data from all CEL files
# Returns ExonFeatureSet
rawData <- read.celfiles(celFiles)
# Extract the perfect match probe-level intensities
str(rawData)
head(rawData)
dim(pm(rawData))
# Platform used for gene information
annotation(rawData)

# I. Pre-processing
# http://www.sthda.com/english/wiki/affymetrix-cel-files
# I.1. Turn probe-level information to gene-level data by normalization using RMA
# Returns ExpressionSet
e <- rma(rawData)
# Remove rawData variable to save space
rm(rawData, celFiles)

# Show details about the preprocessed values
str(e)
head(e)
boxplot(e)


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
diagnosis.fact <- factor(diagnosis, ordered = FALSE)
withDisorder <- diagnosis
# Factor non-control to DIS which means disorder
# CTL - 9, DIS - 1
withDisorder[withDisorder != 9] <- 1
withDisorder.fact <- factor(withDisorder, ordered = FALSE)
withDisorder.fact <- relevel(withDisorder.fact, ref = "9")

# Get expressions
gene_expressions <- exprs(e)
data <- gene_expressions
data <- t(data)
data <- as.data.frame(data)
print(paste("After RMA: ", paste(c("Samples: ", "Expressions: "), dim(data), collapse = " ")))

# I.2. Dimension reduction
# I.2.1 Gene filtering through t-test
alpha <- 0.001
library(genefilter, lib.loc = package_loc, attach.required = TRUE)
# I.2.1.1. T-test to get common genes across all disorder
f1 <- function(x) (t.test(x ~ withDisorder)$p.value < alpha)
ff <- filterfun(f1)
selected <- genefilter(t(data), ff)
# Cuts it down to 169 x 47
data_genefilter.t_test <- data[, which(selected)]
print(paste("After gene filter: ", paste(c("Samples: ", "Expressions: "), dim(data_genefilter.t_test), collapse = " ")))
sig_features <- colnames(data_genefilter.t_test)
print("Significant Features: ")
print(sig_features)
print(with(data_genefilter.t_test, table(withDisorder)))
data_genefilter.t_test <- cbind(data_genefilter.t_test, withDisorder)
data_genefilter.t_test[, 'withDisorder'] <- relevel(factor(data_genefilter.t_test[, 'withDisorder'], ordered = FALSE), ref = "9")

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

# II. Processing
## Metrics for logistic regression and KNN
err_metric <- function(tp, tn, fn, fp) {
  print(paste("Accuracy: ", accuracy <- (tp + tn)/(tp + tn + fp + fn)))
  print(paste("Sensitivity: ", sensitivity <- (tp) / (tp + fn)))
  print(paste("Specificity: ", specificity <- (tn) / (tn + fp)))
  print(paste("Precision: ", precision <- (tp) / (tp + fp)))
  print(paste("F-score: ", fscore <- 2 * ((precision * sensitivity) / (precision + sensitivity))))
}

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

