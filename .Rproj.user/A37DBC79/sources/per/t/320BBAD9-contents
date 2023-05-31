basedir <- "E:/jmcco/Downloads/GSE208338_RAW/"
setwd(basedir)

# Bioconductor Installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = paste(basedir, "lib", sep = ""))

BiocManager::install(c("affy", "caret", "e1071", "genefilter", "ggplot2", "locfit", "oligo", "pd.huex.1.0.st.v2", "BiocGenerics", "wavelets", "withr"), force = TRUE, lib = paste(basedir, "lib", sep = ""))

library(Biobase, lib.loc = paste(basedir, "lib", sep = ""))

library(S4Vectors, lib.loc = paste(basedir, "lib", sep = ""))
library(IRanges, lib.loc = paste(basedir, "lib", sep = ""))
library(XVector, lib.loc = paste(basedir, "lib", sep = ""))
library(GenomeInfoDb, lib.loc = paste(basedir, "lib", sep = ""))
library(Biostrings, lib.loc = paste(basedir, "lib", sep = ""))

library(BiocGenerics, lib.loc = paste(basedir, "lib", sep = ""))

library(oligoClasses, lib.loc = paste(basedir, "lib", sep = ""))

library(memoise, lib.loc = paste(basedir, "lib", sep = ""))
library(pd.huex.1.0.st.v2, lib.loc = paste(basedir, "lib", sep = ""))
library(oligo, lib.loc = paste(basedir, "lib", sep = ""), attach.required = TRUE)

library(locfit, lib.loc = paste(basedir, "lib", sep = ""), attach.required = TRUE)
library(wavelets, lib.loc = paste(basedir, "lib", sep = ""), attach.required = TRUE)

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

# Diagnosis from each sample in GEO
diagnosis <- c(
  # 299-309
  "CTL", "SCZ", "MDD", "SCZ", "BPD", "BPD", "MDD", "SCZ", "SCZ", "CTL", "SCZ",
  # 310-320
  "MDD", "BPD", "SCZ", "CTL", "SCZ", "MDD", "MDD", "BPD", "BPD", "CTL", "BPD",
  # 321-331
  "BPD", "CTL", "BPD", "CTL", "BPD", "CTL", "CTL", "CTL", "SCZ", "MDD", "CTL",
  # 332-342
  "MDD", "CTL", "SCZ", "MDD", "CTL", "MDD", "SCZ", "CTL", "SCZ", "MDD", "MDD",
  # 343-353
  "CTL", "MDD", "CTL", "SCZ", "CTL", "SCZ", "SCZ", "SCZ", "SCZ", "CTL", "SCZ",
  # 354-364
  "CTL", "CTL", "SCZ", "CTL", "SCZ", "SCZ", "SCZ", "CTL", "SCZ", "CTL", "SCZ",
  # 365-375
  "SCZ", "CTL", "CTL", "SCZ", "CTL", "CTL", "CTL", "SCZ", "SCZ", "CTL", "SCZ",
  # 376-386
  "CTL", "SCZ", "SCZ", "CTL", "SCZ", "SCZ", "SCZ", "SCZ", "SCZ", "CTL", "SCZ",
  # 387-397
  "CTL", "SCZ", "SCZ", "CTL", "SCZ", "CTL", "SCZ", "CTL", "SCZ", "SCZ", "SCZ",
  # 398-408
  "CTL", "CTL", "SCZ", "SCZ", "BPD", "CTL", "CTL", "CTL", "CTL", "MDD", "BPD",
  # 409-419
  "CTL", "CTL", "MDD", "BPD", "CTL", "MDD", "BPD", "BPD", "SCZ", "MDD", "BPD",
  # 420-430
  "CTL", "SCZ", "CTL", "SCZ", "MDD", "SCZ", "SCZ", "MDD", "SCZ", "SCZ", "MDD",
  # 431-441
  "SCZ", "MDD", "SCZ", "MDD", "SCZ", "MDD", "SCZ", "MDD", "CTL", "SCZ", "MDD",
  # 442-452
  "CTL", "CTL", "SCZ", "CTL", "CTL", "SCZ", "SCZ", "CTL", "SCZ", "CTL", "SCZ",
  # 453-463
  "SCZ", "SCZ", "CTL", "SCZ", "CTL", "SCZ", "CTL", "CTL", "CTL", "CTL", "CTL",
  # 464-476
  "CTL", "CTL", "SCZ", "SCZ"
)

# Pre-processing
# http://www.sthda.com/english/wiki/affymetrix-cel-files
# Turn probe-level information to gene-level data by normalization using RMA
# Returns ExpressionSet
e <- rma(rawData)
# Remove rawData variable to save space
rm(rawData, celFiles)

# Show details about the preprocessed values
str(e)
head(e)
boxplot(e)


# Get expressions
# 22011 x 169
gene_expressions <- exprs(e)
data <- as.data.frame(gene_expressions)
print(paste("After RMA: ", paste(c("Expressions: ", "Samples: "), dim(data), collapse = " ")))
# Create a matrix of values with the head as the diagnosis
#df <- data
#colnames(df) <- diagnosis   # Replace first row (header)
# Bind diagnosis at the end
#df <- rbind(data, diagnosis)  
#rownames(df)[22012] <- "diagnosis"
# End bind

# Dimension reduction
# Gene filtering
alpha <- 0.001
library(genefilter, lib.loc = paste(basedir, "lib", sep = ""), attach.required = TRUE)
withDisorder <- diagnosis
# Factor non-control to DIS which means disorder
withDisorder[withDisorder != "CTL"] <- "DIS"
withDisorder.fac <- as.factor(withDisorder)
# Normal distribution
f1 <- function(x) (shapiro.test(x)$p.value > alpha)
# Variability
f2 <- function(x) (sd(x)/abs(mean(x)) < 0.1)
# T-Test
f3 <- function(x) (t.test(x ~ withDisorder)$p.value < alpha)
ff <- filterfun(f1,f2,f3)
selected <- genefilter(data, ff)
# Cuts it down to 55 x 169
data_genefilter <- data[which(selected), ]
print(paste("After gene filter: ", paste(c("Expressions: ", "Samples: "), dim(data_genefilter), collapse = " ")))
print(rownames(data_genefilter))

# ANOVA
data_genefilter.vector <- unname(unlist(data_genefilter))
patient.fact <- gl(ncol(data_genefilter), 1, ncol(data_genefilter)*nrow(data_genefilter))
feature.fact <- as.factor(rep(rownames(data_genefilter), each = ncol(data_genefilter)))
diagnosis.fact <- as.factor(rep(diagnosis, nrow(data_genefilter)))
model <- lm(data_genefilter.vector ~ diagnosis.fact)
residuals <- residuals(model)

# Check for independence
gene_expressions <- 1:length(data_genefilter.vector)
plot(gene_expressions, residuals)

# Check for normality
shapiro.test(residuals)

# Check for constancy of variance
library(lmtest)
bptest(model, studentize = FALSE)

# ANOVA
res.aov <- aov(data_genefilter.vector ~ diagnosis.fact)
summary(res.aov)

TukeyHSD(res.aov)

### Create training and test dataset
library(withr, lib.loc = paste(basedir, "lib", sep = ""))
library(ggplot2, lib.loc = paste(basedir, "lib", sep = ""))
library(caret, lib.loc = paste(basedir, "lib", sep = ""))
'%ni%' <- Negate('%in%')  # define 'not in' func
options(scipen=999)  # prevents printing scientific notations.
set.seed(100)
index <- createDataPartition(df['diagnosis',], p=0.75, list = F)
trainset <- data_genefilter[index, ]
testset <- data_genefilter[-index, ]

### Logistic Regression
trainset[nrow(trainset) + 1, ] <- diagnosis
library(e1071, lib.loc = paste(basedir, "lib", sep = ""))
modellogit <- glm(col(trainset, as.factor = TRUE) ~ ., family = binomial(link = "logit"), data = trainset)
summary(modellogit)
print(cbind("coeff" = modellogit$coefficients, "odds ratio" = (exp(modellogit$coefficients) - 1) * 100)) # Odds ratio
pred.modellogit <- predict(modellogit, newdata = testset[, colnames(testset) %ni% "ADDEPEV3"], type = "response") > 0.5
confMatrix.modellogit <- table(testset$ADDEPEV3, pred.modellogit)
print(confMatrix.modellogit)
err_metric(confMatrix.modellogit[1,1], confMatrix.modellogit[2,2], confMatrix.modellogit[1,2], confMatrix.modellogit[2,1])
