# Part 3: Exploration
# Part 3.1: Turn probe-level information to gene-level data by normalization using RMA
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

# Part 3.1.3: Get gene expressions
gene_expressions <- exprs(e)
data <- gene_expressions
remove(gene_expressions, e)

data <- t(data)
data <- as.data.frame(data)
print(paste("After RMA: ", paste(c("Samples: ", "Expressions: "), dim(data), collapse = " ")))
# Append diagnosis.fact factors to the dataset
data.multinomial <- data
data.multinomial[, 'diagnosis'] <- diagnosis.fact
