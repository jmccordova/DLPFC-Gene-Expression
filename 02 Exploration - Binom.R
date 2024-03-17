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
# CTL - 9, BPD, MDD, SCZ - 1
diagnosis <- c(
  # 299-309
  9, 1, 1, 1, 1, 1, 1, 1, 1, 9, 1,
  # 310-320
  1, 1, 1, 9, 1, 1, 1, 1, 1, 9, 1,
  # 321-331
  1, 9, 1, 9, 1, 9, 9, 9, 1, 1, 9,
  # 332-342
  1, 9, 1, 1, 9, 1, 1, 9, 1, 1, 1,
  # 343-353
  9, 1, 9, 1, 9, 1, 1, 1, 1, 9, 1,
  # 354-364
  9, 9, 1, 9, 1, 1, 1, 9, 1, 9, 1,
  # 365-375
  1, 9, 9, 1, 9, 9, 9, 1, 1, 9, 1,
  # 376-386
  9, 1, 1, 9, 1, 1, 1, 1, 1, 9, 1,
  # 387-397
  9, 1, 1, 9, 1, 9, 1, 9, 1, 1, 1,
  # 398-408
  9, 9, 1, 1, 1, 9, 9, 9, 9, 1, 1,
  # 409-419
  9, 9, 1, 1, 9, 1, 1, 1, 1, 1, 1,
  # 420-430
  9, 1, 9, 1, 1, 1, 1, 1, 1, 1, 1,
  # 431-441
  1, 1, 1, 1, 1, 1, 1, 1, 9, 1, 1,
  # 442-452
  9, 9, 1, 9, 9, 1, 1, 9, 1, 9, 1,
  # 453-463
  1, 1, 9, 1, 9, 1, 9, 9, 9, 9, 9,
  # 464-476
  9, 9, 1, 1
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
