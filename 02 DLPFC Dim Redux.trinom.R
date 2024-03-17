# Part 2.5: Create a metadata array to be used for PCA
data.pca.metadata <- matrix(c(transcriptIDs$id_internal_huex, transcriptIDs$transcript_id), nrow = 2, ncol = nrow(transcriptIDs), byrow = TRUE)
data.pca.metadata <- data.frame(data.pca.metadata)
names(data.pca.metadata) <- data.pca.metadata[1, ]
data.pca.metadata <- data.pca.metadata[-1, ]

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

# Part 3.1.3: Get gene expressions
gene_expressions <- exprs(e)
data <- gene_expressions
data <- t(data)
data <- as.data.frame(data)
print(paste("After RMA: ", paste(c("Samples: ", "Expressions: "), dim(data), collapse = " ")))
# Append diagnosis.fact factors to the dataset
data.multinomial <- data
data.multinomial[, 'diagnosis'] <- diagnosis.fact

# Part 3.1.4: Alpha computation
# Perform Bonferroni correction
get_alpha <- function(alpha_orig, n) {
  return(alpha_orig/n)
}
alpha <- 0.05
alpha <- get_alpha(alpha, ncol(data))

# Part 3.1.5: Dimension reduction using Principal Component Analysis
# Step 3.1.5.1: Check for null values (If it returned more than 0, there is a null value)
print(colSums(is.na(data))[colSums(is.na(data)) != 0])
# Step 3.1.5.2: Normalize data
# data.normalized <- scale(data)
data.normalized <- apply(data, 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
# Step 3.1.5.3: Compute correlation matrix
corr_matrix <- cor(data.normalized)
# Step 3.1.5.3.1: Visualize the correlation matrix
qgraph(corr_matrix, minimum = 0.25, cut = 0.4, vsize = 2, legend = TRUE, borders = FALSE)
# Step 3.1.5.4: Computer covariance matrix
cov_matrix <- cov(data.normalized)
# Step 3.1.5.4.1: Check if there are NaN in covariance matrix
print(any(is.nan(cov_matrix)))
# Step 3.1.5.5: Apply PCA
# Step 3.1.5.5.1: Transpose the data first so the probe IDs are the rows
data.t <- data.frame(t(data))
#rownames(data.t) <- transcriptIDs$transcript_id[match(rownames(data.t), transcriptIDs$id_internal_huex)]
# Step 3.1.5.5.2: Check if the data rows and metadata columns match names
all(rownames(data.t) %in% colnames(data.pca.metadata))
all(rownames(data.t) == colnames(data.pca.metadata))
# Step 3.1.5.5.3: Perform PCA
#data.pca <- pca(data.t, metadata = data.pca.metadata, removeVar = 0.1)
data.pca <- pca(data.t, removeVar = 0.1)
# Determine optimum number of PCs using Horn's method
horn <- parallelPCA(data.t)
# Find the elbow point
elbow <- findElbowPoint(data.pca$variance)
# Get the PCs with that can explain at least 50% of the variability
print(paste('Number of PCs to consider:'))
print(paste('Horn:', horn$n, '%Varation:', sum(data.pca$variance[1:horn$n])))
print(paste('Elbow:', elbow, '%Variation:', sum(data.pca$variance[1:elbow])))
data.important <- data.pca$loadings[, 1:horn$n]
write.csv(data.important, paste(datadir, "data.important.csv", sep = ""), row.names = TRUE)
# Show the loadings
# plotloadings(data.pca, labSize = 3)
plotloadings(data.pca,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)
plotloadings(data.pca,
             components = getComponents(data.pca, 1:elbow),
             rangeRetain = 0.1,
             labSize = 4.0,
             absolute = FALSE,
             title = 'Loadings plot',
             subtitle = 'Misc PCs',
             caption = 'Top 10% variables',
             shape = 23, shapeSizeRange = c(1, 16),
             col = c('white', 'pink'),
             drawConnectors = FALSE)
data.pca.vars <- getVars(data.pca)
chooseMarchenkoPastur(data.t, var.explained=data.pca$sdev^2, noise=4)
screeplot(data.pca,
  components = getComponents(data.pca, 1:20),
  vline = c(horn$n, elbow)) +
  geom_label(aes(x = horn$n + 1, y = 50,
    label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
    label = 'Elbow method', vjust = -1, size = 8)
)

biplot(data.pca, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
pairsplot(data.pca)

# Part 3.1.5: Splitting dataset
## For multinomial
'%ni%' <- Negate('%in%')  # define 'not in' func
options(scipen=999)  # prevents printing scientific notations.
set.seed(100)
index.multinomial <- createDataPartition(data.multinomial$diagnosis, p=0.75, list = F)
trainset.multinomial <- data.multinomial[index.multinomial, ]
testset.multinomial <- data.multinomial[-index.multinomial, ]

# I.2.2. Data split
# options(scipen=999)  # prevents printing scientific notations.
# set.seed(100)
# index <- createDataPartition(data_genefilter$withDisorder, p=0.75, list = F)
# trainset <- data_genefilter[index, ]
# testset <- data_genefilter[-index, ]

