# Part 3: Analysis

# Part 3.1.4: Alpha computation
# Perform Bonferroni correction
get_alpha <- function(alpha_orig, n) {
  return(alpha_orig/n)
}
alpha <- 0.05
alpha <- get_alpha(alpha, ncol(data))

# Part 3.1.5: Show correlation matrix
#varcorr.matrix <- rcorr(as.matrix(data))
#ggcorrplot(varcorr.matrix$r, method = "square", type = "lower", lab = TRUE, p.mat = varcorr.matrix$P)

# Part 3.1.5: Splitting dataset
options(scipen=999)  # prevents printing scientific notations.
set.seed(100)
index.multinomial <- createDataPartition(data.multinomial$diagnosis, p=0.75, list = F)
trainset.multinomial <- data.multinomial[index.multinomial, ]
testset.multinomial <- data.multinomial[-index.multinomial, ]
# Perform upsampling so that the predictors are on equal footing
trainset.multinomial <- upSample(x = trainset.multinomial[, colnames(trainset.multinomial) != "diagnosis"],
         y = trainset.multinomial$diagnosis,
         yname = "diagnosis")
# Show the percentage of each predictor
prop.table(table(trainset.multinomial$diagnosis)) * 100
levels(trainset.multinomial$diagnosis) <- levels(trainset.multinomial$diagnosis)
levels(testset.multinomial$diagnosis) <- levels(testset.multinomial$diagnosis)