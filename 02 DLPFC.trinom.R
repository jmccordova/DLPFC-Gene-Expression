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

# Part 4: Machine Learning
# Part 4.1: Naive Bayes
set.seed(100)
model.nb <- train(x = trainset.multinomial[, colnames(trainset.multinomial) != "diagnosis"], 
                  y = trainset.multinomial$diagnosis, 
                  method = "naive_bayes", 
                  trControl = trainControl(method='cv', number=10)
            )
pred.model.nb <- predict(model.nb, newdata = testset.multinomial[, colnames(trainset.multinomial) != "diagnosis"])
confMatrix.model.nb <- confusionMatrix(pred.model.nb, testset.multinomial$diagnosis)
var.model.nb <- varImp(model.nb, useModel = TRUE, nonpara = TRUE, scale = TRUE)

# Part 4.2: K Nearest Neighbors
set.seed(100)
trControl.knn <- trainControl(method='repeatedcv', number = 3, allowParallel = TRUE)
trainset.multinomial.preprocessed <- preProcess(trainset.multinomial[, colnames(trainset.multinomial) != "diagnosis"])
model.knn <- train(x = trainset.multinomial[, colnames(trainset.multinomial) != "diagnosis"], 
                  y = trainset.multinomial$diagnosis, 
                  method = "knn", 
                  trControl = trControl.knn,
                  tuneLength = 20
)
pred.model.knn <- predict(model.knn, newdata = testset.multinomial[, colnames(trainset.multinomial) != "diagnosis"], type = "class")
confMatrix.model.knn <- confusionMatrix(pred.model.knn, testset.multinomial$diagnosis)
var.model.knn <- varImp(model.knn, useModel = TRUE, nonpara = TRUE, scale = TRUE)
# Part 4.3: Decision Tree
model.dt <- rpart(formula = diagnosis ~ .,
                  data = trainset.multinomial,
                  method = "class", 
                  control = rpart.control(minsplit=2, minbucket = 1, cp = 0.001)
            )
pred.model.dt <- predict(model.dt, newdata = testset.multinomial, type = "class")
confMatrix.model.dt <- confusionMatrix(pred.model.dt, testset.multinomial$diagnosis)
var.model.dt<- varImp(model.dt, useModel = TRUE, nonpara = TRUE, scale = TRUE)
summary(model.dt)
rpart.plot(model.dt)
# Part 4.4: SVM
# For SVM and random forest, cut the dataset to 10% of the dataset to make processing quicker
#trainset.cut <- trainset[sample(x = 1:nrow(trainset), size = nrow(trainset) * .10, replace = TRUE), colnames(trainset)] 
#trainset.cut <- upSample(x = trainset.cut[, colnames(trainset.cut) %ni% "ADDEPEV3"], yname = "ADDEPEV3", y = trainset.cut$ADDEPEV3)
for (kernel in c("linear", "polynomial", "radial", "sigmoid")) {
  for (cost in c(0.001, 0.01, 0.1, 1, 5, 10, 100)) {
    print(paste(kernel," @ ", cost))
    model.svm <- svm(x = trainset.multinomial[, colnames(trainset.multinomial) != "diagnosis"], 
                     y = trainset.multinomial$diagnosis, 
                     kernel = kernel, 
                     cost = cost, 
                     scale = TRUE
                  )
    pred.model.svm <- predict(model.svm, testset.multinomial)
    confMatrix.modelsvm <- confusionMatrix(pred.model.svm.testset.multinomial$diagnosis)
    #err_metric(confMatrix.modelsvm$table[1,1], confMatrix.modelsvm$table[2,2], confMatrix.modelsvm$table[1,2], confMatrix.modelsvm$table[2,1])
    # modelsvm <- cbind(modelsvm, c(modelsvm))    
  }
}