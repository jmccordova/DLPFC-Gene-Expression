# Part 4: Machine Learning

# Part 3.1.4: Alpha computation
# Perform Bonferroni correction
get_alpha <- function(alpha_orig, n) {
  return(alpha_orig/n)
}
alpha <- 0.05
alpha <- get_alpha(alpha, ncol(data.binomial) - 1)

# Part 4.1: Naive Bayes
set.seed(100)
model.nb <- train(x = trainset.binomial[, colnames(trainset.binomial) != "diagnosis"], 
                  y = trainset.binomial$diagnosis, 
                  method = "naive_bayes", 
                  trControl = trainControl(method='cv', number=10)
)
pred.model.nb <- predict(model.nb, newdata = testset.binomial[, colnames(trainset.binomial) != "diagnosis"])
confMatrix.model.nb <- confusionMatrix(pred.model.nb, testset.binomial$diagnosis)
var.model.nb <- varImp(model.nb, useModel = TRUE, nonpara = TRUE, scale = TRUE)

# Part 4.2: K Nearest Neighbors
set.seed(100)
trControl.knn <- trainControl(method='repeatedcv', number = 3, allowParallel = TRUE)
trainset.binomial.preprocessed <- preProcess(trainset.binomial[, colnames(trainset.binomial) != "diagnosis"])
model.knn <- train(x = trainset.binomial[, colnames(trainset.binomial) != "diagnosis"], 
                   y = trainset.binomial$diagnosis, 
                   method = "knn", 
                   trControl = trControl.knn,
                   tuneLength = 20
)
pred.model.knn <- predict(model.knn, newdata = testset.binomial)
confMatrix.model.knn <- confusionMatrix(pred.model.knn, testset.binomial$diagnosis)
var.model.knn <- varImp(model.knn, useModel = TRUE, nonpara = TRUE, scale = TRUE)
# Part 4.3: Decision Tree
model.dt <- rpart(formula = diagnosis ~ .,
                  data = trainset.binomial,
                  method = "class", 
                  control = rpart.control(minsplit=2, minbucket = 1, cp = 0.001)
)
pred.model.dt <- predict(model.dt, newdata = testset.binomial)
confMatrix.model.dt <- confusionMatrix(pred.model.dt, testset.binomial$diagnosis)
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
    model.svm <- svm(x = trainset.binomial[, colnames(trainset.binomial) != "diagnosis"], 
                     y = trainset.binomial$diagnosis, 
                     kernel = kernel, 
                     cost = cost,
                     probability = TRUE,
                     scale = TRUE
    )
    pred.model.svm <- predict(model.svm, newdata = testset.binomial[, colnames(trainset.binomial) != "diagnosis"])
    confMatrix.model.svm <- confusionMatrix(pred.model.svm, testset.binomial$diagnosis)
    print(confMatrix.modelsvm)
    #err_metric(confMatrix.modelsvm$table[1,1], confMatrix.modelsvm$table[2,2], confMatrix.modelsvm$table[1,2], confMatrix.modelsvm$table[2,1])
    # modelsvm <- cbind(modelsvm, c(modelsvm))    
  }
}
# Part 4.5: Logistic Regression
#model.logit <- glm(x = trainset.binomial[, colnames(trainset.binomial) != "diagnosis"], 
#                   y = trainset.binomial$diagnosis,
#                   family = binomial(link = "logit"), 
#               )
model.logit <- multinom(diagnosis ~ .,
                        data = trainset.binomial)
pred.model.logit <- predict(model.logit, newdata = testset.binomial[, colnames(trainset.binomial) != "diagnosis"], type = "response") > 0.5
confMatrix.model.logit <- confusionMatrix(pred.model.logit, testset.binomial$diagnosis)
summary(model.logit)
print(cbind("coeff" = model.logit$coefficients, "odds ratio" = (exp(model.logit$coefficients) - 1) * 100)) # Odds ratio
# Part 4.6: Discriminant Analysis
# Part 4.7: Random Forest
model.rf <- randomForest(diagnosis ~ ., 
                         data = trainset.binomial, 
                         ntree = 300, 
                         mtry = 21591
)
pred.model.rf <- predict(model.rf, newdata = testset.binomial)
confMatrix.model.rf <- confusionMatrix(pred.model.rf, testset.binomial$diagnosis)