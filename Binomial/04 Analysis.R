# Part 4: Machine Learning
exportsubdir <- "Step 4 - Analysis"
dir.create(paste(exportdir, exportsubdir, sep = "/"), recursive=TRUE)
  # Part 4.1: Compute the alpha using Bonferroni
  alpha <- get_alpha(alpha, nrow(data.binomial) - 1)
  # Part 4.2: Create a function for each ML
  perform_learning <- function(method, trainset, testset, 
                               svm.kernel = NULL, 
                               svm.cost = NULL,
                               rf.ntree = NULL,
                               rf.mtry = NULL,
                               export.filename = NULL,
                               tune = FALSE) {
    if (method == "NB") {
      # Part 4.1: Naive Bayes
      set.seed(100)
      model.nb <- train(x = trainset[, colnames(trainset) != "diagnosis"], 
                        y = trainset$diagnosis, 
                        method = "naive_bayes", 
                        trControl = trainControl(method='cv', number=10)
      )
      pred.model.nb <- predict(model.nb, newdata = testset[, colnames(trainset) != "diagnosis"])
      roc.model.nb <- auc(actual = testset$diagnosis, 
                          predicted = predict(model.nb, 
                                             newdata = testset[, colnames(testset) != "diagnosis"], 
                                             type = "raw"
                                      )
                          )
      confMatrix.model.nb <- confusionMatrix(pred.model.nb, testset$diagnosis)
      #var.model.nb <- varImp(model.nb, useModel = TRUE, nonpara = TRUE, scale = TRUE)
      var.model.nb <- varImp(model.nb, scale = FALSE)
      return(list(model = model.nb, pred = pred.model.nb, confMatrix = confMatrix.model.nb, var = var.model.nb, roc = roc.model.nb))
    } else if (method == "KNN") {
      # Part 4.2: K Nearest Neighbors
      set.seed(100)
      trControl.knn <- trainControl(method='repeatedcv', number = 3, allowParallel = TRUE)
      trainset.preprocessed <- preProcess(trainset[, colnames(trainset) != "diagnosis"])
      model.knn <- train(x = trainset[, colnames(trainset) != "diagnosis"], 
                         y = trainset$diagnosis, 
                         method = "knn", 
                         trControl = trControl.knn,
                         tuneLength = 20
      )
      pred.model.knn <- predict(model.knn, newdata = testset)
      roc.model.knn <- auc(actual = testset$diagnosis, 
                          predicted = predict(model.knn, 
                                              newdata = testset[, colnames(testset) != "diagnosis"], 
                                              type = "raw"
                                      )
                        )
      confMatrix.model.knn <- confusionMatrix(pred.model.knn, testset$diagnosis)
      #var.model.knn <- varImp(model.knn, useModel = TRUE, nonpara = TRUE, scale = TRUE)
      var.model.knn <- varImp(model.knn, scale = FALSE)
      return(list(model = model.knn, pred = pred.model.knn, confMatrix = confMatrix.model.knn, var = var.model.knn, roc = roc.model.knn))
    } else if (method == "DT") {
      # Part 4.3: Decision Tree
      # Part 4.3.1: Using rpart
      model.dt <- rpart(formula = diagnosis ~ .,
                        data = trainset,
                        method = "class", 
                        control = rpart.control(minsplit=2, minbucket = 1, cp = 0.001)
      )
      #var.model.dt<- varImp(model.dt, useModel = TRUE, nonpara = TRUE, scale = TRUE)
      var.model.dt<- varImp(model.dt, scale = FALSE)
      arrange(var.model.dt, desc(Overall))
      rpart.plot(model.dt)
      pdf(export.filename)
      prp(model.dt, extra=104)
      dev.off()
      
      # Part 4.3.2: Using train
      ## 10-fold CV
      ## repeated ten times
      trControl.dt <- trainControl(
        method = "repeatedcv",
        number = 5,
        repeats = 10)
      model.dt <- train(x = trainset[, colnames(trainset) != "diagnosis"], 
                        y = trainset$diagnosis, 
                        method = "rpart2", 
                        trControl = trControl.dt
      )
      pred.model.dt <- predict(model.dt, newdata = testset)
      roc.model.dt <- auc(actual = testset$diagnosis, 
                          predicted = predict(model.dt, 
                                              newdata = testset[, colnames(testset) != "diagnosis"], 
                                              type = "raw"
                          )
                      )
      confMatrix.model.dt <- confusionMatrix(pred.model.dt, testset$diagnosis)
      return(list(model = model.dt, pred = pred.model.dt, confMatrix = confMatrix.model.dt, var = var.model.dt, roc = roc.model.dt))
    } else if (method == "SVM") {
      # Part 4.4: SVM
      # For SVM and random forest, cut the dataset to 10% of the dataset to make processing quicker
      #trainset.cut <- trainset[sample(x = 1:nrow(trainset), size = nrow(trainset) * .10, replace = TRUE), colnames(trainset)] 
      #trainset.cut <- upSample(x = trainset.cut[, colnames(trainset.cut) %ni% "ADDEPEV3"], yname = "ADDEPEV3", y = trainset.cut$ADDEPEV3)
      model.svm <- rminer::fit(diagnosis ~ ., 
                               data = trainset, 
                               model = "ksvm",
                               kpar = "automatic",
                               task = "class",
                               search = list(search = mparheuristic("ksvm", n = 5))
      )
      pred.model.svm <- predict(model.svm, newdata = testset)
      confMatrix.model.svm <- confusionMatrix(pred.model.svm, testset$diagnosis)
      if (tune) {
        print(paste(kernel," @ ", cost))
        print(confMatrix.model.svm)
      }
      var.model.svm <- Importance(model.svm, data = trainset, method = "DSA", outindex = "diagnosis")
      var.model.svm <- data.frame("feature" = colnames(var.model.svm$data), "Overall" = var.model.svm$imp)
      var.model.svm <- head(var.model.svm, -1)
      var.model.svm <- var.model.svm[order(-var.model.svm$Overall),]
      if (!tune) {
        roc.model.svm <- auc(actual = testset$diagnosis, 
                              predicted = predict(model.svm, 
                                                  newdata = testset[, colnames(testset) != "diagnosis"], 
                                                  type = "raw"
                              )
                          )
      }
      
      if (!tune) {
        return(list(model = model.svm, pred = pred.model.svm, confMatrix = confMatrix.model.svm, var = var.model.svm, roc = roc.model.svm))
      }
    } else if (method == "LOG") {
      model.logit <- multinom(diagnosis ~ .,
                              data = trainset)
      pred.model.logit <- predict(model.logit, newdata = testset[, colnames(testset) != "diagnosis"], type = "class")
      roc.model.logit <- multiclass.roc(response = testset$diagnosis, 
                                        predictor = predict(model.logit, 
                                                            newdata = testset[, colnames(testset) != "diagnosis"], 
                                                            type = "prob"),
                                        percent = TRUE)
      roc.model.logit <- auc(actual = testset$diagnosis, 
                             predicted = predict(model.logit, 
                                                 newdata = testset[, colnames(testset) != "diagnosis"], 
                                                 type = "probs"
                             )
                          )
      confMatrix.model.logit <- confusionMatrix(pred.model.logit, testset$diagnosis)
      #var.model.logit <- varImp(model.logit, useModel = TRUE, nonpara = TRUE, scale = TRUE)
      var.model.logit <- varImp(model.logit, scale = FALSE)
      var.model.logit <- data.frame("feature" = rownames(var.model.logit), "Overall" = var.model.logit$Overall)
      var.model.logit <- var.model.logit[order(-var.model.logit$Overall), ]
      
      #summary(model.logit)
      #View(cbind("coeff" = coef(model.logit), "odds ratio" = (exp(coef(model.logit)) - 1) * 100)) # Odds ratio
      
      return(list(model = model.logit, pred = pred.model.logit, confMatrix = confMatrix.model.logit, var = var.model.logit, roc = roc.model.logit))
    } else if (method == "DA") {
      # Part 4.6: Discriminant Analysis
      trControl.lda <- trainControl(classProbs = TRUE)
      levels(trainset$diagnosis)[match("1",levels(trainset$diagnosis))] <- "DIS"
      levels(trainset$diagnosis)[match("9",levels(trainset$diagnosis))] <- "CTL"
      model.lda <- train(diagnosis ~ ., 
                         data = trainset, 
                         method = "lda",
                         trControl = trControl.lda)
      levels(testset$diagnosis)[match("1",levels(testset$diagnosis))] <- "DIS"
      levels(testset$diagnosis)[match("9",levels(testset$diagnosis))] <- "CTL"
      var.model.lda <- varImp(model.lda, useModel = TRUE, nonpara = TRUE, scale = TRUE)
      
      model.lda <- lda(diagnosis ~ ., 
                       data = trainset,
      )
      pred.model.lda <- predict(model.lda, newdata = testset, type = "response")
      roc.model.lda <- auc(testset$diagnosis, pred.model.lda$class)
      confMatrix.model.lda <- confusionMatrix(pred.model.lda$class, testset$diagnosis)
      #return(list(model = model.lda, pred = pred.model.lda, confMatrix = confMatrix.model.lda, var = var.model.lda))
      return(list(model = model.lda, pred = pred.model.lda, confMatrix = confMatrix.model.lda, var = var.model.lda, roc = roc.model.lda))
    } else if (method == "RF") {
      # Part 4.7: Random Forest
      set.seed(100)
      # tuneRF 
      invisible(capture.output(fgl.res <- tuneRF(x = trainset[, colnames(trainset) != "diagnosis"],
                                                 y = trainset$diagnosis, 
                                                 stepFactor=1.5)
                               )
                )
      
      # choose the best mtry based on the lowest OOB error
      best_mtry <- fgl.res[fgl.res[, 2] == min(fgl.res[, 2]), 1]
      # choose the lowest OOB error
      best_oob  <- fgl.res[fgl.res[, 2] == min(fgl.res[, 2]), 2]
      
      if (tune) {
        ntrees <- seq.int(3001, 4001,  by = 100)
      } else {
        ntrees <- c(rf.ntree)
      }
      
      mtries <- c(best_mtry)
      
      
      for(ntree in ntrees) {
        for(mtry in mtries) {
          model.rf <- randomForest(x = trainset[, colnames(trainset) != "diagnosis"],
                                   y = trainset$diagnosis, 
                                   ntree = ntree, 
                                   mtry = mtry,
                                   oob.error=best_oob
          )
          pred.model.rf <- predict(model.rf, newdata = testset)
          if (!tune) {
            roc.model.rf <- auc(actual = testset$diagnosis, 
                                predicted = predict(model.rf, 
                                                    newdata = testset[, colnames(testset) != "diagnosis"], 
                                                    type = "response"
                                )
                            )
          }
          confMatrix.model.rf <- confusionMatrix(pred.model.rf, testset$diagnosis)
          if (tune) {
            print(paste(ntree," and ", mtry))
            print(confMatrix.model.rf)
          }
          var.model.rf <- varImp(model.rf, useModel = TRUE, nonpara = TRUE, scale = TRUE)
          var.model.rf <- arrange(var.model.rf, desc(Overall))
        }
      }
      
      if (!tune) {
        return(list(model = model.rf, pred = pred.model.rf, confMatrix = confMatrix.model.rf, var = var.model.rf, roc = roc.model.rf))
      }
    } else {
      models <- c("naive","ctree","cv.glmnet","rpart","kknn","ksvm","lssvm","mlp","mlpe", "randomForest","lda","multinom", "naiveBayes","xgboost", "SE")
      if (method == "SOFT") {
        models <- c(models, "SE")
      }
      
      model.auto <- rminer::fit(diagnosis ~ ., 
                                data = trainset, 
                                model = "auto",
                                fdebug = TRUE,
                                feature = "sabs",
                                scale = "inputs",
                                search = list(
                                  search = mparheuristic( 
                                    model = models,
                                    task = "class", 
                                    inputs = ncol(trainset)-1
                                  ),
                                  smethod = "auto",
                                  metric = "AUC",
                                  convex = 0
                                )
      )
      
      pred.model.auto <- predict(model.auto, testset)
      roc.model.auto <- round(mmetric(testset$diagnosis, pred.model.auto, metric="AUC"), 2)
      importanceMethod <- "sens"
      if (model.auto@model == "randomForest") {
        importanceMethod <- "randomForest"
      }
      
      var.model.auto <- Importance(M = model.auto, 
                                   #method = importanceMethod,
                                   #outindex = "diagnosis",
                                   data = trainset
      )
      # Create dataframe that combines feature names and the variable importance
      var.model.auto <- data.frame("feature" = colnames(trainset), "Overall" = var.model.auto$imp)
      # Remove diagnosis row
      var.model.auto <- head(var.model.auto, -1)
      # Remove 0 importance
      var.model.auto <- subset(var.model.auto, Overall != 0)
      # Reorder from highest importance to least
      var.model.auto <- var.model.auto[order(-var.model.auto$Overall),]
      
      # show leaderboard:
      cat("Models  by rank:", model.auto@mpar$LB$model, "\n")
      cat("Validation values:", round(model.auto@mpar$LB$eval,4), "\n")
      cat("Best model:", model.auto@model, "\n")
      cat("AUC", "=", roc.model.auto, "\n")
      
      return(list(model = model.auto, pred = pred.model.auto, confMatrix = rminer::mmetric(testset$diagnosis, pred.model.auto, "ALL"), var = var.model.auto, roc = roc.model.auto))
    }
  }

  # Part 4.3: Perform ML
    # Part 4.3.1: Gene Filtering Dataset
    data.binomial <- createDataset(dataSource = data.pp, feature = features.gf, probeset = huex.probes, filename = "Gene Filtering")
    sets <- buildTrainTest(data.binomial)
    trainset.binomial <- sets$trainset
    testset.binomial <- sets$testset
    remove(sets)
      # Part 4.3.1.1: Perform tuning for SVM and Random Forest
      perform_learning("SVM", trainset.binomial, testset.binomial, tune = TRUE)
      perform_learning("RF", trainset.binomial, testset.binomial, tune = TRUE)
      # Part 4.3.1.2: Perform analysis
        # Part 4.3.1.2.1: Do ensemble of all learning methods
        learn.gf.auto <- perform_learning("auto", trainset.binomial, testset.binomial)
        # Part 4.3.1.2.2: Naive Bayes
        learn.gf.nb <- perform_learning("NB", trainset.binomial, testset.binomial)
        # Part 4.3.1.2.3: KNN
        learn.gf.knn <- perform_learning("KNN", trainset.binomial, testset.binomial)
        # Part 4.3.1.2.4: SVM 
        learn.gf.svm <- perform_learning("SVM", trainset.binomial, testset.binomial, svm.kernel = 'splinedot', svm.cost = 100)
        # Part 4.3.1.2.5: Logistic Regression
        learn.gf.log <- perform_learning("LOG", trainset.binomial, testset.binomial)
        # Part 4.3.1.2.6: Discriminant Analysis
        learn.gf.da <- perform_learning("DA", trainset.binomial, testset.binomial)
        # Part 4.3.1.2.7: Decision Tree
        learn.gf.dt <- perform_learning("DT", trainset.binomial, testset.binomial, export.filename = paste(exportdir, exportsubdir, "Decision Tree (Gene Filter).pdf", sep = "/"))
        # Part 4.3.1.2.8: Random Forest 
        learn.gf.rf <- perform_learning("RF", trainset.binomial, testset.binomial, rf.ntree = 100001)
      # Part 4.3.2: PCA Dataset
      data.binomial <- createDataset(dataSource = data.pp, feature = features.pca, probeset = huex.probes, filename = "PCA")
      sets <- buildTrainTest(data.binomial)
      trainset.binomial <- sets$trainset
      testset.binomial <- sets$testset
      remove(sets)
        # Part 4.3.2.1: Perform tuning for SVM and Random Forest
        perform_learning("SVM", trainset.binomial, testset.binomial, tune = TRUE)
        perform_learning("RF", trainset.binomial, testset.binomial, tune = TRUE)
        # Part 4.3.2.2: Perform analysis
          # Part 4.3.2.2.1: Do ensemble of all learning methods
          learn.pca.auto <- perform_learning("auto", trainset.binomial, testset.binomial)
          # Part 4.3.2.2.2: Naive Bayes
          learn.pca.nb <- perform_learning("NB", trainset.binomial, testset.binomial)
          # Part 4.3.2.2.3: KNN
          learn.pca.knn <- perform_learning("KNN", trainset.binomial, testset.binomial)
          # Part 4.3.2.2.4: SVM (For PCA, SVM tuning had no significant features)
          learn.pca.svm <- perform_learning("SVM", trainset.binomial, testset.binomial, svm.kernel = 'splinedot', svm.cost = 100)
          # Part 4.3.2.2.5: Logistic Regression
          learn.pca.log <- perform_learning("LOG", trainset.binomial, testset.binomial)
          # Part 4.3.2.2.6: Discriminant Analysis
          learn.pca.da <- perform_learning("DA", trainset.binomial, testset.binomial)
          # Part 4.3.2.2.7: Decision Tree
          learn.pca.dt <- perform_learning("DT", trainset.binomial, testset.binomial, export.filename = paste(exportdir, exportsubdir, "Decision Tree (PCA).pdf", sep = "/"))
          # Part 4.3.2.2.8: Random Forest  (For PCA, SVM tuning had no significant features)
          learn.pca.rf <- perform_learning("RF", trainset.binomial, testset.binomial, rf.ntree = 3501, rf.mtry = 5)
    # Part 4.3.3: Combined Dataset
    data.binomial <- createDataset(dataSource = data.pp, feature = features, probeset = huex.probes, filename = "PCA + GF")
    sets <- buildTrainTest(data.binomial)
    trainset.binomial <- sets$trainset
    testset.binomial <- sets$testset
    sets <- buildTrainTest(data.binomial)
    validationset.binomial <- sets$validationset
    remove(sets)
      # Part 4.3.3.1: Perform tuning for SVM and Random Forest
      perform_learning("SVM", trainset.binomial, testset.binomial, tune = TRUE)
      perform_learning("RF", trainset.binomial, testset.binomial, tune = TRUE)
      # Part 4.3.3.2: Perform analysis
        # Part 4.3.3.2.1: Do ensemble of all learning methods
        learn.features.auto <- perform_learning("auto", trainset.binomial, testset.binomial)
        # Part 4.3.3.2.2: Naive Bayes
        learn.features.nb <- perform_learning("NB", trainset.binomial, testset.binomial)
        # Part 4.3.3.2.3: KNN
        learn.features.knn <- perform_learning("KNN", trainset.binomial, testset.binomial)
        # Part 4.3.3.2.4: SVM
        learn.features.svm <- perform_learning("SVM", trainset.binomial, testset.binomial, svm.kernel = 'rbfdot', svm.cost = 0.0170484)
        # Part 4.3.3.2.5: Logistic Regression
        learn.features.log <- perform_learning("LOG", trainset.binomial, testset.binomial)
        # Part 4.3.3.2.6: Discriminant Analysis
        learn.features.da <- perform_learning("DA", trainset.binomial, testset.binomial)
        # Part 4.3.3.2.7: Decision Tree
        learn.features.dt <- perform_learning("DT", trainset.binomial, testset.binomial, export.filename = paste(exportdir, exportsubdir, "Decision Tree (GF + PCA).pdf", sep = "/"))
        # Part 4.3.3.2.8: Random Forest  (For PCA, SVM tuning had no significant features)
        learn.features.rf <- perform_learning("RF", trainset.binomial, testset.binomial, rf.ntree = 3201, rf.mtry = 13)
        # Part 4.3.3.2.9: Soft voting from all the models
        learn.features.soft <- perform_learning("SOFT", trainset.binomial, testset.binomial)
          # Add diagnosis in the important features
          features.important <- c(learn.features.soft$var$feature, "diagnosis")
          # Remove the unimportant features
          data.binomial.important <- trainset.binomial[colnames(trainset.binomial) %in% features.important]
          # Rebuild the dataset
          alpha <- get_alpha(alpha, nrow(data.binomial.important) - 1)
          sets <- buildTrainTest(data.binomial.important)
          trainset.binomial.important <- sets$trainset
          testset.binomial.important <- sets$testset
          sets <- buildTrainTest(data.binomial.important)
          validationset.binomial.important <- sets$validationset
          remove(sets)
        