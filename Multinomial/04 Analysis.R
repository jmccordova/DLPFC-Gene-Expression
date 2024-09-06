# Part 4: Machine Learning
exportsubdir <- "Step 4 - Analysis"
dir.create(paste(exportdir, exportsubdir, sep = "/"), recursive=TRUE)
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
      roc.model.nb <- multiclass.roc(response = testset$diagnosis, 
                                     predictor = predict(model.nb, 
                                                         newdata = testset[, colnames(testset) != "diagnosis"], 
                                                         type = "prob"),
                                     percent = TRUE)
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
      roc.model.knn <- multiclass.roc(response = testset$diagnosis, 
                                     predictor = predict(model.knn, 
                                                         newdata = testset[, colnames(testset) != "diagnosis"], 
                                                         type = "prob"),
                                     percent = TRUE)
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
        var.model.dt <- varImp(model.dt, scale = FALSE)
        var.model.dt <- data.frame("feature" = rownames(var.model.dt), "Overall" = var.model.dt$Overall)
        var.model.dt <- var.model.dt[var.model.dt$Overall > 0, ]
        var.model.dt <- var.model.dt[order(-var.model.dt$Overall), ]
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
        roc.model.dt <- multiclass.roc(response = testset$diagnosis, 
                                        predictor = predict(model.dt, 
                                                            newdata = testset[, colnames(testset) != "diagnosis"], 
                                                            type = "prob"),
                                        percent = TRUE)
        confMatrix.model.dt <- confusionMatrix(pred.model.dt, testset$diagnosis)
        return(list(model = model.dt, pred = pred.model.dt, confMatrix = confMatrix.model.dt, var = var.model.dt, roc = roc.model.dt))
    } else if (method == "SVM") {
      # Part 4.4: SVM
      # For SVM and random forest, cut the dataset to 10% of the dataset to make processing quicker
      #trainset.cut <- trainset[sample(x = 1:nrow(trainset), size = nrow(trainset) * .10, replace = TRUE), colnames(trainset)] 
      #trainset.cut <- upSample(x = trainset.cut[, colnames(trainset.cut) %ni% "ADDEPEV3"], yname = "ADDEPEV3", y = trainset.cut$ADDEPEV3)
      if (tune) {
        kernels <- c("rbfdot", "polydot", "tanhdot", "vanilladot", "laplacedot", "besseldot", "anovadot", "splinedot")
        costs <- c(0.001, 0.01, 0.1, 1, 5, 10, 100)
      } else {
        kernels <- c(svm.kernel)
        costs <- c(svm.cost)
      }
      
      for (kernel in kernels) {
        for (cost in costs) {
          model.svm <- rminer::fit(diagnosis ~ ., 
                         data = trainset, 
                         model = "svm",
                         kernel = kernel,
                         kpar = "automatic", 
                         C = cost,
                         task = "class"
                       )
          pred.model.svm <- predict(model.svm, newdata = testset)
          confMatrix.model.svm <- confusionMatrix(pred.model.svm, testset$diagnosis)
          var.model.svm <- Importance(model.svm, data = trainset, method = "DSA", outindex = "diagnosis")
          var.model.svm <- data.frame("feature" = colnames(var.model.svm$data), "Overall" = var.model.svm$imp)
          var.model.svm <- head(var.model.svm, -1)
          var.model.svm <- var.model.svm[order(-var.model.svm$Overall),]
          
          #if (!tune) {
          roc.model.svm <- multiclass.roc(response = testset$diagnosis, 
                                          predictor = predict(rminer::fit(diagnosis ~ ., 
                                                                          data = trainset, 
                                                                          model = "svm",
                                                                          kernel = kernel,
                                                                          kpar = "automatic", 
                                                                          C = cost
                                                              ), 
                                                              newdata = testset[, colnames(testset) != "diagnosis"], 
                                                              type = "prob"),
                                          percent = TRUE)
          #}
          
          
          #if (tune && confMatrix.model.svm$overall['AccuracyPValue'] < 0.05) {
          if (tune) {
            print(paste(kernel," @ ", cost))
            print(confMatrix.model.svm$overall['Accuracy'])
            print(confMatrix.model.svm$overall['AccuracyPValue'])
            print(roc.model.svm$auc)
          }
        }
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
      levels(trainset$diagnosis)[match("1",levels(trainset$diagnosis))] <- "BPD"
      levels(trainset$diagnosis)[match("2",levels(trainset$diagnosis))] <- "MDD"
      levels(trainset$diagnosis)[match("3",levels(trainset$diagnosis))] <- "SCZ"
      levels(trainset$diagnosis)[match("9",levels(trainset$diagnosis))] <- "CTL"
      model.lda <- train(diagnosis ~ ., 
                         data = trainset, 
                         method = "mda",
                         trControl = trControl.lda)
      levels(testset$diagnosis)[match("1",levels(testset$diagnosis))] <- "BPD"
      levels(testset$diagnosis)[match("2",levels(testset$diagnosis))] <- "MDD"
      levels(testset$diagnosis)[match("3",levels(testset$diagnosis))] <- "SCZ"
      levels(testset$diagnosis)[match("9",levels(testset$diagnosis))] <- "CTL"
      #var.model.lda <- varImp(model.lda, useModel = TRUE, nonpara = TRUE, scale = TRUE)
      var.model.lda <- varImp(model.lda, scale = FALSE)
      
      model.lda <- lda(diagnosis ~ ., 
                       data = trainset,
      )
      pred.model.lda <- predict(model.lda, newdata = testset)
      roc.model.lda <- multiclass.roc(response = testset$diagnosis, 
                                      predictor = predict(model.lda, 
                                                          newdata = testset[, colnames(testset) != "diagnosis"], 
                                                          type = "prob")$posterior,
                                      percent = TRUE)
      confMatrix.model.lda <- confusionMatrix(pred.model.lda$class, testset$diagnosis)
      #return(list(model = model.lda, pred = pred.model.lda, confMatrix = confMatrix.model.lda, var = var.model.lda))
      return(list(model = model.lda, pred = pred.model.lda, confMatrix = confMatrix.model.lda, var = var.model.lda, roc = roc.model.lda))
    } else if (method == "RF") {
      # Part 4.7: Random Forest
      set.seed(100)
      if (tune) {
        mtries <- sort.int(sample(ncol(trainset)-1, 5))
        ntrees <- c(201, 501, 1501, 2501, 3501)
      } else {
        mtries <- c(rf.mtry)
        ntrees <- c(rf.ntree)
      }
      
      for(ntree in ntrees) {
        for(mtry in mtries) {
          model.rf <- randomForest(x = trainset[, colnames(trainset) != "diagnosis"],
                                   y = trainset$diagnosis, 
                                   ntree = ntree, 
                                   mtry = mtry
          )
          pred.model.rf <- predict(model.rf, newdata = testset)
          #if (!tune) {
          roc.model.rf <- multiclass.roc(response = testset$diagnosis, 
                                          predictor = predict(model.rf, 
                                                              newdata = testset[, colnames(testset) != "diagnosis"], 
                                                              type = "prob"),
                                          percent = TRUE)
          #}
          confMatrix.model.rf <- confusionMatrix(pred.model.rf, testset$diagnosis)
          
          #if (tune && confMatrix.model.rf$overall['AccuracyPValue'] < 0.05) {
          if (tune) {
            print(paste(ntree," and ", mtry))
            print(confMatrix.model.rf$overall['Accuracy'])
            print(confMatrix.model.rf$overall['AccuracyPValue'])
            print(roc.model.rf$auc)
          }
          #var.model.rf <- varImp(model.rf, useModel = TRUE, nonpara = TRUE, scale = TRUE)
          var.model.rf <- varImp(model.rf, scale = FALSE)
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
      #roc.model.auto <- multiclass.roc(response = testset$diagnosis, 
      #                                predictor = predict(model.auto, 
      #                                                    newdata = testset[, colnames(testset) != "diagnosis"], 
      #                                                    type = "prob"),
      #                                percent = TRUE)
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
      #var.model.auto <- data.frame("feature" = colnames(var.model.auto$data), "Overall" = var.model.auto$imp)
      #var.model.auto <- head(var.model.auto, -1)
      #var.model.auto <- var.model.auto[order(-var.model.auto$Overall),]
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
    data.multinomial <- createDataset(dataSource = data.pp, feature = features.gf, probeset = huex.probes, filename = "Gene Filtering")
    # Part 4.1: Compute the alpha using Bonferroni
    alpha <- get_alpha(alpha, nrow(data.multinomial) - 1)
    sets <- buildTrainTest(data.multinomial)
    trainset.multinomial <- sets$trainset
    testset.multinomial <- sets$testset
    sets <- buildTrainTest(data.multinomial)
    validationset.multinomial <- sets$validationset
    remove(sets)
      # Part 4.3.1.1: Perform tuning for SVM and Random Forest
      perform_learning("SVM", trainset.multinomial, validationset.multinomial, tune = TRUE)
      perform_learning("RF", trainset.multinomial, validationset.multinomial, tune = TRUE)
      # Part 4.3.1.2: Perform analysis
        # Part 4.3.1.2.1: Do ensemble of all learning methods
        learn.gf.auto <- perform_learning("auto", trainset.multinomial, testset.multinomial)
        # Part 4.3.1.2.2: Naive Bayes
        learn.gf.nb <- perform_learning("NB", trainset.multinomial, testset.multinomial)
        # Part 4.3.1.2.3: KNN
        learn.gf.knn <- perform_learning("KNN", trainset.multinomial, testset.multinomial)
        # Part 4.3.1.2.4: SVM
        learn.gf.svm <- perform_learning("SVM", trainset.multinomial, testset.multinomial, svm.kernel = 'laplacedot', svm.cost = 100)
        # Part 4.3.1.2.5: Logistic Regression
        learn.gf.log <- perform_learning("LOG", trainset.multinomial, testset.multinomial)
        # Part 4.3.1.2.6: Discriminant Analysis
        learn.gf.da <- perform_learning("DA", trainset.multinomial, testset.multinomial)
        # Part 4.3.1.2.7: Decision Tree
        learn.gf.dt <- perform_learning("DT", trainset.multinomial, testset.multinomial, export.filename = paste(exportdir, exportsubdir, "Decision Tree (Gene Filter).pdf", sep = "/"))
        # Part 4.3.1.2.8: Random Forest
        learn.gf.rf <- perform_learning("RF", trainset.multinomial, testset.multinomial, rf.ntree = 201, rf.mtry = 10)
    # Part 4.3.2: PCA Dataset
    data.multinomial <- createDataset(dataSource = data.pp, feature = features.pca, probeset = huex.probes, filename = "PCA")
    # Part 4.1: Compute the alpha using Bonferroni
    alpha <- get_alpha(alpha, nrow(data.multinomial) - 1)
    sets <- buildTrainTest(data.multinomial)
    trainset.multinomial <- sets$trainset
    testset.multinomial <- sets$testset
    sets <- buildTrainTest(data.multinomial)
    validationset.multinomial <- sets$validationset
    remove(sets)
      # Part 4.3.2.1: Perform tuning for SVM and Random Forest
      perform_learning("SVM", trainset.multinomial, validationset.multinomial, tune = TRUE)
      perform_learning("RF", trainset.multinomial, validationset.multinomial, tune = TRUE)
      # Part 4.3.2.2: Perform analysis
        # Part 4.3.2.2.1: Do ensemble of all learning methods
        learn.pca.auto <- perform_learning("auto", trainset.multinomial, testset.multinomial)
        # Part 4.3.2.2.2: Naive Bayes
        learn.pca.nb <- perform_learning("NB", trainset.multinomial, testset.multinomial)
        # Part 4.3.2.2.3: KNN
        learn.pca.knn <- perform_learning("KNN", trainset.multinomial, testset.multinomial)
        # Part 4.3.2.2.4: SVM (For PCA, SVM tuning had no significant features)
        learn.pca.svm <- perform_learning("SVM", trainset.multinomial, testset.multinomial, svm.kernel = 'rbfdot', svm.cost = 1)
        # Part 4.3.2.2.5: Logistic Regression
        learn.pca.log <- perform_learning("LOG", trainset.multinomial, testset.multinomial)
        # Part 4.3.2.2.6: Discriminant Analysis
        learn.pca.da <- perform_learning("DA", trainset.multinomial, testset.multinomial)
        # Part 4.3.2.2.7: Decision Tree
        learn.pca.dt <- perform_learning("DT", trainset.multinomial, testset.multinomial, export.filename = paste(exportdir, exportsubdir, "Decision Tree (PCA).pdf", sep = "/"))
        # Part 4.3.2.2.8: Random Forest  (For PCA, SVM tuning had no significant features)
        learn.pca.rf <- perform_learning("RF", trainset.multinomial, testset.multinomial, rf.ntree = 501, rf.mtry = 6)
    # Part 4.3.3: Combined Dataset
    data.multinomial <- createDataset(dataSource = data.pp, feature = features, probeset = huex.probes, filename = "PCA + GF")
    # Part 4.1: Compute the alpha using Bonferroni
    alpha <- get_alpha(alpha, nrow(data.multinomial) - 1)
    sets <- buildTrainTest(data.multinomial)
    trainset.multinomial <- sets$trainset
    testset.multinomial <- sets$testset
    sets <- buildTrainTest(data.multinomial)
    validationset.multinomial <- sets$validationset
    remove(sets)
      # Part 4.3.3.1: Perform tuning for SVM and Random Forest
      perform_learning("SVM", trainset.multinomial, validationset.multinomial, tune = TRUE)
      perform_learning("RF", trainset.multinomial, validationset.multinomial, tune = TRUE)
      # Part 4.3.3.2: Perform analysis
        # Part 4.3.3.2.1: Do ensemble of all learning methods
        learn.features.auto <- perform_learning("auto", trainset.multinomial, testset.multinomial)
        # Part 4.3.3.2.2: Naive Bayes
        learn.features.nb <- perform_learning("NB", trainset.multinomial, testset.multinomial)
        # Part 4.3.3.2.3: KNN
        learn.features.knn <- perform_learning("KNN", trainset.multinomial, testset.multinomial)
        # Part 4.3.3.2.4: SVM
        learn.features.svm <- perform_learning("SVM", trainset.multinomial, testset.multinomial, svm.kernel = 'laplacedot', svm.cost = 0.1)
        # Part 4.3.3.2.5: Logistic Regression
        learn.features.log <- perform_learning("LOG", trainset.multinomial, testset.multinomial)
        # Part 4.3.3.2.6: Discriminant Analysis
        learn.features.da <- perform_learning("DA", trainset.multinomial, testset.multinomial)
        # Part 4.3.3.2.7: Decision Tree
        learn.features.dt <- perform_learning("DT", trainset.multinomial, testset.multinomial, export.filename = paste(exportdir, exportsubdir, "Decision Tree (GF + PCA).pdf", sep = "/"))
        # Part 4.3.3.2.8: Random Forest  
        learn.features.rf <- perform_learning("RF", trainset.multinomial, testset.multinomial, rf.ntree = 1501, rf.mtry = 10)
        
        learn.features.soft <- perform_learning("SOFT", trainset.multinomial, testset.multinomial)
        