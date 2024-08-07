# Part 5 - Ranking
exportsubdir <- "Step 5 - Ranking"
dir.create(paste(exportdir, exportsubdir, sep = "/"), recursive=TRUE)
  # Step 5.1. Features in GF
    # Step 5.1.1. Build list
    ranking.gf <- list(
      auto.gf = learn.gf.auto$var$feature,
      naive_bayes.gf = rownames(learn.gf.nb$var$importance),
      knn.gf = rownames(learn.gf.knn$var$importance),
      svm.gf = learn.gf.svm$var$feature,
      logistic_regression.gf = gsub("`", "", rownames(learn.gf.log$var), fixed = T),
      discriminant_analysis.gf = rownames(learn.gf.da$var$importance),
      decision_tree.gf = rownames(learn.gf.dt$var),
      random_forest.gf = rownames(learn.gf.rf$var)
    )
    # Step 5.1.2: Export
    write.csv(ranking.gf, paste(exportdir, exportsubdir, "Ranking (GF).csv", sep = "/"), row.names = TRUE)
    
    par(mar=c(5, 5, 1, 1)); barplot(learn.gf.auto$var$Overall, horiz = TRUE, names.arg = learn.gf.auto$var$feature, las = 2)
  # Step 5.2. Features in PCA
    # Part 5.2,1. Build list
    ranking.pca <- list(
      naive_bayes.pca = rownames(learn.pca.nb$var$importance)[1:20],
      knn.pca = rownames(learn.pca.knn$var$importance)[1:20],
      svm.pca =learn.pca.svm$var$feature[1:20],
      logistic_regression.pca = gsub("`", "", rownames(learn.pca.log$var), fixed = T)[1:20],
      discriminant_analysis.pca = rownames(learn.pca.da$var$importance)[1:20],
      decision_tree.pca = rownames(learn.pca.dt$var)[1:20],
      random_forest.pca = rownames(learn.pca.rf$var)[1:20],
    )
    # Step 5.2.2: Export
    write.csv(ranking.pca, paste(exportdir, exportsubdir, "Ranking.csv", sep = "/"), row.names = TRUE)
  # Step 5.3. Features combined
    # Part 5.3.1. Build list
    ranking.features <- list(
      auto.features = learn.gf.auto$var$feature,
      naive_bayes.features = rownames(learn.gf.nb$var$importance),
      knn.features = rownames(learn.gf.knn$var$importance),
      svm.features = learn.gf.svm$var$feature,
      logistic_regression.features = gsub("`", "", rownames(learn.gf.log$var), fixed = T),
      discriminant_analysis.features = rownames(learn.gf.da$var$importance),
      decision_tree.features = rownames(learn.gf.dt$var),
      random_forest.features = rownames(learn.gf.rf$var)
    )
    # Step 5.3.2: Export
    write.csv(ranking.features, paste(exportdir, exportsubdir, "Ranking.csv", sep = "/"), row.names = TRUE)
  # Step 5.4. All but only top 20 features
    # Part 5.4.1. Build list
    ranking <- list(
      naive_bayes.gf = rownames(learn.gf.nb$var$importance)[1:20],
      naive_bayes.pca = rownames(learn.pca.nb$var$importance)[1:20],
      naive_bayes.features = rownames(learn.features.nb$var$importance)[1:20],
      knn.gf = rownames(learn.gf.knn$var$importance)[1:20],
      knn.pca = rownames(learn.pca.knn$var$importance)[1:20],
      knn.features = rownames(learn.features.knn$var$importance)[1:20],
      svm.gf = learn.gf.svm$var$feature[1:20],
      svm.pca =learn.pca.svm$var$feature[1:20],
      svm.features = learn.features.svm$var$feature[1:20],
      logistic_regression.gf = gsub("`", "", rownames(learn.gf.log$var), fixed = T)[1:20],
      logistic_regression.pca = gsub("`", "", rownames(learn.pca.log$var), fixed = T)[1:20],
      logistic_regression.features = gsub("`", "", rownames(learn.pca.log$var), fixed = T)[1:20],
      discriminant_analysis.gf = rownames(learn.gf.da$var$importance)[1:20],
      discriminant_analysis.pca = rownames(learn.pca.da$var$importance)[1:20],
      discriminant_analysis.features = rownames(learn.features.da$var$importance)[1:20],
      decision_tree.gf = rownames(learn.gf.dt$var)[1:20],
      decision_tree.pca = rownames(learn.pca.dt$var)[1:20],
      decision_tree.features = rownames(learn.features.dt$var)[1:20],
      random_forest.gf = rownames(learn.gf.rf$var)[1:20],
      random_forest.pca = rownames(learn.pca.rf$var)[1:20],
      random_forest.features = rownames(learn.features.rf$var)[1:20]
    )
    # Step 5.4.2: Export
    write.csv(ranking, paste(exportdir, exportsubdir, "Ranking.csv", sep = "/"), row.names = TRUE)
  
  # Step 5.3: Get data from annotation
  write.csv(huex.probes[, which(huex.probes$probeset_id %in% features)], paste(exportdir, exportsubdir, "Annotated Features.csv", sep = "/"), row.names = TRUE)
  
  # Step 5.4: Use logistic regression values for interpretation
    # Step 5.4.1: Rename the diagnosis
    levels(trainset.multinomial$diagnosis)[match("1",levels(trainset.multinomial$diagnosis))] <- "BPD"
    levels(trainset.multinomial$diagnosis)[match("2",levels(trainset.multinomial$diagnosis))] <- "MDD"
    levels(trainset.multinomial$diagnosis)[match("3",levels(trainset.multinomial$diagnosis))] <- "SCZ"
    levels(trainset.multinomial$diagnosis)[match("9",levels(trainset.multinomial$diagnosis))] <- "CTL"
    # Step 5.4.2: Get the model
    fit <- vglm(diagnosis ~ ., family = multinomial(refLevel="CTL"), data = trainset.multinomial) # vglm = vector GLM
    # Step 5.4.3: Get the coefficients
    coef(fit, matrix = TRUE)
    # Step 5.4.4: Get logit
    exp(coef(fit))
    # Step 5.4.5: Export the logit
    write.csv(exp(coef(fit)), paste(exportdir, exportsubdir, "LR Log-likelihood.csv", sep = "/"), row.names = TRUE)
    
  # Step 5.5: Get the coefficients in discriminant analysis
  write.csv(coef(learn.features.da$model), paste(exportdir, exportsubdir, "DA (GF + PCA).csv", sep = "/"), row.names = TRUE)
  
  # Logistic Regression tuning for parsimony
  matrix <- c()
  #for(i in seq(1, length(features))) {
  for(i in seq(1, 3)) {
    print(i)
    combinations <- t(data.frame(combn(features, i)))
    for (j in seq(1, nrow(combinations))) {
      cols <- c(combinations[j, ], 'diagnosis')
      trainset <- trainset.multinomial[, cols]
      testset <- testset.multinomial[, cols]
      
      model.logit <- multinom(diagnosis ~ .,
                              data = trainset,
                              trace = FALSE)
      if (i == 1) {
        testset.newdata <- as.data.frame(testset[, combinations[j, ]])
        colnames(testset.newdata) <- combinations[j, ]
      } else {
        testset.newdata <- testset[, combinations[j, ]]
      }
      
      roc.model.logit <- multiclass.roc(response = testset$diagnosis, 
                                        predictor = predict(model.logit, 
                                                            newdata = testset.newdata, 
                                                            type = "prob"),
                                        percent = TRUE)
      
      # Wald's Z-score for the model
      #model.logit.z <- summary(model.logit)$coefficients/summary(model.logit)$standard.errors
      # 2-tailed z test
      #model.logit.p <- (1 - pnorm(abs(model.logit.z), 0, 1)) * 2
      
      # Test the goodness of fit
      #model.logit.chisq <- chisq.test(trainset$diagnosis, predict(model.logit))
      
      # Compute R squared to get model fit
      model.logit.r2 <- PseudoR2(model.logit, which = c("CoxSnell","Nagelkerke","McFadden"))
      
      # Log Likelihood Ratio
      # lmtest::lrtest(model.logit, "<spacific feature>')
      
      matrix <- rbind(matrix, 
                      c(
                        roc.model.logit$auc[1], 
                        model.logit.r2,
                        paste(cols, sep = '', collapse = ' ')
                      )
                )
    }
  }
  matrix <- as.data.frame(matrix)
  colnames(matrix) <- c("AUC", "CoxSnell","Nagelkerke","McFadden", "Features")
  write.csv(matrix, paste(exportdir, exportsubdir, "LR Tuning.csv", sep = "/"), row.names = TRUE)
  remove(matrix, roc.model.logit, model.logit, model.logit.z, model.logit.p, model.logit.r2)
  
  