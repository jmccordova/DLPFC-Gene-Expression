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
  
  # Step 5.6: Listed all feature ranking from Importance, Auto mode, and Softvoting mode
    # Step 5.6.1: Create matrix containing all features
    matrix <- c(features,
                c(sub("feature_", "", learn.features.auto$var$feature), rep("", length(features) - nrow(learn.features.auto$var))),
                c(sub("feature_", "", learn.features.soft$var$feature), rep("", length(features) - nrow(learn.features.soft$var)))
              )
    matrix <- as.data.frame(matrix(matrix, ncol = 3))
    colnames(matrix) <- c("Features", "Auto", "Soft")
    # Step 5.6.2: Write into a CSV files
    write.csv(matrix, paste(exportdir, exportsubdir, "Rminer features.csv", sep = "/"), row.names = TRUE)
  
  # Step 5.7: Display common features between multinomial and binomial
  ggvenn(list('Multinomial' = c("3157311", "3744589", "3352040", "2844453", "3632862", "3134828", "3500787", "2773545", "2406361", "3890913", "2968401", "3005266", "3354095", "3475038", "3337168", "3346453", "2568968", "2674432", "2601230", "3405587"), 
              'Binomial' = c("3354095", "3302740", "2831519", "3428268", "2804221", "3500787", "3843934", "3155489", "2536508", "3749734", "3699810", "2672467", "2626097", "3157311", "3980614", "3017068", "3048778", "3096512", "2674432", "3556202")
              ),
  show_elements = TRUE,
  label_sep = "\n",
  text_size = 2.5
  )
  
  # Step 5.8: Display common features across the disorders
  ggvenn(list('BD' = c("3157311", "3352040", "3556094", "3134828", "3744589", "2385696", "3809690", "3337168", "2844453", "3005266", "2968401", "3405587", "2568968", "2674432", "2406361", "3500787", "3890913", "3780271", "3354095", "3818376"), 
              'MDD' = c("3157311", "3556094", "3134828", "3352040", "3744589", "3337168", "3890913", "3475038", "3780271", "2844453", "2385696", "2968401", "3632862", "2406361", "3809690", "2531129", "3585272", "3346453", "2601230", "3005266"), 
              'SCZ' = c("2844453", "3808600", "2674432", "3890913", "2359439", "3780271", "3475038", "2674317", "2568968", "3407583", "2531129", "3556274", "3352040", "2773545", "3104698", "3556202", "3744589", "3708874", "2601230", "3061438"
              )
        ),
        show_elements = TRUE,
        label_sep = "\n",
        text_size = 2.5
  )
  
  