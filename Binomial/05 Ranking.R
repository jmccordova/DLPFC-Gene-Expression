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
      auto.features = learn.features.auto$var$feature,
      naive_bayes.features = rownames(learn.features.nb$var$importance),
      knn.features = rownames(learn.features.knn$var$importance),
      svm.features = learn.features.svm$var$feature,
      logistic_regression.features = gsub("`", "", rownames(learn.features.log$var), fixed = T),
      discriminant_analysis.features = rownames(learn.features.da$var$importance),
      decision_tree.features = rownames(learn.features.dt$var),
      random_forest.features = rownames(learn.features.rf$var)
    )
    # Step 5.3.2: Export
    write.csv(ranking.features, paste(exportdir, exportsubdir, "Ranking.csv", sep = "/"), row.names = TRUE)
  # Part 5.1. Build list
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
  par(mar=c(5, 5, 1, 1)); barplot(learn.gf.auto$var$Overall, horiz = TRUE, names.arg = learn.gf.auto$var$feature, las = 2)
  # Step 5.2: Export
  write.csv(ranking, paste(exportdir, exportsubdir, "Ranking.csv", sep = "/"), row.names = TRUE)

  # Step 5.3: Get data from annotation
  write.csv(huex.probes[, which(huex.probes$probeset_id %in% features)], paste(exportdir, exportsubdir, "Annotated Features.csv", sep = "/"), row.names = TRUE)
  