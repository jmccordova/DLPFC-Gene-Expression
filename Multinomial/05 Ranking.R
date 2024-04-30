#ranking <- list(
#  naive_bayes = list(
#    gf = rownames(learn.gf.nb$var$importance),
#    pca = rownames(learn.pca.nb$var$importance),
#    features = rownames(learn.features.nb$var$importance)
#  ),
#  knn = list(
#    gf = rownames(learn.gf.knn$var$importance),
#    pca = rownames(learn.pca.knn$var$importance),
#    features = rownames(learn.features.knn$var$importance)
#  ),
#  svm = list(
#    gf = rownames(learn.gf.svm$var$importance),
#    pca = c(),
#    features = rownames(learn.features.svm$var$importance)
#  ),
#  logistic_regression = list(
#    gf = rownames(learn.gf.log$var$importance),
#    pca = rownames(learn.pca.log$var$importance),
#    features = rownames(learn.features.log$var$importance)
#  ),
#  discriminant_analysis = list(
#    gf = rownames(learn.gf.da$var$importance),
#    pca = rownames(learn.pca.da$var$importance),
#    features = rownames(learn.features.da$var$importance)
#  ),
#  decision_tree = list(
#    gf = rownames(learn.gf.dt$var$importance),
#    pca = rownames(learn.pca.dt$var$importance),
#    features = rownames(learn.features.dt$var$importance)
#  ),
#  random_forest = list(
#    gf = rownames(learn.gf.rf$var$importance),
#    pca = c(),
#    features = c()
#  )
#)

ranking <- list(
  naive_bayes.gf = rownames(learn.gf.nb$var$importance)[1:20],
  naive_bayes.pca = rownames(learn.pca.nb$var$importance)[1:20],
  naive_bayes.features = rownames(learn.features.nb$var$importance)[1:20],
  knn.gf = rownames(learn.gf.knn$var$importance)[1:20],
  knn.pca = rownames(learn.pca.knn$var$importance)[1:20],
  knn.features = rownames(learn.features.knn$var$importance)[1:20],
  svm.gf = rownames(learn.gf.svm$var$importance)[1:20],
  svm.pca = rep('NA', 20),
  svm.features = rownames(learn.features.svm$var$importance)[1:20],
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
  random_forest.pca = rep('NA', 20),
  random_forest.features = rep('NA', 20)
)

write.csv(ranking, paste(datadir, "../Export/Ranking.csv", sep = ""), row.names = TRUE)
