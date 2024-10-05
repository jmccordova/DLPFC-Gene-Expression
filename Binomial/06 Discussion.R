# Part 6 - Discussion
exportsubdir <- "Step 6 - Discussion"
dir.create(paste(exportdir, exportsubdir, sep = "/"), recursive=TRUE)
# Part 6.1. Make the diagnosis column a named one
# Part 6.1.1: Create a function that transforms numeric factors to strings
labelDiagnosis <- function(data) {
  a <- data$diagnosis
  a <- sapply(a, as.character)
  a[which(a == "1")] <- "NCTL"
  a[which(a == "9")] <- "CTL"
  return(factor(a))
}
# Part 6.1.2: Replace the diagnosis column of trainset and testset
trainset.binomial.important$diagnosis <- labelDiagnosis(trainset.binomial.important)
testset.binomial.important$diagnosis <- labelDiagnosis(testset.binomial.important)

# Part 6.2: Using LIME
# Part 6.2.1: Create folder containing the LIME
dir.create(paste(exportdir, exportsubdir, "LIME", sep = "/"), recursive=TRUE)
# Part 6.2.2: Since LIME only accept train object, I need to rebuild SVM with radial using the same values
explainer.binomial.model <- caret::train(diagnosis ~ ., 
                                            data = trainset.binomial.important,
                                            method = "svmRadial",
                                            preProc = c("center", "scale", "nzv"), 
                                            family = "binomial",
                                         type = "prob",
                                            trControl = trainControl(
                                              method = "cv", 
                                              savePred = TRUE, 
                                              classProb = TRUE,
                                              verboseIter = FALSE
                                            )
)

explainer.binomial.model <- kernlab::ksvm(
  diagnosis ~ .,
  data = trainset.binomial.important,
  C = 0.0170484,
  prob.model = TRUE
)

# Part 6.2.3: Create an explainer object
explainer.binomial <- lime(trainset.binomial.important, explainer.binomial.model)

# Part 6.2.4: Create a list of number of features to play around
set.seed(12465)
explainer.binomial.n_features <- sort.int(sample(ncol(trainset.binomial.important) - 1, 5))

model_type.ksvm <- function(x) {
  return("classification")
}

predict_model.ksvm <- function(x, newdata, type = "prob") {
  # return classification probabilities only   
  # “response”, “probabilities”, “votes”, “decision”
  res <- predict(x, newdata, type = "probabilities")
  print(head(res))
  # Error in res[, c("model_type", "case", "label", "label_prob", "model_r2",  : 
  # Add those columns maybe?
  return(res)
}

# Part 6.2.5: For each disorder, create an explanation
for(n_feature in explainer.binomial.n_features) {
  # Part 6.2.5.1: Create an explanation object for LIME for a specific disorder (labels) and a number of features to consider (n_features); Onlt take a subset as LIME forbids using all to avoid memory exhaustion
  explanation.binomial <- lime::explain(
    x = testset.binomial.important[15:25, ], 
    explainer = explainer.binomial, 
    n_labels = 1, 
    n_features = n_feature
  )
  
  # Part 6.2.5.2: Export to a file
  ggsave2(
    filename = paste(exportdir, exportsubdir, i, paste("LIME Binomial ", n_feature, ".png", sep = ""), sep = "/"), 
    device = "png",
    plot = lime::plot_features(explanation.binomial),
    width = 20, 
    height = 18, 
    dpi = 150, 
    units = "in"
  )
}

# Part 6.3: Using SHAP
# Part 6.3.1: Create a ranger object similar to RF since fastshap only takes ranger
ranger.binomial <- ranger::ranger(
  diagnosis ~ ., 
  data = trainset.binomial.important, 
  probability = TRUE
)

# Part 6.3.2: Define the function for prediction wrapper
pfun <- function(object, newdata) {
  diagnosis <- "9" # Replace diagnosis
  unname(predict(object, data = newdata)$predictions[, diagnosis])
}

# Part 6.3.3: Explain each features in a global scope
set.seed(101) # for reproducibility
shap <- fastshap::explain(
  ranger.binomial, 
  X = trainset.binomial.important[, colnames(trainset.binomial.important) != "diagnosis"], 
  pred_wrapper = pfun,
  nsim = 100, 
  adjust = TRUE,
  shap_only = FALSE
)

# Part 6.3.4: Transform to tibble
tibble::as_tibble(shap$shapley_values)
# Part 6.3.5: Visualize the SHAP features
shv.global <- shapviz(shap)
# Part 6.3.6: Plot the important features
sv_importance(shv.global)  
