# Part 6 - Discussion
exportsubdir <- "Step 6 - Discussion"
dir.create(paste(exportdir, exportsubdir, sep = "/"), recursive=TRUE)
  # Part 6.1. Make the diagnosis column a named one
    # Part 6.1.1: Create a function that transforms numeric factors to strings
    labelDiagnosis <- function(data) {
      a <- data$diagnosis
      a <- sapply(a, as.character)
      a[which(a == "1")] <- "BPD"
      a[which(a == "2")] <- "MDD"
      a[which(a == "3")] <- "SCZ"
      a[which(a == "9")] <- "CTL"
      return(factor(a))
    }
    # Part 6.1.2: Replace the diagnosis column of trainset and testset
    #trainset.multinomial.important$diagnosis <- labelDiagnosis(trainset.multinomial.important)
    #testset.multinomial.important$diagnosis <- labelDiagnosis(testset.multinomial.important)

  # Part 6.2: Using LIME
    # Part 6.2.1: Create folder containing the LIME
    dir.create(paste(exportdir, exportsubdir, "LIME", sep = "/"), recursive=TRUE)
    # Part 6.2.2: Since LIME only accept train object, I need to rebuild RF using the same values
    explainer.multinomial.model <- caret::train(diagnosis ~ ., 
                 data = trainset.multinomial.important,
                 method = "rf",
                 trControl = trainControl(
                   method = "repeatedcv", 
                   number = 10,
                   repeats = 5, 
                   verboseIter = FALSE
                 )
    )
    # Part 6.2.3: Create an explainer object
    explainer.multinomial <- lime(trainset.multinomial.important, explainer.multinomial.model)
    # Part 6.2.4: Create a list of number of features to play around
    set.seed(1241)
    explainer.multinomial.n_features <- sort.int(sample(ncol(trainset.multinomial.important) - 1, 5))
    # Part 6.2.5: For each disorder, create an explanation
    for(i in c(1, 2, 3)) {
      for(n_feature in explainer.multinomial.n_features) {
        # Part 6.2.5.1: Create an explanation object for LIME for a specific disorder (labels) and a number of features to consider (n_features); Onlt take a subset as LIME forbids using all to avoid memory exhaustion
        explanation.multinomial <- lime::explain(
                                x = testset.multinomial.important[15:20, ], 
                                explainer = explainer.multinomial, 
                                labels = c(i),
                                #n_labels = 1, 
                                n_features = n_feature
                        )
        
        # Part 6.2.5.2: Export to a file
        ggsave2(
            filename = paste(exportdir, exportsubdir, i, paste("LIME Multinomial ", n_feature, ".png", sep = ""), sep = "/"), 
            device = "png",
            plot = lime::plot_features(explanation.multinomial),
            width = 20, 
            height = 18, 
            dpi = 150, 
            units = "in"
        )
      }
    }

# Part 6.3: Using SHAP
  # Part 6.3.1: Create a ranger object similar to RF since fastshap only takes ranger
  ranger.multinomial <- ranger::ranger(
    diagnosis ~ ., 
    data = trainset.multinomial.important, 
    probability = TRUE
  )

  # Part 6.3.2: Define the function for prediction wrapper
  pfun <- function(object, newdata) {
    diagnosis <- "1" # Replace diagnosis
    unname(predict(object, data = newdata)$predictions[, diagnosis])
  }

  # Part 6.3.3: Explain each features in a global scope
  set.seed(101) # for reproducibility
  shap <- fastshap::explain(
    ranger.multinomial, 
    X = trainset.multinomial.important[, colnames(trainset.multinomial.important) != "diagnosis"], 
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
