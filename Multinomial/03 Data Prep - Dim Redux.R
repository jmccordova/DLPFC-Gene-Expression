# Part 3: Data Preparation
  # Part 3.1: Gene Filtering
    # Part 3.1.1: Coefficient of Variation (cv)
    # cv = 1 - standard deviation equals the mean, so that the experimental effect is small relative to the precision of measurement. If,
    # cv < 0.2 - the mean is ¯ve times larger than the standard deviation, so that both the experimental effect and the measurement precision are large.
    cvval <- apply(exprs(data.pp), 1, function(x){sd(x)/abs(mean(x))})
    features.gf.1 <- rownames(data[cvval < 0.2, ])
    remove(cvval)
    
    # Part 3.1.2: Normality and T-test
    f1 <- function(x) (shapiro.test(x)$p.value > alpha)
    f2 <- function(x) (t.test(x ~ diagnosis.binom.fact)$p.value < alpha)    # Ha - True differences between group means is zero
    sel1 <- genefilter(exprs(e)[, which(diagnosis.fact != 9)], filterfun(f1))   # Select diagnosed patients
    sel2 <- genefilter(exprs(e)[, which(diagnosis.fact == 9)], filterfun(f1))   # Select control patients
    sel3 <- genefilter(exprs(e), filterfun(f2))
    selected <- sel1 & sel2 & sel3
    features.gf.2 <- rownames(data[selected, ])
    remove(f1, f2, sel1, sel2, sel3, selected)
    
    # Part 3.1.3: ANOVA
    # Get all with diagnosis
    diagnosis.fact.named <- diagnosis
    diagnosis.fact.named[which(diagnosis.fact.named == 1)] <- "BPD"
    diagnosis.fact.named[which(diagnosis.fact.named == 2)] <- "MDD"
    diagnosis.fact.named[which(diagnosis.fact.named == 3)] <- "SCZ"
    diagnosis.fact.named[which(diagnosis.fact.named == 9)] <- "CTL"
    data.pp.diagnosed <- data.pp[, which(diagnosis.fact.named != "CTL")]
    
    panova <- apply(exprs(data.pp), 1, function(x) anova(lm(x ~ diagnosis.fact))$Pr[1])
    features.gf.3 <- featureNames(data.pp)[panova < alpha]
    remove(panova, data.pp.diagnosed)
  
    # Part 3.1.4: Combine similar features
    features.gf <- Reduce(intersect, List(features.gf.1, features.gf.2, features.gf.3))
    # Part 3.1.5: Check the collinearity of each predictor
      # Part 3.1.5.1: Look for perfect collinearity
      data.gf <- as.matrix(exprs(data.pp)[features.gf, ])
      show_perfect_collinearity(data.gf)
    # Part 3.1.6: Since there's no perfect collinearity, proceed
    data.gf <- as.matrix(exprs(data.pp)[features.gf, ])
    data.gf <- rbind(data.gf, diagnosis.binom)
    data.gf <- as.data.frame(t(data.gf))
    model.gf <- lm(diagnosis.binom ~ ., data = data.gf)
      # Part 3.1.5.1: Get only the factors with no to moderate collinearity (VIF <= 5)
      features.gf <- names(vif(model.gf)[vif(model.gf) <= 5])
      features.gf <- gsub("`", "", features.gf, fixed = T)
      remove(data.gf, model.gf)
    # Part 3.1.6: Put to HTML the names of the features
    atab <- aafTableAnn(features.gf, "pd.huex.1.0.st.v2", aaf.handler()[c(1:3,8:9,11:13)])
    saveHTML(atab, file=paste(datadir, "../Export/Gene Filtering Probe Names.html", sep = ""))
    # Part 3.1.7: Create a correlation matrix across each features
    features.gf.corr <- rcorr(as.matrix(t(exprs(data.pp)[features.gf, ])))
    corrplot(features.gf.corr$r)
    
    # Part 3.1.6: create Venn diagram and display all sets 
    ggvenn(list('Coefficient of Variance' = features.gf.1, 
                'T-Test' = features.gf.2, 
                'ANOVA' = features.gf.3),
           digits = 2
    )
    
  # Part 3.2: Principal Component Analysis
    # Step 3.2.1: Check for null values (If it returned more than 0, there is a null value)
    print(colSums(is.na(exprs(data.pp)))[colSums(is.na(exprs(data.pp))) != 0])
    # Step 3.2.3: Perform PCA
    model.pca <- pca(exprs(data.pp), removeVar = 0.1)
    # Step 3.2.4: Determine optimum number of PCs 
      # Step 3.2.4.1: Horn's method
      pca.loadings.horn <- parallelPCA(exprs(data.pp))
      # Step 3.2.4.2: Find the elbow point
      pca.loadings.elbow <- findElbowPoint(model.pca$variance)
      # Step 3.2.4.3: Get the PCs with that can explain at least 50% of the variability
      print(paste('Number of PCs to consider:'))
      print(paste('Horn:', pca.loadings.horn$n, '%Varation:', sum(model.pca$variance[1:pca.loadings.horn$n])))
      print(paste('Elbow:', pca.loadings.elbow, '%Variation:', sum(model.pca$variance[1:pca.loadings.elbow])))
      pca.loadings.min <- min(pca.loadings.horn$n, pca.loadings.elbow)
    # Step 3.2.5: Get the loadings
    model.pca.important <- as.data.frame(model.pca$loadings[, 1:pca.loadings.min])
    # Step 3.2.6: Export
    write.csv(model.pca.important, paste(datadir, "../Export/model.pca.important.csv", sep = ""), row.names = TRUE)
    # Step 3.2.7: Visualize PCA
      # Step 3.2.7.1: Scree plot to see where to cut the data
      screeplot(model.pca,
                components = getComponents(model.pca, 1:20),
                vline = c(pca.loadings.horn$n, pca.loadings.elbow)) +
      geom_label(aes(x = pca.loadings.horn$n + 1, y = 50,
                label = 'Horn\'s', vjust = -1, size = 8)) +
      geom_label(aes(x = pca.loadings.elbow + 1, y = 50,
                label = 'Elbow method', vjust = -1, size = 8))
      # Step 3.2.7.2: Pairs plot to compare one PC against another across all 5 PCs.
      pairsplot(model.pca)
      # Step 3.2.7.3: Biplot
      biplot(model.pca, showLoadings = TRUE, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
      # Step 3.2.7.4: Top 1%
      plotloadings(model.pca,
                   rangeRetain = 0.01,
                   labSize = 4.0,
                   title = 'Loadings plot',
                   subtitle = 'PC1, PC2, PC3, PC4, PC5',
                   caption = 'Top 1% variables',
                   shape = 24,
                   col = c('limegreen', 'black', 'red3'),
                   drawConnectors = TRUE)
      # Step 3.2.7.5: Top 10%
      features.pca.1 <- plotloadings(model.pca,
                   components = getComponents(model.pca, 1:pca.loadings.min),
                   rangeRetain = 0.1,
                   labSize = 4.0,
                   absolute = FALSE,
                   title = 'Loadings plot',
                   subtitle = 'Misc PCs',
                   caption = 'Top 10% variables',
                   shape = 23, shapeSizeRange = c(1, 16),
                   col = c('white', 'pink'),
                   drawConnectors = FALSE)
      # Step 3.2.7.6: Get the features of Top 10%
      features.pca.1 <- features.pca.1$data$var
    # Step 3.2.8: Transform to prcomp
    #model.pca.prcomp <- BiocSingular::runPCA(t(exprs(e)), rank = pca.loadings.min, scale = TRUE)
    model.pca.prcomp <- list(sdev = model.pca$sdev,
                          rotation = data.matrix(model.pca$loadings),
                          x = data.matrix(model.pca$rotated),
                          center = TRUE, scale = TRUE
                        )
    class(model.pca.prcomp) <- 'prcomp'
    # Step 3.2.9: Visualization of PCs
      # Step 3.2.9.1: Graph of individuals. Individuals with a similar profile are grouped together.
      fviz_pca_ind(model.pca.prcomp,
                    label="ind", 
                    habillage=diagnosis.fact.named,
                    addEllipses=TRUE, 
                    ellipse.level=0.95,
                    repel = TRUE
                  )
      # Step 3.2.9.2: Graph of variables. Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.
      features.pca.2 <- fviz_pca_var(model.pca.prcomp, 
                   alpha.var="contrib",
                   col.var = "red3",
                   repel = TRUE,
                   select.var = list(contrib = length(features.pca.1)),
                   ) + theme_minimal()
      features.pca.2 <- as.character(features.pca.2$data[, 'name'])
      # Step 3.2.9.3: Biplot of individuals and variables
      fviz_pca_biplot(model.pca.prcomp, repel = TRUE,
                      col.var = "#2E9FDF", # Variables color
                      col.ind = "#696969",  # Individuals color,
                      select.var = list(contrib = round(nrow(data.pp) * 0.01))
      )
    # Step 3.2.10: For each of the principal component, get the variable with highest magnitude of eigenvalues
    features.pca <- Reduce(intersect, List(features.pca.1, features.pca.2))
    # Step 3.2.11: Check for multicollinearity
      data.pca <- as.matrix(exprs(data.pp)[features.pca, ])
      show_perfect_collinearity(data.pca)
      # Step 3.2.11.1: From a correlation testing, features "2908474" "2908505" have perfect collinearity; remove them
      features.pca <- features.pca[features.pca != c("2908474", "2908505")]
      # Step 3.2.11.2: Proceed with getting VIF
      data.pca <- as.matrix(exprs(data.pp)[features.pca, ])
      data.pca <- rbind(data.pca, diagnosis.binom)
      data.pca <- as.data.frame(t(data.pca))
      model.pca <- lm(diagnosis.binom ~ ., data = data.pca)
      # Part 3.2.11.3: Get only the factors with no to moderate collinearity (VIF <= 5)
      features.pca <- names(vif(model.pca)[vif(model.pca) <= 5])
      features.pca <- gsub("`", "", features.pca, fixed = T)
      remove(data.pca, model.pca)
    # Part 3.2.13: Put to HTML the names of the features
    atab <- aafTableAnn(features.pca, "pd.huex.1.0.st.v2", aaf.handler()[c(1:3,8:9,11:13)])
    saveHTML(atab, file=paste(datadir, "../Export/PCA Probe Names.html", sep = ""))
    # Part 3.2.14: Create a correlation matrix across each features
    features.pca.corr <- rcorr(as.matrix(t(exprs(data.pp)[features.pca, ])))
    corrplot(features.pca.corr$r)
      
    # Step 3.2.15: create Venn diagram and display all sets 
    ggvenn(list('PCATools' = features.pca.1, 
                'Factoextra' = features.pca.2
                ),
           digits = 2
    )
    
  # Step 3.3: Combine the features from Gene Filtering and PCA
  features <- unique(append(features.gf, features.pca))
    # Step 3.3.1: Look for perfect collinearity
    data.features <- as.matrix(exprs(data.pp)[features, ])
    show_perfect_collinearity(data.features)
    # Step 3.3.2: Since there is no perfect collinearity, proceed
    data.features <- as.matrix(exprs(data.pp)[features, ])
    data.features <- rbind(data.features, diagnosis.binom)
    data.features <- as.data.frame(t(data.features))
    model.features <- lm(diagnosis.binom ~ ., data = data.features)
    # Part 3.3.3: Get only the factors with no to moderate collinearity (VIF <= 5)
    features <- names(vif(model.features)[vif(model.features) <= 5])
    features <- gsub("`", "", features, fixed = T)
    remove(data.features, model.features)
    # Part 3.3.4: Put to HTML the names of the features
    atab <- aafTableAnn(features.pca, "pd.huex.1.0.st.v2", aaf.handler()[c(1:3,8:9,11:13)])
    saveHTML(atab, file=paste(datadir, "../Export/GF + PCA Probe Names.html", sep = ""))
    # Part 3.2.14: Create a correlation matrix across each features
    features.corr <- rcorr(as.matrix(t(exprs(data.pp)[features, ])))
    corrplot(features.corr$r)
    
    # Step 3.2.15: create Venn diagram and display all sets 
    ggvenn(list('Gene Filter' = features.gf, 
                'PCA' = features.pca
    ),
    digits = 2
    )
    
  # Part 3.4: Get probe profiles and annotations
    # This part was done using Ensemble Biomart with following specifications:
    # Dataset
    #   Human genes (GRCh38.p14)
    # Filters
    #   AFFY HuEx 1 0 st v2 probe ID(s) [e.g. 4037576]: [ID-list specified]
    #Attributes
    #   Gene stable ID
    #   Gene stable ID version
    #   Transcript stable ID
    #   Transcript stable ID version
    #   Exon stable ID
    #   Gene description
    #   GO term accession

  # Part 3.4: Splitting dataset
  options(scipen=999)  # prevents printing scientific notations.
  set.seed(100)
    # Part 3.4.1: Get only the chosen features
    data.multinomial <- as.matrix(exprs(data.pp)[features, ])
    # Part 3.4.2: Keep only selected probes
    huex.probes <- huex.probes[which(huex.probes$probeset_id %in% features), ]
      # Part 3.4.2.1: Show probes not in the probeset annotation
      features.missed <- features[which(features %ni% huex.probes$probeset_id)]
      print(features.missed)
    # Part 3.4.3: Insert the diagnosis factor in the dataframe
    data.multinomial <- rbind(data.multinomial, diagnosis)
    data.multinomial <- as.data.frame(t(data.multinomial))
    data.multinomial$diagnosis <- factor(data.multinomial$diagnosis, ordered = FALSE)
    # Part 3.4.4: Choose the index on which to choose as training
    index.multinomial <- createDataPartition(data.multinomial$diagnosis, p = 0.75, list = F)
    # Part 3.4.5: To make the training better, perform upsampling
    trainset.multinomial <- data.multinomial[index.multinomial, ]
    trainset.multinomial <- upSample(trainset.multinomial[, names(trainset.multinomial) %ni% c("diagnosis")], trainset.multinomial$diagnosis, yname = "diagnosis")
    # Part 3.4.6: Correct the factors in testing
    testset.multinomial <- data.multinomial[-index.multinomial, ]
    testset.multinomial$diagnosis <- factor(testset.multinomial$diagnosis, ordered = FALSE)
    
    # Step 3.4.7: Export
    write.csv(data.multinomial, paste(datadir, "../Export/Chosen Dataset.csv", sep = ""), row.names = TRUE)
    write.csv(huex.probes, paste(datadir, "../Export/Chosen Probes.csv", sep = ""), row.names = TRUE)
    
  remove(e, data)

  results <- c()
  for ( probe in features.missed ){
    q <- paste(sep="","https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&input=affyid&inputValues=",probe,"&outputs=genesymbol&taxonId=9606&format=row")
    results <- rbind(results,fromJSON(txt=q))
  }
  