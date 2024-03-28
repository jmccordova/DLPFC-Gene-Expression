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
    # The below is not useful since these diseases are not progression but different types
    #design.ma <- model.matrix(~ 0 + factor(diagnosis.fact.named[diagnosis.fact.named != "CTL"]))
    #colnames(design.ma) <- c("BPD","MDD","SCZ")
    #fit <- lmFit(data.pp.diagnosed, design.ma)
    #fit <- eBayes(fit)
    #toptab <- topTable(fit, coef=2,5,adjust.method="fdr")
    #print(toptab[,1:5],digits=4)
    #cont.ma <- makeContrasts(BPD-MDD, MDD-SCZ, levels=factor(diagnosis.fact.named[diagnosis.fact.named != "CTL"]))
    # End
    #panova <- apply(exprs(data.pp.diagnosed), 1, function(x) anova(lm(x ~ factor(diagnosis.fact.named[diagnosis.fact.named != "CTL"])))$Pr[1])
    #genenames <- featureNames(data.pp.diagnosed)[panova < 0.000001]
    panova <- apply(exprs(data.pp), 1, function(x) anova(lm(x ~ diagnosis.fact))$Pr[1])
    features.gf.3 <- featureNames(data.pp)[panova < alpha]
    atab <- aafTableAnn(features.gf.3, "pd.huex.1.0.st.v2", aaf.handler()[c(1:3,8:9,11:13)])
    saveHTML(atab, file=paste(datadir, "Export/ANOVA Probe Names.html", sep = ""))
    remove(panova, atab)
    
    
  # Part 3.1: Principal Component Analysis
    # Part 3.1.5: Dimension reduction using Principal Component Analysis
      # Step 3.1.5.1: Check for null values (If it returned more than 0, there is a null value)
      print(colSums(is.na(data))[colSums(is.na(data)) != 0])
      # Step 3.1.5.2: Normalize data
      # data.normalized <- scale(data)
      data.normalized <- apply(data, 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
      # Step 3.1.5.3: Compute correlation matrix
      corr_matrix <- cor(data.normalized)
        # Step 3.1.5.3.1: Visualize the correlation matrix
        qgraph(corr_matrix, minimum = 0.25, cut = 0.4, vsize = 2, legend = TRUE, borders = FALSE)
      # Step 3.1.5.4: Computer covariance matrix
      cov_matrix <- cov(data.normalized)
        # Step 3.1.5.4.1: Check if there are NaN in covariance matrix
        print(any(is.nan(cov_matrix)))
      # Step 3.1.5.5: Apply PCA
        # Step 3.1.5.5.1: Transpose the data first so the probe IDs are the rows
        data.t <- data.frame(t(data))
        #rownames(data.t) <- transcriptIDs$transcript_id[match(rownames(data.t), transcriptIDs$id_internal_huex)]
        # Step 3.1.5.5.2: Check if the data rows and metadata columns match names
        all(rownames(data.t) %in% colnames(data.pca.metadata))
        all(rownames(data.t) == colnames(data.pca.metadata))
        # Step 3.1.5.5.3: Perform PCA
        #data.pca <- pca(data.t, metadata = data.pca.metadata, removeVar = 0.1)
        data.pca <- pca(data.t, removeVar = 0.1)
        # Step 3.1.5.5.4: Determine optimum number of PCs using Horn's method
        horn <- parallelPCA(data.t)
        # Step 3.1.5.5.5: Find the elbow point
        elbow <- findElbowPoint(data.pca$variance)
        # Get the PCs with that can explain at least 50% of the variability
        print(paste('Number of PCs to consider:'))
        print(paste('Horn:', horn$n, '%Varation:', sum(data.pca$variance[1:horn$n])))
        print(paste('Elbow:', elbow, '%Variation:', sum(data.pca$variance[1:elbow])))
        # Step 3.1.5.5.6: Put the important PCAs in a CSV
        data.important <- data.pca$loadings[, 1:horn$n]
        write.csv(data.important, paste(datadir, "data.important.csv", sep = ""), row.names = TRUE)
        # Step 3.1.5.5.7: Visualize the loadings
        plotloadings(data.pca,
                     rangeRetain = 0.01,
                     labSize = 4.0,
                     title = 'Loadings plot',
                     subtitle = 'PC1, PC2, PC3, PC4, PC5',
                     caption = 'Top 1% variables',
                     shape = 24,
                     col = c('limegreen', 'black', 'red3'),
                     drawConnectors = TRUE)
        plotloadings(data.pca,
                     components = getComponents(data.pca, 1:elbow),
                     rangeRetain = 0.1,
                     labSize = 4.0,
                     absolute = FALSE,
                     title = 'Loadings plot',
                     subtitle = 'Misc PCs',
                     caption = 'Top 10% variables',
                     shape = 23, shapeSizeRange = c(1, 16),
                     col = c('white', 'pink'),
                     drawConnectors = FALSE)
        screeplot(data.pca,
          components = getComponents(data.pca, 1:20),
          vline = c(horn$n, elbow)) +
          geom_label(aes(x = horn$n + 1, y = 50,
            label = 'Horn\'s', vjust = -1, size = 8)) +
          geom_label(aes(x = elbow + 1, y = 50,
            label = 'Elbow method', vjust = -1, size = 8)
        )

        biplot(data.pca, showLoadings = TRUE,
               labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
        pairsplot(data.pca)

# Part 3.2: Factor Analysis
nScree(cov(data))
data.corr <- foreach(i = seq_len(ncol(data)),
             .combine = rbind,
             .multicombine = TRUE,
             .inorder = FALSE,
             .packages = c('data.table', 'doParallel')) %dopar% {
                print(i)
                cor(data[,i], data, method = 'pearson')
             }
write.csv(data.corr, file = paste(datadir, '../Export/corr_matrix.csv', sep = ''), row.names = TRUE)
data.eigen <- foreach(i = seq_len(ncol(data)),
             .combine = rbind,
             .multicombine = TRUE,
             .inorder = FALSE,
             .packages = c('data.table', 'doParallel')) %dopar% {
               print(i)
               cor(data[,i], data, method = 'pearson')
             }
eigen(data.corr)



raw <- RAWPAR(data, 
            randtype = "permuted", 
            extraction = "PCA", 
            Ndatasets = round(nrow(data) * 0.25),
            percentile = 95, 
            corkind = "pearson", 
            corkindRAND = NULL, 
            Ncases = NULL, 
            verbose = TRUE
          )


data.factors.n <- raw$values
# Using kaiser method, get eigen > 1
factanal(data, factors=data.factors.n, rotation="varimax") 
factanal(data, factors=data.factors.n, rotation="oblimin")
factanal(data, factors=data.factors.n, rotation="none")









# Part 3.1.5: Splitting dataset
## For multinomial
'%ni%' <- Negate('%in%')  # define 'not in' func
options(scipen=999)  # prevents printing scientific notations.
set.seed(100)
index.multinomial <- createDataPartition(data.multinomial$diagnosis, p=0.75, list = F)
trainset.multinomial <- data.multinomial[index.multinomial, ]
testset.multinomial <- data.multinomial[-index.multinomial, ]

