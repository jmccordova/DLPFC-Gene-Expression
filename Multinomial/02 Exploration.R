# Part 2: Exploration
exportsubdir <- "Part 2 - Exploration"
  # Part 2.1: Get information about the dataset
  slotNames(rawData)
  str(rawData)
  head(rawData)
  # Part 2.2: Show how many gene expressions in the raw data
  dim(exprs(rawData))   # 6553600     169
  # Part 2.3: Platform used for gene information
  annotation(rawData)   # "pd.huex.1.0.st.v2"
  # Part 2.4: Show all the probes used for this dataset
  unique(probeNames(rawData))
  # Part 2.5:
  matplot(pm(rawData, subset = NULL, target='core'), type="l", xlab="Probe No.", ylab="PM Probe intensity")    # variability of the data within and between probes
  #hist(rawData)   # Density plots of the log of the probe values
  #par(mar=c(1,1,1,1))
  #MAplot(rawData, pairs=TRUE, plot.method= "smoothScatter")
  # image(rawData)    # Display an image of the expressions for each sample

  # Part 2.6: Turn probe-level information to gene-level data by normalization using RMA
  # RMA - Robust Multi-Array Average (http://www.sthda.com/english/wiki/affymetrix-cel-files)
  # Returns ExpressionSet
  e <- rma(rawData)
  # Part 2.7: Show details about the preprocessed values
  str(e)
  head(e)
  
  # Part 2.8: Rename the columns of exprs(e)
  colnames.e <- colnames(e)
  for (i in 1:length(colnames.e)) {
    colnames.e[i] <- strsplit(colnames.e[i], split = "_")[[1]][1]
  }
  colnames(e) <- colnames.e
  data <- as.data.frame(exprs(e))
  
  # Part 2.8: View data as box and whiskers
  boxplot(data)
  
  # Part 2.9: Check if expression set is normally distributed
  p.vals <- rbind(
              apply(data, 1, function(x) shapiro.test(as.matrix(x))$p.value), 
              apply(data, 1, function(x) kruskal.test(as.numeric(x) ~ as.factor(x))$p.value)
            )
  p.vals <- as.data.frame(p.vals)
  row.names(p.vals) <- c("Normality: Shapiro-Wilk", "Population Equality: Kruskal-Wallis")
  print(p.vals[1:10])
  
  # Part 2.10: Create factors for each diagnosis
    # Part 2.10.1: Replace string to integer classification
    # Diagnosis from each sample in GEO
    # CTL - 9, BPD - 1, MDD - 2, SCZ - 3
    diagnosis <- c(
      # 299-309
      9, 3, 2, 3, 1, 1, 2, 3, 3, 9, 3,
      # 310-320
      2, 1, 3, 9, 3, 2, 2, 1, 1, 9, 1,
      # 321-331
      1, 9, 1, 9, 1, 9, 9, 9, 3, 2, 9,
      # 332-342
      2, 9, 3, 2, 9, 2, 3, 9, 3, 2, 2,
      # 343-353
      9, 2, 9, 3, 9, 3, 3, 3, 3, 9, 3,
      # 354-364
      9, 9, 3, 9, 3, 3, 3, 9, 3, 9, 3,
      # 365-375
      3, 9, 9, 3, 9, 9, 9, 3, 3, 9, 3,
      # 376-386
      9, 3, 3, 9, 3, 3, 3, 3, 3, 9, 3,
      # 387-397
      9, 3, 3, 9, 3, 9, 3, 9, 3, 3, 3,
      # 398-408
      9, 9, 3, 3, 1, 9, 9, 9, 9, 2, 1,
      # 409-419
      9, 9, 2, 1, 9, 2, 1, 1, 3, 2, 1,
      # 420-430
      9, 3, 9, 3, 2, 3, 3, 2, 3, 3, 2,
      # 431-441
      3, 2, 3, 2, 3, 2, 3, 2, 9, 3, 2,
      # 442-452
      9, 9, 3, 9, 9, 3, 3, 9, 3, 9, 3,
      # 453-463
      3, 3, 9, 3, 9, 3, 9, 9, 9, 9, 9,
      # 464-476
      9, 9, 3, 3
    )
    diagnosis.binom <- diagnosis
    diagnosis.binom[which(diagnosis.binom != 9)] <- 1
    # Part 2.10.2: Create factors for each diagnosis
    diagnosis.fact <- factor(diagnosis, ordered = FALSE)
    diagnosis.binom.fact <- factor(diagnosis.binom, ordered = FALSE)
    
  # Part 2.11: Alpha computation
  # Perform Bonferroni correction
  get_alpha <- function(alpha_orig, n) {
    return(alpha_orig/n)
  }
  alpha <- 0.05

  # Part 2.11: Check the normality of gene expressions per patient
    # Part 2.11.1: Pick 5000 features; 5000 is the limit for shapiro.test()
    testgenes.patient <- data[sample(rownames(data), 5000), ]
    # Part 2.11.2: Perform tests on the column (patient)
    p.vals.patient <- rbind(
        apply(testgenes.patient, 2, function(x) shapiro.test(as.matrix(x))$p.value), 
        apply(testgenes.patient, 2, function(x) kruskal.test(as.numeric(x) ~ as.factor(x))$p.value)
      )
    row.names(p.vals.patient) <- c("Normality: Shapiro-Wilk", "Population Equality: Kruskal-Wallis")
    p.vals.patient <- as.data.frame(p.vals.patient)
    # Part 2.11.3: Check which of the patients have non-normal distribution of gene expressions
    # Ho: The gene expressions are normally distributed.
    # Ha: The gene expressions are not normally distributed.
    print(colnames(p.vals.patient[1, which(p.vals.patient[1, ] < alpha)]))
    # Since all patients' gene expressions have p-values less than alpha, we can safely reject Ho and accept Ha: 
    # the gene expressions of each patient are not normally distributed.

  # "In case the gene expression values over the patients are non-normally distributed one may want 
  # to subtract the median and divide by the MAD." (Krijnen)
  data.pp <- data1 <- e[, colnames(p.vals.patient[1, which(p.vals.patient[1, ] < alpha)])]   # Get only those with
  mads <- apply(exprs(data1), 2, mad)
  meds <- apply(exprs(data1), 2, median)
  dat <- sweep(exprs(data1), 2, meds)
  exprs(data.pp) <- sweep(dat, 2, mads, FUN="/")
  boxplot(exprs(data.pp))
  
  # Remove rawData variable to save space
  rm(rawData, celFiles, data1)