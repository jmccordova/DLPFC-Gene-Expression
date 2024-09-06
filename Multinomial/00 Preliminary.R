# Part 0: Preliminary
# Part 0.1: Indicate your home folder
userdir <- "E:/jmcco/"
#userdir <- "C:/Users/jmcco/"

# Part 0.2: Indicate the directory where you put the GEO dataset
datadir <- "Downloads/BNF 300.2 Data/GSE208338_RAW/"
datadir <- paste(userdir, datadir, sep = "")
# Part 0.3: Indicate the directory of the probe annotations
probedir <- "Downloads/BNF 300.2 Data/HuEx-1_0-st-v2-na36-hg19 Probeset/"
probedir <- paste(userdir, probedir, sep = "")

# Part 0.4: Sets the location of the data to be used and where the packages should be put
setwd(datadir)
install.packages("rstudioapi"); library("rstudioapi");
exportdir <- paste(dirname(rstudioapi::getSourceEditorContext()$path), "/Export", sep = "")
package_loc <- paste(datadir, "lib", sep = "")
