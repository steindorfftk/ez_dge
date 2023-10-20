#Set working directory
setwd("~/Desktop/Coding/deaR")

#Install and load required libraries
#if (!require("BiocManager", quietly = TRUE))
#	install.packages("BiocManager")
#BiocManager::install("edgeR")
library("edgeR")
# BiocManager::install("sva")
library("sva")
# BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

#Select files
cond1_dir <- "input/cond1" #Set condition 1 directory path
cond2_dir <- "input/cond2" #Set condition 2 directory path
cond1_files <- list.files(cond1_dir, pattern = "\\.tabular$", full.names = TRUE) #Save cond1 files
cond2_files <- list.files(cond2_dir, pattern = "\\.tabular$", full.names = TRUE) #Save cond2 files
all_files <- c(cond1_files, cond2_files) #Combine sets of files
DG <- readDGE(all_files, header=T)







