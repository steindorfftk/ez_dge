#Set working directory
setwd("~/Desktop/Coding/deaR")

#Install and load required libraries
#if (!require("BiocManager", quietly = TRUE))
#	install.packages("BiocManager")
#BiocManager::install("edgeR")
library("edgeR")
# BiocManager::install("sva")
#library("sva")
# BiocManager::install("org.Hs.eg.db")
#library("org.Hs.eg.db")

#01 - readDGE
cond1_dir <- "input/cond1" #Set condition 1 directory path
cond2_dir <- "input/cond2" #Set condition 2 directory path
cond1_files <- list.files(cond1_dir, pattern = "\\.tabular$", full.names = TRUE) #Save cond1 files
cond2_files <- list.files(cond2_dir, pattern = "\\.tabular$", full.names = TRUE) #Save cond2 files
all_files <- c(cond1_files, cond2_files) #Combine sets of files
DG <- readDGE(all_files, header=T) #Run readDGE

#02 - Create data frame
#02.1 - Load files
# Load files for cond1
for (i in seq_along(cond1_files)) {
  assign(paste0("cond_1_", i), read.table(cond1_files[i], sep = '\t', header = TRUE))
}

# Load files for cond2
for (i in seq_along(cond2_files)) {
  assign(paste0("cond_2_", i), read.table(cond2_files[i], sep = '\t', header = TRUE))
}

#02.2 - Create data frame
all_vars <- ls() #List all variables
cond_vars <- grep("^cond_", all_vars, value = TRUE) #Filter desired data variables names
cond_data <- mget(cond_vars) #Gets data for each file
pre_data_frame_a <- do.call(cbind, cond_data) #Combines data in a data frame
pre_data_frame_b <- pre_data_frame_a[, seq(2, ncol(pre_data_frame_a), by = 2)] #Filters just the second column of the data frame
geneCounts <- as.data.frame(pre_data_frame_b) #Saves all as the geneCounts data frame
row.names(geneCounts) <- cond_1_1[,1] #Recovers row names
oldcolnames <- colnames(geneCounts) #Saves old col names in a vector
newcolnames <- sub("\\..*", "", oldcolnames) #Modifies col names to a cleaner version
colnames(geneCounts) <- newcolnames #Assign new col names
sizeGeneCounts <- dim(geneCounts)
geneCounts <- geneCounts[1:(sizeGeneCounts[1]-5),]



