#Definitions
study_name = "example"
organism = "Mus musculus" #Homo sapiens or Mus musculus

#Set working directory
setwd("~/Desktop/Coding/deaR")


#Install and load required libraries
#if (!require("BiocManager", quietly = TRUE))
#	install.packages("BiocManager")
#BiocManager::install("edgeR")
library("edgeR")
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
  assign(paste0("cond_1_", i), read.table(cond1_files[i], sep = '\t', header = FALSE))
}

# Load files for cond2
for (i in seq_along(cond2_files)) {
  assign(paste0("cond_2_", i), read.table(cond2_files[i], sep = '\t', header = FALSE))
}

#02.2 - Create data frame
all_vars <- ls() #List all variables
cond_vars <- grep("^cond_", all_vars, value = TRUE) #Filter desired data variables names
cond_data <- mget(cond_vars) #Gets data for each file
pre_data_frame_a <- do.call(cbind, cond_data) #Combines data in a data frame
pre_data_frame_b <- pre_data_frame_a[, seq(2, ncol(pre_data_frame_a), by = 2)] #Filters just the second column of the data frame
geneCounts <- as.data.frame(pre_data_frame_b) #Saves all as geneCounts data frame
row.names(geneCounts) <- cond_1_1[,1] #Recovers row names
oldcolnames <- colnames(geneCounts) #Saves old col names in a vector
newcolnames <- sub("\\..*", "", oldcolnames) #Modifies col names to a cleaner version
colnames(geneCounts) <- newcolnames #Assign new col names
sizeGeneCounts <- dim(geneCounts)
geneCounts <- geneCounts[1:(sizeGeneCounts[1]-5),]

#03 Rename the columns and set conditions
colNames <- colnames(geneCounts) #Saves the col names
colNames <- gsub('cond_1_', 'ctr', colNames) #Changes 'cond_1_' for 'ctr' in col names
colNames <- gsub('cond_2_', 'case', colNames) #Changes 'cond_2_' for 'case' in col names
colnames(geneCounts) <- colNames #Applies new col names for geneCounts
numCtr <- sum(grepl('ctr', colNames)) #Counts the number of controls
numCase <- sum(grepl('case', colNames)) #Counts the number of cases
condition <- c(rep('ctr', numCtr), rep('case', numCase)) #Sets conditions vector
sampleNames <- colNames #Saves sample names in a vector


#04 Construct Linear Generalized Model 
dge <- DGEList(counts=geneCounts, group=condition) #Create DGEList
design <- model.matrix(~condition+0, data=dge$samples) #Create design matrix
colnames(design) <- gsub("condition","",colnames(design)) #Rename design matrix col names

#05 TMM Normalization
dge <- calcNormFactors(dge) #Calculate normalization factors
plotMDS(dge) #Plots normalization factors
norm_counts <- cpm(dge,log = TRUE, prior.count = 3) #Adds 3 units to each count, calculates counts per million and transform into log.
exp <- as.data.frame(norm_counts) #Saves cpm as data frame
plotNameA <- paste0("output/plots/01_", study_name , "_MDSplot.jpeg" ) #Defines first plot name
jpeg(file=plotNameA, width=5000, height=5000, units="px", res=300) #Defines first plot specs
plotMDS(dge) #Saves first plot as jpeg
dev.off() #Clears plot panel

#06 Dispersion estimation and plotting
disp <- estimateGLMCommonDisp(dge, design) #Estimate commom dispersion
disp <- estimateGLMTrendedDisp(disp, design) #Estimate trended dispersion
disp <- estimateGLMTagwiseDisp(disp, design) #Estimate tagwise dispersion
plotBCV(disp) #Plot dispersions
plotNameB <- paste0("output/plots/02_", study_name , "_Dispersion_BCVplot.jpeg" ) #Defines second plot name
jpeg(file=plotNameB, width=5000, height=5000, units="px", res=300) #Defines second plot specs
plotBCV(disp) #Saves second plot as jpeg
dev.off() #Clears plot panel

#07 SVA Normalization
# BiocManager::install("sva")
library("sva")
design0 <- as.data.frame(design) #Create a second design data frame for SVA
design0 <- model.matrix(~1, data=design0) #Assigns 1 for all conditions in the new data frame 
n.sv <- num.sv(norm_counts,design,method="leek") #Estimate number of latent factors
svobj <- sva(norm_counts,design,design0,n.sv=n.sv) #Identify latent variation factors
svobj.df <- data.frame(svobj$sv) #Saves identified latent factors as data frame
cleaningP <- function(y, design, svaobj,  P=ncol(design)) {
  X=cbind(design,svaobj$sv)
  Hat=solve(t(X)%*%X)%*%t(X)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
} #Defines function for removing latent variation factors as described in GitHub:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0808-5
cleanp <- cleaningP(norm_counts,design,svobj) #Removes latent variation factors from data
pca <- prcomp(t(cleanp)) #Calculates principal components for the new data frame
plot(pca) #Plots the PCA
plotNameC <- paste0("output/plots/03_", study_name , "_PCA_after_SVA.jpeg" ) #Defines third plot name
jpeg(file=plotNameC, width=5000, height=5000, units="px", res=300) #Defines third plot specs
plot(pca) #Saves third plot
dev.off() #Clears plot panel

