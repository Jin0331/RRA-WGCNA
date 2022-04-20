library(WGCNA)
library(GEOquery)
library(tidyverse)

# function
source("src/function.R")

# load data
load(file = "RData/RobustDEGs.RData")
robustdegs <- up_down_rra_gene %>% pull(1)

# base data
gse_name <- "GSE14520"
gse_data <- getGeneExpressionFromGEO(datasetGeoCode = gse_name, 
                                     retrieveGeneSymbols = TRUE, 
                                     verbose = TRUE)
save(gse_data, file = paste0("GSE/", gse_name, ".RData"))

geneExpression <- gse_data$gene_expression
pheno <- gse_data$pheno@data

# RUN WGCNA
# preprocessing----
# checking missing value for Gene and Samples
robustdeg_ge <- lapply(X = robustdegs, FUN = function(deg){
  error <- FALSE
  tryCatch(
    expr = {
      tmp <- geneExpression[deg, ]
    },
    error = function(e) {
      error <<- TRUE
    }
  )
  if(error){
    return(NULL)
  } else {
    tmp <- as.matrix(tmp) %>% t()
    rownames(tmp) <- deg
    return(tmp)
  }}) %>% do.call(rbind, .) %>% 
  t() %>% 
  as.data.frame()

gsg <- goodSamplesGenes(robustdeg_ge, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(robustdeg_ge)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(robustdeg_ge)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  robustdeg_ge <-  robustdeg_ge[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <-  hclust(dist(robustdeg_ge), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# clinical trait 
featch_suppl <- getGEOSuppFiles(gse_name, fetch_files = FALSE) # 확인 후, load
clinical_trait <- read_delim(file = "RData/GSE14520_Extra_Supplement.txt", delim = "\t") %>% 
  as.data.frame()

expression_sample <- rownames(robustdeg_ge)
traitRows = match(expression_sample, clinical_trait$Affy_GSM)
data_trait <- clinical_trait[traitRows, -3]
rownames(data_trait) <- clinical_trait[traitRows, 3]


# Automatic network construction and module detection ----
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(robustdeg_ge, powerVector = powers, verbose = 5)

# plot
{
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
}

softPower <-  7 #Chosen in the graphs before
net <- blockwiseModules(datExpr = robustdeg_ge, 
                        power = softPower,
                       # corType = "pearson", 
                       TOMType = "unsigned", 
                       minModuleSize = 30,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.15,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "GSE14520",
                       verbose = 3)

moduleLabels <- net$colors
moduleColors %>% unique()
moduleColors <- labels2colors(net$colors)



nGenes <- ncol(robustdeg_ge)
nSamples <- nrow(robustdeg_ge)
MEs0 = moduleEigengenes(robustdeg_ge, moduleColors)$eigengenes
