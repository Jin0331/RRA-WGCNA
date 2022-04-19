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
geneExpression <- gse_data$gene_expression
pheno <- gse_data$pheno@data
gse_data$pheno@varMetadata
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

robustdeg_ge %>% t() %>% View()

# RUN WGCNA

# gene expression ----
# checking missing value for Gene and Samples
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

# clinical trait ----
featch_suppl <- getGEOSuppFiles(gse_name, fetch_files = FALSE) # 확인 후, load
clinical_trait <- read_delim(file = "RData/GSE14520_Extra_Supplement.txt", delim = "\t")



