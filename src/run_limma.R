# load library
library(GEOquery)
library(limma)
library(tidyverse)

# function
source("src/function.R")

# 1 - gene expression, 2 - phenotype
multiple_limma <- list()
gse_name <- "GSE57957"
gse_data <- getGeneExpressionFromGEO(datasetGeoCode = gse_name, 
                                           retrieveGeneSymbols = TRUE, 
                                           verbose = TRUE)
save(gse_data, file = paste0("GSE/", gse_name, ".RData"))


geneExpression <- gse_data$gene_expression
pheno <- gse_data$pheno@data

# sample selection
View(pheno)
grp <- pheno %>% 
  pull(`description`) %>% 
  lapply(X = ., FUN = tolower) %>% 
  unlist() %>% 
  as.factor()

grp %>% unique()
design <- model.matrix(~0 + grp)
head(grp)
colnames(design) <- c("NT","TP") # 데이터에 맞추어 manual로 설정해야 됨

# Limma
multiple_limma[[gse_name]] <- run_limma(ge = geneExpression, de = design)

# save
save(multiple_limma, file = "GEO_integrated.RData")
