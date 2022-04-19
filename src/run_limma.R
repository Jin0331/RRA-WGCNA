# load library
library(GEOquery)
library(limma)
library(tidyverse)

# function
source("src/function.R")


multiple_limma <- list()

# 1 - gene expression, 2 - phenotype
gse_name <- "GSE14520"
gse_data <- getGeneExpressionFromGEO(datasetGeoCode = gse_name, 
                                           retrieveGeneSymbols = TRUE, 
                                           verbose = TRUE)

geneExpression <- gse_data$gene_expression
# rownames(geneExpression) <- geneExpression$GeneSymbol
pheno <- gse_data$pheno@data

# sample selection
View(pheno)

# tissue:ch1, source_name_ch1
# sample_select <- pset@data %>% 
#   filter(str_detect(characteristics_ch1, "HCC")) %>% 
#   pull(`source_name_ch1`)
grp <- pheno %>% 
  # filter(str_detect(characteristics_ch1, "HCC")) %>% 
  pull(`characteristics_ch1`) %>% 
  lapply(X = ., FUN = tolower) %>% 
  unlist() %>% 
  as.factor()

# duplicate remove probe / gene
gene_name <- geneExpression$GeneSymbol
sample_name <- colnames(geneExpression)[1:(length(colnames(geneExpression)) - 1)]
probe_MAD <- geneExpression %>% 
  select(-GeneSymbol) %>% 
  as.matrix() %>% 
  apply(.,1,mad) %>% 
  tibble(id = seq(length(.)), gene_name = gene_name, MAD = .) %>% 
  arrange(desc(MAD))
probe_MAD_dup <- probe_MAD[which(!duplicated(probe_MAD$gene_name)),] %>% 
  filter(gene_name != "") %>% 
  arrange(id)
geneExpression <- geneExpression[probe_MAD_dup$id, ]

geneExpression_dup <- geneExpression %>% 
  select(-GeneSymbol) %>% 
  as.matrix()
rownames(geneExpression_dup) <- geneExpression$GeneSymbol
colnames(geneExpression_dup) <- colnames(geneExpression)[1:(length(geneExpression) - 1)]

# Limma
grp %>% unique()
design <- model.matrix(~0 + grp)
head(grp)
colnames(design) <- c("NT","TP") # 데이터에 맞추어 manual로 설정해야 됨

fit <- lmFit(eset,design)
cont <- makeContrasts(TP-NT,levels=design) # 데이터에 맞추어 manual로 설정해야 됨
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont)
res <- topTable(fit.cont,number=Inf) 
target <- res[res$P.Value <0.05,] %>% 
  rownames_to_column() %>% as_tibble()

# save
multiple_limma[[gse_name]] <- target
save(multiple_limma, file = "GEO_integrated.RData")
