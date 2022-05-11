# load library


# function
source("src/r-function.R")

# 1 - gene expression, 2 - phenotype
multiple_limma <- list()
gse_name <- "GSE36376"

# file load
if(file.exists(paste0("GSE/", gse_name, ".RData"))){
  load(file = paste0("GSE/", gse_name, ".RData"))
} else {
  gse_data <- getGeneExpressionFromGEO(datasetGeoCode = gse_name, 
                                       retrieveGeneSymbols = TRUE, 
                                       verbose = TRUE)
  save(gse_data, file = paste0("GSE/", gse_name, ".RData"))
}

geneExpression <- gse_data$gene_expression
boxplot(geneExpression[1:50, 1:50])
geneExpression <- log2(geneExpression)

pheno <- gse_data$pheno@data

# sample selection
View(pheno)
grp <- pheno %>% 
  # filter(str_detect(`source_name_ch1`, "HCC")) %>%
  # pull(`characteristics_ch1`) %>%
  # pull(`characteristics_ch1.1`) %>%
  pull(`source_name_ch1`) %>%
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
save(multiple_limma, file = "RData/HCC_GEO_integrated_norm.RData")

# rra & robust DEGs
robust_degs <- rra_analysis(m_list = multiple_limma, logfc = 0, fdr = 0.05)

