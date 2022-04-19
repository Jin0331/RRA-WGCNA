# load library
library(GEOquery)
library(limma)
library(tidyverse)

# function
source("src/function.R")


multiple_limma <- list()


# geoquery
gse_name <- "GSE14520"
geneExpression <- getGeneExpressionFromGEO(datasetGeoCode = gse_name, 
                                           retrieveGeneSymbols = TRUE, 
                                           verbose = TRUE)


gse <- getGEO(gse_name,GSEMatrix=TRUE)
gse <- gse[[1]]

eset <- exprs(gse)
fset <- fData(gse)

colnames(fset) #"Symbol" 또는 "Gene Symbol"가 있는지 확인
rownames(eset) <- fset[,"Gene Symbol"]

# sample selection 
pset <- phenoData(gse)
View(pset@data)

# tissue:ch1, source_name_ch1
sample_select <- pset@data %>% 
  # filter(str_detect(characteristics_ch1, "HCC")) %>% 
  pull(`source_name_ch1`)
grp <- pset@data %>% 
  # filter(str_detect(characteristics_ch1, "HCC")) %>% 
  pull(`source_name_ch1`) %>% 
  lapply(X = ., FUN = tolower) %>% 
  unlist() %>% 
  as.factor()

# eset <- eset[, sample_select]


# duplicate remove probe / gene
probe_MAD <- apply(eset,1,mad) %>% 
  tibble(id = seq(length(.)), gene_name = names(.), MAD = .) %>% 
  arrange(desc(MAD))
probe_MAD_dup <- probe_MAD[which(!duplicated(probe_MAD$gene_name)),] %>% 
  filter(gene_name != "") %>% 
  arrange(id)
eset <- eset[probe_MAD_dup$id, ]

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
