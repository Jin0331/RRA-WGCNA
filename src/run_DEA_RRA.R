# function
source("src/r-function.R")

time_stamp <- Sys.time()
pr_name <- "LIHC"
base_dir <- paste("WGCNA_RRA_RESULT", pr_name, time_stamp,sep = "/")
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# ======= MANUAL# ======= 
# 1 - gene expression, 2 - phenotype
multiple_limma <- list()
gse_name <- readline('enter GSE assesion : ')

# file load
if(file.exists(paste0("GSE/", gse_name, ".RData"))){
  load(file = paste0("GSE/", gse_name, ".RData"))
} else {
  dir.create("GSE/", showWarnings = FALSE, recursive = TRUE)
  gse_data <- getGeneExpressionFromGEO(datasetGeoCode = gse_name, 
                                       retrieveGeneSymbols = TRUE, 
                                       verbose = TRUE)
  save(gse_data, file = paste0("GSE/", gse_name, ".RData"))
}

geneExpression <- gse_data$gene_expression
boxplot(geneExpression[1:50, 1:50])

# run log2-transformation
log2_trans <- readline('Run log2-transformation? [y]/[n]  : ')
if(tolower(log2_trans) == "yes" | tolower(log2_trans) == "y"){
  log2_trans <- TRUE
} else {
  log2_trans <- FALSE
}

if(log2_trans){
  geneExpression <- log2(geneExpression)
}

# phenotype selection
pheno <- gse_data$pheno@data
pheno_check <- lapply(X = names(pheno), FUN = function(col_name){
  df <- pheno[col_name]
  df_el <- df %>% pull(1) %>% unique()
  
  if(length(df_el) < 2 | length(df_el) > 5){
    return(NULL) 
  } else {
    df_el <- c(length(df_el), df_el)
    df_el <- paste0(df_el, collapse = " / ")
    tibble(col_name = col_name, factor = df_el) %>% 
      return()
  }
}) %>% compact() %>% bind_rows()
print(pheno_check)


# sample selection
selected_pheno <- readline('enter phenotype : ')
grp <- pheno %>% 
  select(starts_with(selected_pheno)) %>% 
  pull(1) %>% 
  # filter(str_detect(`source_name_ch1`, "HCC")) %>%
  lapply(X = ., FUN = tolower) %>% 
  unlist() %>% 
  as.factor()
print(grp %>% unique())
design <- model.matrix(~0 + grp)

nt_tp_order <- readline("enter NT-TP order : ") %>% 
  str_split(pattern = " ") %>% unlist()
colnames(design) <- nt_tp_order # 데이터에 맞추어 manual로 설정해야 됨

# Limma
multiple_limma[[gse_name]] <- run_limma(ge = geneExpression, de = design)


# ======= RUN DEA & RRA ======= 
# save
save(multiple_limma, file = "RData/HCC_GEO_integrated_norm.RData")

# rra & robust DEGs

load("RData/HCC_GEO_integrated_norm.RData")
robustdegs <- rra_analysis(m_list = multiple_limma, logfc = 0, fdr = 0.05)

