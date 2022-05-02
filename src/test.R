# install.packages("reticulate")
library(reticulate)

use_condaenv(condaenv = "multiomics-cpu")
source_python("src/py-function.py")

robustdeg_ge # TCGA gene expression  - count
data_trait

# intra-extra analysis for top hub gene
top_hub_gene <- extra_analysis_hub %>% pull(preferredName_A) %>% unique()

# data-preprocessing for X, Y making
Y_col_name <- "sample_type"
DF <- data_trait %>% 
  select(starts_with(Y_col_name)) %>% 
  rownames_to_column(var = "sample") %>% 
  inner_join(x = ., y = robustdeg_ge %>% rownames_to_column(var = "sample"),
             by = "sample") %>% 
  column_to_rownames(var = "sample") %>% 
  select(1, all_of(top_hub_gene))

y_df <- DF %>% select_at(1)
x_df <- DF %>% select_at(-1)



lasso_coef <- df_test(x_df, y_df)
lasso_coef[lasso_coef > 0]

pkl_path, raw_path, cancer_type

cancer_fpkm <- load_tcga_dataset(pkl_path = "PKL/", raw_path = "RAW_DATA/", cancer_type = "LIHC") %>% 
  rownames_to_column(var = "id")

gene_probe_mapping <- read_delim(file = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap",
                                 delim = "\t") %>% 
  select(id, gene) %>% 
  inner_join(x = ., y = cancer_fpkm, by = "id") %>% 
  select(-id)

gene_probe_mapping_tr <- gene_probe_mapping %>% 
  mutate_if(is.numeric, .funs = function(value){
    2^value - 0.001 %>% return()
  }) %>% 
  mutate_if(is.numeric, .funs = function(value){
    ifelse(value < 0, 0, value) %>% return()
  }) %>% 
  mutate_if(is.numeric, .funs = function(value){
    log2(value + 1) %>% return()
  }) %>% as.matrix.data.frame()



