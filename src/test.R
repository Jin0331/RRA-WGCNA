# install.packages("reticulate")
library(reticulate)

use_condaenv(condaenv = "multiomics-cpu")
source_python("src/py-function.py")

robustdeg_ge # TCGA gene expression  - count
data_trait

# intra-extra analysis for top hub gene
# top_hub_gene <- extra_analysis_hub %>% pull(preferredName_A) %>% unique()
top_hub_gene <- intra_analysis_hub %>% unique()

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



lasso_coef <- feature_selection_LASSO(x_df, y_df)
svm_rfe_binary <- feature_selection_svm_rfecv(x_df, y_df)

lasso_validation_hub_gene <- x_df[,which(lasso_coef > 0)] %>% colnames()



