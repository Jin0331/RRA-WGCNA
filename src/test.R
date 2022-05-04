# install.packages("reticulate")
library(reticulate)

# load py-function
use_condaenv(condaenv = "geo-py")
source_python("src/py-function.py")

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

DF %>% rownames_to_column(var = "sample_barcode") %>% write_delim(file = "test.txt", delim = "\t", )

y_df <- DF %>% select_at(1)
x_df <- DF %>% select_at(-1)

lasso_coef <- feature_selection_LASSO(x_df, y_df)
svm_rfe_binary <- feature_selection_svm_rfecv(x_df, y_df)

lasso_validation_hub_gene <- x_df[ ,which(lasso_coef > 0)] %>% colnames()

svm_rfe_validation_hub_gene <- x_df[ ,which(svm_rfe_binary == TRUE)] %>% colnames()
intersect_gene_selected <- intersect(lasso_validation_hub_gene, svm_rfe_validation_hub_gene)

# intersection DF
DF %>% select(1, all_of(interse))


# intersect venn
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

ggVennDiagram(x = list(LASS = lasso_validation_hub_gene, `SVM-RFE` = svm_rfe_validation_hub_gene)) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
