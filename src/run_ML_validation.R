# function
source("src/r-function.R")

# intra-extra analysis for top hub gene
total_keyhub <- key_hub_gene$total_keyhub
total_keyhub_list <- key_hub_gene$total_keyhub_merge
clinical_trait <- key_hub_gene$clinical_trait
geneExpression <- key_hub_gene$network$deg

# gene selection
selected_gene <- gene_selection(total_keyhub_list = total_keyhub_list)




ml_validation <- function(pr_name, selected_gene, time_stamp){
  
  log_save <- paste("ML_LOG", pr_name, time_stamp, sep = "/")
  dir.create(log_save, recursive = T, showWarnings = FALSE)
  
  trait_names <- selected_gene %>% names()
  
  for(trait_name in trait_names){
    Y_col_name <- trait_name
    DF <- clinical_trait %>% 
      select(Y_col_name) %>%
      rownames_to_column(var = "sample") %>% 
      inner_join(x = ., y = geneExpression %>% 
                   rownames_to_column(var = "sample") %>% 
                   select(sample, all_of(selected_gene[[Y_col_name]])),
                 by = "sample") %>% 
      column_to_rownames(var = "sample")
    
    # X, y 분리
    # y_df <- DF %>% select(-all_of(total_keyhub_list[[Y_col_name]]))
    # x_df <- DF %>% select_at(all_of(total_keyhub_list[[Y_col_name]]))
    
    # gene selection
    lasso_coef <- feature_selection_LASSO(x_df, y_df)
    lasso_selection_gene <- x_df %>% select(which(abs(lasso_coef) > 0)) %>% colnames()
    
    gene_selection_list[[trait_name]] <- lasso_selection_gene
  }
  
  # venn diagram
  seleted_gene_intersection_plot(gene_selection_list, log_save)
  
  
  return(gene_selection_list)
  
}











