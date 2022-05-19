# function
source("src/r-function.R")

# intra-extra analysis for top hub gene
total_keyhub <- key_hub_gene$total_keyhub
total_keyhub_list <- key_hub_gene$total_keyhub_merge
clinical_trait <- key_hub_gene$clinical_trait
geneExpression <- key_hub_gene$network$deg


# gene selection


gene_selection <- function(total_keyhub_list){
  trait_name <- total_keyhub_list %>% names()
  gene_selection_list <- list()
  
  for(trait_name in total_keyhub_list %>% names()){
    Y_col_name <- trait_name
    DF <- clinical_trait %>% 
      select(Y_col_name) %>%
      rownames_to_column(var = "sample") %>% 
      inner_join(x = ., y = geneExpression %>% 
                   rownames_to_column(var = "sample") %>% 
                   select(sample, all_of(total_keyhub_list[[Y_col_name]])),
                 by = "sample") %>% 
      column_to_rownames(var = "sample")
    
    # X, y 분리
    y_df <- DF %>% select(-all_of(total_keyhub_list[[Y_col_name]]))
    x_df <- DF %>% select_at(all_of(total_keyhub_list[[Y_col_name]]))
    
    # gene selection
    lasso_coef <- feature_selection_LASSO(x_df, y_df)
    lasso_selection_gene <- x_df %>% select(which(abs(lasso_coef) > 0)) %>% colnames()
    
    gene_selection_list[[trait_name]] <- lasso_selection_gene
  }
  
  return(gene_selection_list)
  
}



# venn diagram
ggVennDiagram(x = gene_selection_list) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")









