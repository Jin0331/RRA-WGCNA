# function
source("src/r-function.R")

# intra-extra analysis for top hub gene
total_keyhub <- key_hub_gene$total_keyhub
total_keyhub_list <- key_hub_gene$total_keyhub_merge
clinical_trait <- key_hub_gene$clinical_trait
geneExpression <- key_hub_gene$network$deg

# gene selection
selected_gene <- gene_selection(total_keyhub_list = total_keyhub_list, pr_name = pr_name, time_stamp = time_stamp)
final_candidatte_gene <- ml_validation(selected_gene = selected_gene, pr_name = pr_name, time_stamp = time_stamp) %>% 
  auc_cutoff(sg_list = ., auc_cutoff = 0.7)














