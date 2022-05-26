# function
source("src/r-function.R")

# load RData
if(!exists("key_hub_gene")){
  load(paste0(base_dir, "/Step3_KEY_HUB_GENE.RData"))
}

# intra-extra analysis for top hub gene
total_keyhub <- key_hub_gene$total_keyhub
total_keyhub_list <- key_hub_gene$total_keyhub_merge
clinical_trait <- key_hub_gene$clinical_trait
geneExpression <- key_hub_gene$network$deg

# gene selection
selected_gene <- gene_selection(base_dir = base_dir,
                                total_keyhub_list = total_keyhub_list, 
                                over_sampling = TRUE)
save(robustdegs, file = paste0(base_dir, "/Step4_SELECTED_GENE.RData"))

final_candidate_gene <- ml_validation(base_dir = base_dir, selected_gene = selected_gene,
                                       over_sampling=TRUE, cv = TRUE) %>% 
  auc_cutoff(sg_list = ., auc_cutoff = 0.7)

# final write
final_candidate_gene %>% 
  tibble(GENE_NAME = .) %>% 
  arrange(GENE_NAME) %>% 
  write_delim(file = paste0(base_dir, "/",pr_name, "_final_candidate.csv"), delim = ",")












