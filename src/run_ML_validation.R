# function
source("src/r-function.R")

# load RData
if(!exists("key_hub_gene")){
  load(paste0(base_dir, "/Step3_KEY_HUB_GENE.RData"))
}

# intra-extra analysis for top hub gene

module_candidate <- list()
for(m_name in names(key_hub_gene$intra_module)){
  total_keyhub <- key_hub_gene$intra_module[[m_name]]
  clinical_trait <- key_hub_gene$clinical_trait
  geneExpression <- key_hub_gene$network$deg
  
  # gene selection
  selected_gene <- gene_selection(base_dir = base_dir,
                                  total_keyhub_list = total_keyhub, 
                                  over_sampling = TRUE) %>% 
    compact()
  
  # save(selected_gene, file = paste0(base_dir, "/Step4_SELECTED_GENE.RData"))
  
  final_candidate_gene <- ml_validation(base_dir = base_dir, selected_gene = selected_gene,
                                        over_sampling=TRUE, cv = TRUE, module_name = m_name) %>% 
    auc_cutoff(sg_list = ., selected_gene = selected_gene, auc_cutoff = 0.75)
  
  # final write
  module_candidate[[m_name]] <- tibble(GENE_NAME = final_candidate_gene$auc_cutoff) %>% 
    mutate(TRAIT = paste0(final_candidate_gene$trait, collapse = ";"),
           MODULE = m_name) %>% 
    arrange(GENE_NAME)
}

save(module_candidate, file = paste0(base_dir, "/Step4_MODULE_SELECTED_GENE.RData"))
