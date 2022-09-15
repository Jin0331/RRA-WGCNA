# function
source("src/r-function.R")

# load RData
if(!exists("key_hub_gene")){
  load(paste0(base_dir, "/Step3_KEY_HUB_GENE.RData"))
}

# Machine learning validation
module_candidate <- list()
clinical_trait <- key_hub_gene$clinical_trait
geneExpression <- key_hub_gene$network$deg

for(m_name in names(key_hub_gene$intra_module)){
  print(m_name)
  total_keyhub <- key_hub_gene$intra_module[[m_name]]
  
  # gene selection
  selected_gene <- gene_selection(base_dir = base_dir,
                                  total_keyhub_list = total_keyhub, 
                                  over_sampling = TRUE) %>% compact()
  final_candidate_gene <- ml_validation(base_dir = base_dir, selected_gene = selected_gene,
                                        over_sampling=TRUE, cv = TRUE, module_name = m_name) %>% 
    auc_cutoff(sg_list = ., selected_gene = selected_gene, auc_cutoff = 0.7)
  
  # final write
  module_candidate[[m_name]] <- tibble(GENE_NAME = final_candidate_gene$auc_cutoff) %>% 
    mutate(TRAIT = paste0(final_candidate_gene$trait, collapse = ";"),
           MODULE = m_name) %>% 
    arrange(GENE_NAME)
}

save(module_candidate, file = paste0(base_dir, "/Step4_MODULE_SELECTED_GENE.RData"))


# STRING analysis
string_filtering <- string_analysis(mc = module_candidate, base_dir = base_dir)

# Survival analysis
survival_filtering <- survival_analysis(base_dir = base_dir, 
                                        geneExpression = geneExpression, 
                                        mc = module_candidate)

# DNA metylation analysis
methylation_filtering <- methylation_analysis(base_dir = base_dir)

