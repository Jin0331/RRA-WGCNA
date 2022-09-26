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
string_filtering <- string_analysis(mc = module_candidate, base_dir = base_dir) %>% 
  lapply(X = ., FUN = function(df){
    df %>% mutate(STRING = TRUE)
  }) %>% bind_rows()

# Survival analysis
survival_filtering <- survival_analysis(base_dir = base_dir, 
                                        geneExpression = geneExpression, 
                                        mc = module_candidate) %>% 
  bind_rows()

# DNA metylation analysis
methylation_filtering <- methylation_analysis(base_dir = base_dir)
methylation_filtering_ <- lapply(X = names(methylation_filtering), FUN = function(method){
    GENE_NAME <- methylation_filtering[[method]]$Symbol
    P.adjust <- methylation_filtering[[method]]$Pe
    
    tmp_df <- tibble(GENE_NAME = GENE_NAME, Methylated = method, P.adjust = P.adjust) %>% 
      distinct(GENE_NAME, .keep_all = TRUE)
    colnames(tmp_df) <- c("GENE_NAME", paste0('Methylated_', method),
                          paste0('P.adjust_', method))
    return(tmp_df)
})
names(methylation_filtering_) <- names(methylation_filtering)


# analysis merge
module_candidate %>% bind_rows() %>% 
  left_join(x = ., y = string_filtering, by = c("GENE_NAME", "TRAIT", "MODULE")) %>% 
  left_join(x = ., y = survival_filtering, by = c("GENE_NAME", "TRAIT", "MODULE")) %>% 
  left_join(x = ., y = methylation_filtering_$hypo, by = "GENE_NAME") %>% 
  left_join(x = ., y = methylation_filtering_$hyper, by = "GENE_NAME") %>% 
  unite(Methylated, Methylated_hypo, Methylated_hyper, sep = "/", na.rm = TRUE) %>% 
  unite(Methylated.Padj, P.adjust_hypo, P.adjust_hyper, sep = "/", na.rm = TRUE)

# Enrichment analysis
ora_go_kegg(geneName = module_candidate$turquoise$GENE_NAME, module_name = "turquoise",
            base_dir = base_dir)

lapply(X = names(module_candidate), FUN = function(module_name){
  ora_go_kegg(geneName = module_candidate[[module_name]]$GENE_NAME,
              module_name = module_name,
              base_dir = base_dir)
})


