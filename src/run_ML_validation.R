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
                                  over_sampling = TRUE) %>% 
    compact()
  
  # save(selected_gene, file = paste0(base_dir, "/Step4_SELECTED_GENE.RData"))
  
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


# 
function(mc, base_dir, score = 400){
  log_save <- paste(base_dir, "ANALYSIS/STRING", sep = "/")
  dir.create(paste(base_dir, "ANALYSIS/STRING", sep = "/"), showWarnings = FALSE, recursive = TRUE)
  
  
  tmp <- lapply(X = names(mc), FUN = function(mc_name){
    g <- mc[[mc_name]]
    proteins_mapped <- rba_string_map_ids(ids = g$GENE_NAME, species = 9606)
    int_net <- rba_string_interactions_network(ids = proteins_mapped$stringId,
                                               species = 9606,
                                               required_score = score) # minimum confidence
    
    rba_string_network_image(ids = proteins_mapped$stringId,
                             image_format = "highres_image",
                             species = 9606,
                             save_image = paste0(getwd(),
                                                 "/", base_dir, "/ANALYSIS/STRING/", mc_name, ".png"),
                             required_score = 400,
                             network_flavor = "confidence", 
                             hide_disconnected_nodes = TRUE)
    
    
    string_filtered_g <- c(int_net$preferredName_A, int_net$preferredName_B) %>% 
      unique() %>% 
      tibble(GENE_NAME = .) %>% 
      inner_join(x = g, y = ., by = "GENE_NAME")
    
    return(string_filtered_g)
  })
  
}

graph_1 <- rba_string_network_image(ids = proteins_mapped$stringId,
                                    image_format = "image",
                                    species = 9606,
                                    save_image = paste0(base_dir, "/ANALYSIS/STRING/", mc_name),
                                    required_score = 400,
                                    network_flavor = "confidence", 
                                    hide_disconnected_nodes = TRUE)

# primary tumor
survival_filtering <- survival_analysis(base_dir = base_dir, 
                                        geneExpression = geneExpression, 
                                        string_filtering = string_filtering)

survival_filtering %>% bind_rows() %>% View()

