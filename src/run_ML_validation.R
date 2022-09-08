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

for(m_name in names(key_hub_gene$intra_module)[6:7]){
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
    auc_cutoff(sg_list = ., selected_gene = selected_gene, auc_cutoff = 0.75)
  
  # final write
  module_candidate[[m_name]] <- tibble(GENE_NAME = final_candidate_gene$auc_cutoff) %>% 
    mutate(TRAIT = paste0(final_candidate_gene$trait, collapse = ";"),
           MODULE = m_name) %>% 
    arrange(GENE_NAME)
}
save(module_candidate, file = paste0(base_dir, "/Step4_MODULE_SELECTED_GENE.RData"))


# STRING, Survival analysis
string_filtering <- module_candidate %>% bind_rows() %>% 
  pull(1) %>% 
  string_network() %>% 
  filter(score >= 0.7) 
string_filtering <- c(string_filtering$preferredName_A, string_filtering$preferredName_B) %>% unique()

# Survival analysis
survival_analysis <- function(base_dir, geneExpression){
  log_save <- paste(base_dir, "OS", sep = "/")
  dir.create(paste(base_dir, "OS", sep = "/"), showWarnings = FALSE, recursive = TRUE)
  
  suv_exp <- geneExpression %>% rownames_to_column(var = "sample") %>% 
    select(sample, all_of(string_filtering))
  deg_list <- suv_exp %>% colnames() %>% .[-1]
  
  survival_trait <- read_delim(paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2F",pr_name,"_survival.txt"),
                               delim = "\t", show_col_types = FALSE, progress = FALSE) %>% 
    dplyr::select(sample, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time)
  suv_exp_trait <- inner_join(x = suv_exp, y = survival_trait, by = "sample")
  
  os_list <- list()
  for(gene_name in deg_list){
    print(gene_name)
    exp_median_group <- suv_exp_trait %>% dplyr::mutate(group = ifelse(.[[gene_name]] > median(.[[gene_name]]), "H", "L"))
    exp_median_group$group <- exp_median_group$group %>% factor(labels = c("High-exp", "Low-exp"))
    
    sfit <- survfit(Surv(OS.time, OS) ~ group, data = exp_median_group)
    sdf <- survdiff(Surv(OS.time, OS) ~ group, data = exp_median_group)
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    
    os_list[[gene_name]] <- tibble(GENE = gene_name, p.value = p.val)
    
    ggsurvplot(sfit, title = paste0(gene_name, "-Expression High-Low(Median)"), pval = T)
    
    ggsave(filename = paste0(log_save, "/", gene_name, "-OS.png"))
    dev.off()
  }
  os_list %>% bind_rows() %>% return()
}
