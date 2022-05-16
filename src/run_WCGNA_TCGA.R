# function
source("src/r-function.R")

# LOCAL VARIABLE
# rdata_path <- "RData/"
load(file = "RData/HCC_GEO_RobustDEGs_norm.RData")
pr_name <- "LIHC"

# weighted gene co-expression network costruction
network_list <- list()
max_module_cnt <- 0
max_module <- NULL

# module cnt test
mch_test <- purrr::quietly(mergecutheight_test)(pr_name = "LIHC", robustdegs = robustdegs)
network <- network_preprocessing(pr_name = "LIHC", robustdegs = robustdegs, 
                                 mch = mch_test$result$max_mch)

module_cluster_plot(network = network)

# network variable
moduleLabels <- network[[2]]$colors
moduleColors <-  labels2colors(network[[2]]$colors)
MEs <- network[[2]]$MEs;
geneTree <- network[[2]]$dendrograms[[1]]

# Calculate ME
MEs0 <- moduleEigengenes(network[[1]], moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
    

find_key_modulegene <- function(network, MEs, select_clinical=NULL, mm=0.85, gs=0.2){
  
  # variable
  expression_sample <- rownames(network[[1]])
  nGenes <- ncol(network[[1]])
  nSamples <- nrow(network[[1]])
  
  # clinical feature
  # default_clinical <- c('sample_type', 'OS', 'OS.time', 'age_at_initial_pathologic_diagnosis', 'pathologic_T', 
  #                       'pathologic_M', 'pathologic_N', 'pathologic_stage', 'child_pugh_classification_grade', 
  #                       'fibrosis_ishak_score')
  
  default_clinical <- c('sample_type', 'OS.time', 'pathologic_T', 'pathologic_M', 'pathologic_N', 'pathologic_stage')
  # default_clinical <- c('OS.time', 'pathologic_T', 'pathologic_M', 'pathologic_N', 'pathologic_stage')
  
  use_clinical <- c(default_clinical, select_clinical)
  
  # UCSCXena clinical
  clinical_trait <- read_delim("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FLIHC_clinicalMatrix",
                               delim = "\t", show_col_types = FALSE, progress = FALSE)
  survival_trait <- read_delim("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FLIHC_survival.txt",
                               delim = "\t", show_col_types = FALSE, progress = FALSE) %>% 
    dplyr::select(sample, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time)
  
  clinical_trait <- left_join(x = clinical_trait, y = survival_trait, by = c("sampleID" = "sample")) 
  
  # clinical trait preprocessing
  traitRows <- match(expression_sample, clinical_trait$sampleID)
  data_trait <- clinical_trait[traitRows, ] %>% 
    column_to_rownames(var = "sampleID") %>% 
    dplyr::select(all_of(use_clinical)) %>% 
    # mutate(sample_type = ifelse(sample_type == "Primary Tumor", 1, 0)) %>%  # sample type 한정
    mutate_if(is.character, as.factor) %>% 
    mutate_all(as.numeric)
  data_trait[is.na(data_trait)] <- 0
  
  # module relation calculation 
  moduleTraitCor <-  WGCNA::cor(MEs, data_trait, use = "p")
  moduleTraitPvalue <-  corPvalueStudent(moduleTraitCor, nSamples)
  module_trait_plot(moduleTraitCor = moduleTraitCor, moduleTraitPvalue = moduleTraitPvalue,
                    data_trait = data_trait, MEs = MEs)
  
  
  # top 3 sign. module
  signModule <- lapply(X = 1:nrow(moduleTraitPvalue), FUN = function(row_index){
    trait_cnt <- moduleTraitPvalue[row_index, ] %>% .[. < 0.05] %>% length()
    row_name <- rownames(moduleTraitPvalue)[row_index] %>% 
      substring(., 3)
    if(str_detect(row_name, "grey"))
      return(NULL)
    else
      return(tibble(MM = row_name, signTrait = trait_cnt))
  }) %>% 
    bind_rows() %>% arrange(desc(signTrait)) %>% 
    dplyr::pull(1) %>% 
    .[1:3]
  
  gene_module_key <- network[[2]]$colors %>% names() %>% tibble(gene = .) %>% 
    bind_cols(., tibble(module = moduleColors)) %>% 
    filter(module %in% signModule)
  
  # Module membership
  MM <- as.data.frame(cor(network[[1]], MEs, use = "p")) 
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(MM), nSamples))
  MM <- MM %>% 
    rownames_to_column(var = "gene")
  
  # Gene significance
  # gene significance
  GS <- lapply(X =  use_clinical, FUN = function(s_type){
    trait <- data_trait[s_type]
    names(trait) <- s_type
    
    # gene significance
    geneTraitSignificance <- as.data.frame(cor(network[[1]], trait, use = "p")) %>%
      mutate_all(abs) %>%
      rownames_to_column("gene")
    
  }) %>% purrr::reduce(., inner_join, by = "gene") %>%
    column_to_rownames("gene") %>%
    bind_cols(., apply(., 1, median) %>% tibble(GS_median = .)) %>%
    # select(GS_median) %>%
    rownames_to_column("gene")
  
  # module 별 clinical trait 각각 calulation
  intra_module <- lapply(X = signModule, FUN = function(sig_module){
    module_gene <- gene_module_key %>% 
      filter(module == sig_module) %>% 
      dplyr::pull(gene)
    
    module_MM <- MM %>% 
      select(gene, MM = ends_with(sig_module)) %>% 
      filter(gene %in% module_gene, abs(MM) > mm)
    module_MM_gene <- module_MM %>% dplyr::pull(gene)
    
    module_MM_GS_filtered <- lapply(X = use_clinical, FUN = function(t){
      GS %>% 
        select(gene, GS = t) %>% 
        filter(gene %in% module_MM_gene, abs(GS) > gs)})
    names(module_MM_GS_filtered) <- use_clinical
    
    return(module_MM_GS_filtered)
  })
  
  total_keyhub <- list()
  for(index in use_clinical){
    tmp <- lapply(X = intra_module, FUN = function(trait){
      trait[[index]]
    }) %>% bind_rows() %>% dplyr::pull(gene)
    
    if(length(tmp) != 0)
      total_keyhub[[index]] <- tmp
  }

  p <- ggVennDiagram::ggVennDiagram(x = total_keyhub) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none")
  
  ggsave(p, filename = "test.png")
  
  total_keyhub %>% unlist() %>% unname() %>% 
    return()
  
}
# modules Hub gene
intra_analysis_hub <- intra_module %>% 
  bind_rows() %>% 
  pull(1)
    
# extramodular analysis
# extra_analysis_hub <- string_network(hub_gene = intra_analysis_hub) %>% 
#   filter(score > 0.9) 
# 
# ## cytoscape (?)
# extra_analysis_hub %>% write_delim(file = "Cytoscape/tophubgene.txt", delim = "\t")

# collection of plot
{
  # plot - GS & MM correlation 
  {
    sizeGrWindow(7, 7);
    par(mfrow = c(1,1));
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "GS for sample type",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  }
