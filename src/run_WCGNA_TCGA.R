# function
source("src/r-function.R")

# LOCAL VARIABLE
# rdata_path <- "RData/"
load(file = "RData/HCC_GEO_RobustDEGs_norm.RData")
pr_name <- "LIHC"

# weighted gene co-expression network costruction
network_list <- list()

# module cnt test
for(mch_value in seq(from = 0.04, to = 0.5, by = 0.02)){
  network <- network_preprocessing(pr_name = "LIHC", robustdegs = robustdegs, mch = 0.2)
  module_cnt <- length(network[[2]]$colors %>% unique()) - 1
  
  network_list[[as.character(mch_value)]] <- list(network, module_cnt)
}

network <- network_preprocessing(pr_name = "LIHC", robustdegs = robustdegs, mch = 0.2)
module_cluster_plot(network = network)

# network variable
moduleLabels <- network[[2]]$colors
moduleColors <-  labels2colors(network[[2]]$colors)
MEs <- network[[2]]$MEs;
geneTree <- network[[2]]$dendrograms[[1]]

# Calculate ME
MEs0 <- moduleEigengenes(network[[1]], moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
    

find_key_modulegene <- function(network, MEs, select_clinical=NULL){
  
  # variable
  expression_sample <- rownames(network[[1]])
  nGenes <- ncol(network[[1]])
  nSamples <- nrow(network[[1]])
  
  # clinical feature
  # default_clinical <- c('sample_type', 'OS', 'OS.time', 'age_at_initial_pathologic_diagnosis', 'pathologic_T', 
  #                       'pathologic_M', 'pathologic_N', 'pathologic_stage', 'child_pugh_classification_grade', 
  #                       'fibrosis_ishak_score')
  
  default_clinical <- c('sample_type', 'OS.time', 'pathologic_T', 'pathologic_M', 'pathologic_N', 'pathologic_stage')
  
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
    dplyr::select(sample_type, all_of(use_clinical)) %>% 
    mutate(sample_type = ifelse(sample_type == "Primary Tumor", 1, 0)) %>%  # sample type 한정
    mutate_if(is.character, as.factor) %>% 
    mutate_all(as.numeric)
  data_trait[is.na(data_trait)] <- 0
  
  # multiple class to binary class
  # traitRows <- match(expression_sample, clinical_trait$sampleID)
  # data_trait <- clinical_trait[traitRows, ] %>%
  #   column_to_rownames(var = "sampleID") %>%
  #   dplyr::select(sample_type, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time,age_at_initial_pathologic_diagnosis,
  #                 pathologic_T, pathologic_M, pathologic_N, pathologic_stage, child_pugh_classification_grade,
  #                 fibrosis_ishak_score)
  # # mutate_if(is.character, as.factor)
  # data_trait[is.na(data_trait)] <- 0
  # ###### using model matrix
  # data_trait <- model.matrix( ~ child_pugh_classification_grade - 1, data_trait) %>%
  #   as.data.frame()
  
  
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
  
  # # gene significance
  # GS <- lapply(X =  use_clinical, FUN = function(s_type){
  #     trait <- data_trait[s_type]
  #     names(trait) <- s_type
  #     
  #     # gene significance
  #     geneTraitSignificance <- as.data.frame(cor(network[[1]], trait, use = "p")) %>% 
  #       mutate_all(abs) %>% 
  #       rownames_to_column("gene")
  #     
  #   }) %>% purrr::reduce(., inner_join, by = "gene") %>% 
  #   column_to_rownames("gene") %>% 
  #   bind_cols(., apply(., 1, median) %>% tibble(GS_median = .)) %>% 
  #   select(GS_median) %>% 
  #   rownames_to_column("gene")
  # 
  # MM <- lapply(X =  use_clinical, FUN = function(s_type){
  #   trait <- data_trait[s_type]
  #   names(trait) <- s_type
  #   
  #   # membership module
  #   geneModuleMembership <- as.data.frame(cor(network[[1]], MEs, use = "p")) %>%
  #     mutate_all(abs) %>% 
  #     rownames_to_column("gene")
  #   
  # }) %>% purrr::reduce(., inner_join, by = "gene") %>% 
  #   column_to_rownames("gene") %>% 
  #   bind_cols(., apply(., 1, median) %>% tibble(MM_median = .)) %>% 
  #   select(MM_median) %>% 
  #   rownames_to_column("gene")
  
  gene_module_key <- network[[2]]$colors %>% names() %>% tibble(gene = .) %>% 
    bind_cols(., tibble(module = moduleColors))
  
  
  gene_module_key %>% filter(module %in% signModule) %>% View()
  
  
  trait <- data_trait[s_type]
  names(trait) <- s_type

  geneTraitSignificance <- as.data.frame(cor(network[[1]], trait, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  
  # Module membership
  geneModuleMembership <- as.data.frame(cor(network[[1]], MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  # intramodular analysis
  # module <- "turquoise"
  geneGS <- geneTraitSignificance %>% rownames_to_column()
  geneMM <- geneModuleMembership %>% rownames_to_column()
  
  intra_module <- lapply(X = signModule, FUN = function(module){
    module_name <- paste0("MM", module)
    moduleGenes <- (moduleColors==module)
    
    geneGS_MM <- inner_join(x = geneGS, y = geneMM, by = "rowname") %>%
      dplyr::select(rowname, starts_with(sig_trait), starts_with(module_name)) %>% 
      .[moduleGenes,]
    colnames(geneGS_MM) <- c("GENE", "GS", "MM")
    
    geneGS_MM %>% 
      filter(abs(GS) > 0.3 & abs(MM) > 0.85) %>% 
      return()
  })
  
  
  names(intra_module) <- signModule
  
  # modules Hub gene
  intra_analysis_hub <- intra_module %>% 
    bind_rows() %>% 
    dplyr::pull(1)
  
  
  }




    

# ## Gene significance
# sig_trait <- "sample_type"
# s_type <- as.data.frame(data_trait$sample_type)
# names(s_type) <- sig_trait
# 
# modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(robustdeg_ge, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# 
# # module membership
# names(geneModuleMembership) <- paste("MM", modNames, sep="")
# names(MMPvalue) <- paste("p.MM", modNames, sep="")
# 
geneTraitSignificance <- as.data.frame(cor(robustdeg_ge, s_type, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
    
# intramodular analysis
# module <- "turquoise"
geneGS <- geneTraitSignificance %>% rownames_to_column()
geneMM <- geneModuleMembership %>% rownames_to_column()

intra_module <- lapply(X = signModule, FUN = function(module){
  module_name <- paste0("MM", module)
  moduleGenes <- (moduleColors==module)
  
  geneGS_MM <- inner_join(x = geneGS, y = geneMM, by = "rowname") %>%
    dplyr::select(rowname, starts_with(sig_trait), starts_with(module_name)) %>% 
    .[moduleGenes,]
  colnames(geneGS_MM) <- c("GENE", "GS", "MM")
  
  geneGS_MM %>% 
    filter(abs(GS) > 0.3 & abs(MM) > 0.85) %>% 
    return()
})
names(intra_module) <- signModule

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
