# function
source("src/r-function.R")

# LOCAL VARIABLE
# rdata_path <- "RData/"
load(file = "RData/HCC_GEO_RobustDEGs_norm.RData")
pr_name <- "LIHC"

# weighted gene co-expression network costruction
network <- network_preprocessing(pr_name = "LIHC", robustdegs = robustdegs)

# network variable
moduleLabels <- network[[2]]$colors
moduleColors <-  labels2colors(network[[2]]$colors)
MEs <- network[[2]]$MEs;
geneTree <- network[[2]]$dendrograms[[1]]

# Calculate ME
MEs0 <- moduleEigengenes(robustdeg_ge, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
    

find_key_modulegene <- function(network, MEs, select_clinical=NULL,){
  
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
  
  geneTraitSignificance <- as.data.frame(cor(robustdeg_ge, s_type, use = "p"))
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
  # plot - sample cluster
  {
    sampleTree <-  hclust(dist(robustdeg_ge), method = "average");
    # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
    # The user should change the dimensions if the window is too large or too small.
    sizeGrWindow(12,9)
    #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
    par(cex = 0.6);
    par(mar = c(0,4,2,0))
    plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2)
  }
  # plot - module cluster
  {
    sizeGrWindow(12, 9)
    # Convert labels to colors for plotting
    mergedColors <- labels2colors(net$colors)
    # Plot the dendrogram and the module colors underneath
    plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  }
  # plot - module-trait relationships
  {
    sizeGrWindow(10,6)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3));
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = names(data_trait),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = greenWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
  }
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