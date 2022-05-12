# function
source("src/r-function.R")

use_condaenv(condaenv = "geo-py")
source_python("src/py-function.py")

# LOCAL VARIABLE
# rdata_path <- "RData/"
load(file = "RData/HCC_GEO_RobustDEGs_norm.RData")
robustdegs <- robust_degs
pr_name <- "LIHC"


# weighted gene co-expression network costruction

## GDC input data prerpocessing

network_preprocessing <- function(pr_name, robustdegs){
  cancer_fpkm <- load_tcga_dataset(pkl_path = "PKL/", raw_path = "RAW_DATA/", cancer_type = pr_name) %>% 
    rownames_to_column(var = "id")
  
  gene_probe_mapping <- read_delim(file = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap",
                                   delim = "\t") %>% 
    select(id, gene) %>% 
    inner_join(x = ., y = cancer_fpkm, by = "id") %>% 
    select(-id)
  
  dataFilt <- gene_probe_mapping %>% 
    mutate_if(is.numeric, .funs = function(value){
      2^value - 0.001 %>% return()
    }) %>% 
    mutate_if(is.numeric, .funs = function(value){
      ifelse(value < 0, 0, value) %>% return()
    }) %>% 
    mutate_if(is.numeric, .funs = function(value){
      log2(value + 1) %>% return()
    }) %>% 
    distinct(gene, .keep_all = TRUE) %>% 
    column_to_rownames(var = "gene") %>% 
    as.matrix()
  
  
  # row to col AND DESeq2 normalized log2(x+1)
  robustdeg_ge <- lapply(X = robustdegs, FUN = function(deg){
    error <- FALSE
    tryCatch(
      expr = {
        tmp <- dataFilt[deg, ]
      },
      error = function(e) {
        error <<- TRUE
      }
    )
    if(error){
      return(NULL)
    } else {
      tmp <- as.matrix(tmp) %>% t()
      rownames(tmp) <- deg
      return(tmp)
    }}) %>% do.call(rbind, .) %>% 
    t() %>% 
    as.data.frame()
  
  # split rowname
  rownames(robustdeg_ge) <- rownames(robustdeg_ge) %>% as_tibble() %>% 
    separate(col = value, into = c("A","B","C","D")) %>% 
    unite(col = id, sep = "-") %>% dplyr::pull(1)
  # gsub('.{1}$', '', .)
  
  return(robustdeg_ge)
}

robustdeg_ge <- network_preprocessing(pr_name = "LIHC")


    
    
    
    # sample & gene filtering
    gsg <- goodSamplesGenes(robustdeg_ge, verbose = 3)
    if (!gsg$allOK){
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(robustdeg_ge)[!gsg$goodGenes], collapse = ", ")));
      if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(robustdeg_ge)[!gsg$goodSamples], collapse = ", ")));
      # Remove the offending genes and samples from the data:
      robustdeg_ge <-  robustdeg_ge[gsg$goodSamples, gsg$goodGenes]
    }
    


    # power calculation
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))
    sft <- pickSoftThreshold(robustdeg_ge, powerVector = powers, verbose = 5)
    
    # net construct
    net <- blockwiseModules(datExpr = robustdeg_ge, 
                            power = 9,
                            corType = "pearson",
                            TOMType = "unsigned", 
                            minModuleSize = 30,
                            reassignThreshold = 0, 
                            mergeCutHeight = 0.20,
                            numericLabels = TRUE, 
                            pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            saveTOMFileBase = "TCGA-LIHC",
                            verbose = 3)
    
    moduleLabels <- net$colors
    moduleColors <-  labels2colors(net$colors)
    MEs <- net$MEs;
    geneTree <- net$dendrograms[[1]]
    
    # clinical trait
    # TCGAbiolinks clinical data
    # clinical_trait <- GDCquery_clinic(project = paste0("TCGA-", pr_name), type = "clinical")
    
    # UCSCXena clinical
    clinical_trait <- read_delim("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FLIHC_clinicalMatrix",
                                 delim = "\t")
    survival_trait <- read_delim("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FLIHC_survival.txt",
                                 delim = "\t") %>% 
      dplyr::select(sample, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time)
    
    clinical_trait <- left_join(x = clinical_trait, y = survival_trait, by = c("sampleID" = "sample"))
    # write_delim(clinical_trait, file = "TCGA-LIHC_Clinical_TCGAibolinks.txt", delim = "\t")
    
    expression_sample <- rownames(robustdeg_ge)
    
    # Factor 
    # Sample Type
    ## group 0 - long survival group (wild type)
    ## group 1 - short survival group (mu type)
    traitRows <- match(expression_sample, clinical_trait$sampleID)
    data_trait <- clinical_trait[traitRows, ] %>% 
      column_to_rownames(var = "sampleID") %>% 
      dplyr::select(sample_type, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time,age_at_initial_pathologic_diagnosis, 
             pathologic_T, pathologic_M, pathologic_N, pathologic_stage, child_pugh_classification_grade, 
             fibrosis_ishak_score) %>% 
      mutate(sample_type = ifelse(sample_type == "Primary Tumor", 1, 0)) %>%  # sample type 한정
      mutate_if(is.character, as.factor) %>% 
      mutate_all(as.numeric)
    data_trait[is.na(data_trait)] <- 0
    # write_delim(data_trait, file = "TCGA-LIHC_Clinical_impute0.txt", delim = "\t")
    
    # # dummy variable (binary trait)
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
    
    
    
    # relating modules to external clinical traits ----
    nGenes <- ncol(robustdeg_ge)
    nSamples <- nrow(robustdeg_ge)
    MEs0 <- moduleEigengenes(robustdeg_ge, moduleColors)$eigengenes
    MEs <- orderMEs(MEs0)
    
    # trait cor
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
    }) %>% bind_rows() %>% arrange(desc(signTrait)) %>% 
      pull(1) %>% 
      .[1:3]
    
    
    ## Gene significance
    sig_trait <- "sample_type"
    s_type <- as.data.frame(data_trait$sample_type)
    names(s_type) <- sig_trait
    
    modNames <- substring(names(MEs), 3)
    geneModuleMembership <- as.data.frame(cor(robustdeg_ge, MEs, use = "p"))
    MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
    
    # module membership
    names(geneModuleMembership) <- paste("MM", modNames, sep="")
    names(MMPvalue) <- paste("p.MM", modNames, sep="")
    
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
  # })  
  
  # return(tcga_deseq_result_tidy)
# }