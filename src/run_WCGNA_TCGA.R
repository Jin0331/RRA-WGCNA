library(DESeq2)
library(impute)
library(WGCNA)

rdata_path <- "RData/"
load(file = "RData/HCC_GEO_RobustDEGs_norm.RData")
robustdegs <- up_down_rra_gene %>% pull(1)

run_deseq_normal <- function(pr_name, rdata_path, deg_path, batch_removal){
  register(MulticoreParam(20))
  suppressMessages({
    if((!file.exists(paste0(rdata_path, pr_name, "_normal.RData"))) | 
       (!file.exists(paste0(rdata_path, pr_name, "_RnaseqSE_normal.RData")))){
      query <- GDCquery(project = paste0("TCGA-", pr_name), 
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        experimental.strategy = "RNA-Seq",
                        platform = "Illumina HiSeq",
                        file.type = "results",
                        sample.type = c("Primary Tumor", "Solid Tissue Normal"), 
                        legacy = TRUE)
      
      GDCdownload(query)
      RnaseqSE <- GDCprepare(query)
      
      save(RnaseqSE, file = paste0(rdata_path, pr_name, "_RnaseqSE_normal.RData"))
      
      Rnaseq_CorOutliers <- assay(RnaseqSE, "raw_count") # to matrix
      
      # normalization of genes, # quantile filter of genes
      dataNorm <- TCGAanalyze_Normalization(tabDF = Rnaseq_CorOutliers, geneInfo =  geneInfo)
      dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                        method = "quantile", 
                                        qnt.cut =  0.25)
      
      dataFilt <- varianceStabilizingTransformation(dataFilt)
      
      save(dataFilt, file = paste0(rdata_path, pr_name, "_normal.RData"))
    } else {
      load(paste0(rdata_path, pr_name, "_RnaseqSE_normal.RData"))
      load(paste0(rdata_path, pr_name, "_normal.RData"))
    }
    
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
    
    # plot
    sampleTree <-  hclust(dist(robustdeg_ge), method = "average");
    # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
    # The user should change the dimensions if the window is too large or too small.
    sizeGrWindow(12,9)
    #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
    par(cex = 0.6);
    par(mar = c(0,4,2,0))
    plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2)
    
    # power calculation
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))
    sft <- pickSoftThreshold(robustdeg_ge, powerVector = powers, verbose = 5)
    
    # net construct
    net <- blockwiseModules(datExpr = robustdeg_ge, 
                            power = 14,
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
    
    # plot
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
    
    moduleLabels <- net$colors
    moduleColors <-  labels2colors(net$colors)
    MEs <- net$MEs;
    geneTree <- net$dendrograms[[1]]
    
    # clinical trait
    # clinical_trait <- GDCquery_clinic(project = paste0("TCGA-", pr_name), type = "clinical")
    # expression_sample <- rownames(robustdeg_ge) %>% as_tibble() %>% 
    #   separate(col = value, into = c("A","B","C")) %>% 
    #   unite(col = id, sep = "-") %>% pull(1)
    clinical_trait <- read_delim("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FLIHC_clinicalMatrix",
                                 delim = "\t")
    write_delim(clinical_trait, file = "TCGA-LIHC_Clinical_TCGAibolinks.txt", delim = "\t")
    
    expression_sample <- rownames(robustdeg_ge) %>% as_tibble() %>% 
      separate(col = value, into = c("A","B","C","D")) %>% 
      unite(col = id, sep = "-") %>% pull(1) %>% 
      gsub('.{1}$', '', .)
    
    traitRows <- match(expression_sample, clinical_trait$sampleID)
    data_trait <- clinical_trait[traitRows, ]
    data_trait %>% select(sample_type)
    
    rownames(data_trait) <- clinical_trait[traitRows, 3]
    
    expression_sample <- rownames(robustdeg_ge)
    traitRows = match(expression_sample, clinical_trait$Affy_GSM)
    data_trait <- clinical_trait[traitRows, -3]
    rownames(data_trait) <- clinical_trait[traitRows, 3]
  })  
  
  return(tcga_deseq_result_tidy)
}