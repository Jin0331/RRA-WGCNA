# library ----
suppressMessages({
  library(GEOquery)
  library(limma)
  library(RobustRankAggreg)
  library(pheatmap)
  library(DESeq2)
  library(impute)
  library(TCGAbiolinks)
  library(reticulate)
  library(ggVennDiagram)
  library(tidyverse)
  library(WGCNA)
  library(CorLevelPlot)
  library(caret)
  
  use_condaenv(condaenv = "geo-py")
  source_python("src/py-function.py")
})


# function ----
#' Function that returns numeric values with 2 decimal numbers.
#' 
#' @param x input numeric value with N decimal numbers.
#' @return a numeric value with 2 decimal numbers.
#' @examples
#' aaa <- dec_two(8.31232)
dec_two <- function(x) {
  return (format(round(x, 2), nsmall = 2));
}

# mapping function ====
symbol_mapping <- function(ge, col_name, platform_ann_df){
  thisGene <- rownames(ge) %>%
    lapply(FUN = function(value){
      platform_ann_df %>% 
        dplyr::filter(ID == value) %>% 
        pull(col_name)
    }) %>% 
    do.call(c, .)
  
  if(col_name != "gene_assignment"){
    
    return(thisGene)
    
  } else { # gene assignment
    
    SIZE_SPLIT_STRING <- 2
    GENE_SYMBOL_INDEX <- 2
    
    thisGene %>% 
      lapply(X = ., FUN = function(value){
        split_string <- strsplit(value, "//")
        if (length(split_string[[1]]) >= SIZE_SPLIT_STRING) {
          thisGeneSymbol_temp <- NULL
          thisGeneSymbol  <- NULL
          thisGeneSymbol_temp <- split_string[[1]][GENE_SYMBOL_INDEX]
          thisGeneSymbol <- gsub(" ", "", thisGeneSymbol_temp, fixed = TRUE)
          thisGeneSymbol %>% return()
        } else {
          value %>% return()
        }
      }) %>% do.call(c, .) %>% 
      return()
  }
}
biodbnet_db2db <- function(id){
  base_url <- "https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/"
  json_url <- paste0(base_url, "biodbnetRestApi.json?")
  mapping_result <- list()
  id_chunk <- split(id, ceiling(seq_along(id)/1000))
  
  for(index in 1:length(id_chunk)){
    parameters <- list(method="db2db", 
                       inputValues=id_chunk[[index]],
                       input="genbanknucleotideaccession",
                       outputs=c("genesymbol"),
                       taxonId="9606",
                       format="row")
    
    mapping_result[[index]] <- rcurl_request(json_url, parameters)
  }
  
  mapping_result <- mapping_result %>% 
    bind_rows() %>% 
    filter(`Gene Symbol` != "-")
  
  return(mapping_result)
}


# GEO download function  ====
#' Function that reads in a URL to check and verifies if it exists (function taken from https://stackoverflow.com/a/12195574 )
#' 
#' @param url the URL of a webpage
#' @return the output of a webpage verification check
#' @examples y <- readUrl("http://stat.ethz.ch/R-manual/R-devel/library/base/html/connections.html")
readUrl <- function(url) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      
      
      readLines(con=url, warn=FALSE) 
      # The return value of `readLines()` is the actual value 
      # that will be returned in case there is no condition 
      # (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", url))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return("EMPTY_STRING")
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", url))
      message(cond)
      # Choose a return value in case of warning
      return("EMPTY_STRING")
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      message(paste("Processed URL:", url))
    }
  )    
  return(out)
}
rcurl_request <- function(service_url, parameters) {
  # Collapse all parameters into one string
  all_parameters <- paste(
    sapply(names(parameters), 
           FUN=function(param_name, parameters) {
             paste(param_name, paste(parameters[[param_name]], collapse=','), collapse='', sep='=')
           }, 
           parameters),
    collapse="&")
  
  # Paste base URL and parameters
  requested_url <- paste0(service_url, all_parameters)
  
  # Encode URL (in case there would be any space character for instance)
  # requested_url <- URLencode(requested_url)
  
  # Start request to service
  response <- GET(requested_url)
  
  raise <- content(response, as="text")
  #parse JSON
  new <- fromJSON(raise)
  
  return(new)
}
#' Function that reads in the GEO code of a dataset, and returns the gene expression dataframe.
#' 
#' @param datasetGeoCode the GEO code of a dataset.
#' @param retrieveGeneSymbols a boolean flag stating if the function should retrieve the gene symbols or not.
#' @param verbose a boolean flag stating if helping messages should be printed or not
#' @return a gene expression dataset.
#' @examples
#' geneExpressionDF1 <- getGeneExpressionFromGEO("GSE3268", FALSE, FALSE)
getGeneExpressionFromGEO <- function(datasetGeoCode, retrieveGeneSymbols, verbose = FALSE) {
  
  GSE_code <- datasetGeoCode
  
  
  # r <- NULL
  # attempt <- 1
  # while( is.null(r) && attempt <= 3 ) {
  #   attempt <- attempt + 1
  #   try(
  #     r <- some_function_that_may_fail()
  #   )
  # } 
  
  # check   URL
  checked_html_text <- "EMPTY_STRING"
  checked_html_text <- retry(expr =  xml2::read_html("https://ftp.ncbi.nlm.nih.gov/geo/series/"),
                             maxErrors=10, 
                             sleep=2)
  
  checked_html_text_url <- "EMPTY_STRING"
  url_to_check <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", datasetGeoCode)
  GSE_code_for_url <- GSE_code
  GSE_code_for_url <- substr(GSE_code_for_url,1,nchar(GSE_code_for_url)-3)
  GSE_code_for_url <- paste0(GSE_code_for_url, "nnn")
  complete_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", GSE_code_for_url, "/", GSE_code)
  
  checked_html_text_url <- lapply(complete_url, readUrl)
  
  
  
  #             
  if(all(checked_html_text == "EMPTY_STRING")) {
    
    cat("The web url https://ftp.ncbi.nlm.nih.gov/geo/series/ is unavailable right now. Please try again later. The function will stop here\n")
    return(NULL)
    
  } else if(all(checked_html_text_url == "EMPTY_STRING" | is.null(checked_html_text_url[[1]]) )) {
    
    cat("The web url ", complete_url," is unavailable right now (Error 404 webpage not found). The GEO code might be wrong. The function will stop here\n", sep="")
    return(NULL)        
    
  } else {
    
    gset <- retry(expr = GEOquery::getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE), maxErrors = 10, sleep = 2)
    
    thisGEOplatform <- toString((gset)[[1]]@annotation)
    
    if (length(gset) > 1) 
      idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    gene_expression <- as.data.frame(Biobase::exprs(gset))
    
    
    if(retrieveGeneSymbols == TRUE) {
      gene_expression$GeneSymbol <- ""
      
      
      # we retrieve the platform details
      platform_ann <- retry(expr = annotate::readGEOAnn(GEOAccNum = thisGEOplatform), maxErrors = 10, sleep = 2)
      platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
      
      if (verbose == TRUE) {
        print("sort(names(platform_ann_df))")
        print(sort(names(platform_ann_df)))
      }
      
      # "gene_assignment
      platformsWithGene_assignmentField <- c("GPL11532", "GPL23126", "GPL6244", "GPL17586", "GPL5175")
      
      # "Gene Symbol"
      platformsWithGeneSpaceSymbolField <- c("GPL80", "GPL8300", "GPL80", "GPL96", "GPL570", "GPL571", "GPL3921")
      
      # "gene_symbol"
      platformsWithGene_SymbolField <- c("GPL20115")
      
      # "HUGOname"
      platformsWithHUGOnamelField <- c("GPL5918")
      
      # "symbol"
      platformsWithSymbolField <- c("GPL1293", "GPL6102", "GPL6104", "GPL6883", "GPL6884", "GPL10558")
      
      # "GENE_SYMBOL
      platformsWith_GENE_SYMBOL_Field <- c("GPL13497", "GPL14550", "GPL17077", "GPL6480")
      
      # if symbol
      if(thisGEOplatform %in% c(platformsWithHUGOnamelField, platformsWithGene_assignmentField, platformsWithGeneSpaceSymbolField, platformsWithGene_SymbolField, platformsWithSymbolField, platformsWith_GENE_SYMBOL_Field)   ) {
        
        emptyGeneSymbol <- ""
        FIRST_GENE_EXPRESSION_INDEX <- 2
        
        if (verbose == TRUE)    
          cat("\n[start] loop for the association of the gene symbols to the probeset ID's\n", sep="")
        
        if(thisGEOplatform %in% platformsWithGeneSpaceSymbolField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "Gene Symbol", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWithGene_SymbolField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "GeneSymbol", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWithSymbolField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "Symbol", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWith_GENE_SYMBOL_Field)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "GENE_SYMBOL", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWithGene_assignmentField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "gene_assignment", platform_ann_df = platform_ann_df)
        if(thisGEOplatform %in% platformsWithHUGOnamelField)
          thisSymbol <- symbol_mapping(ge = gene_expression, col_name = "HUGOname", platform_ann_df = platform_ann_df)
        gene_expression$GeneSymbol <- thisSymbol
        
        if (verbose == TRUE) 
          cat("\n [end] loop for the association of the gene symbols to the probeset ID's \n ", sep="")
        
      }  else{
        if (verbose == TRUE) 
          cat("\n\n[Impossible to perform gene symbol retrieval]\n")
        if (verbose == TRUE) 
          cat("We're sorry but the indicated platform (", thisGEOplatform, ") is not among the platforms included in this function.\nThe gene symbol retrieval cannot be performed.\nPlease visit the https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", thisGEOplatform, " website for more information about this platform.\n\n", sep="")
        gene_expression$GeneSymbol <- NULL
      }
      
    }
    
    # gene expression duplicated probe remove using mdedin max
    # duplicate remove probe / gene
    gene_name <- gene_expression$GeneSymbol
    sample_name <- colnames(gene_expression)[1:(length(colnames(gene_expression)) - 1)]
    probe_MAD <- gene_expression %>% 
      select(-GeneSymbol) %>% 
      as.matrix() %>% 
      apply(.,1,mad) %>% 
      tibble(id = seq(length(.)), gene_name = gene_name, MAD = .) %>% 
      arrange(desc(MAD))
    probe_MAD_dup <- probe_MAD[which(!duplicated(probe_MAD$gene_name)),] %>% 
      filter(gene_name != "") %>% 
      arrange(id)
    gene_expression <- gene_expression[probe_MAD_dup$id, ]
    
    geneExpression_dup <- gene_expression %>% 
      select(-GeneSymbol) %>% 
      as.matrix()
    rownames(geneExpression_dup) <- gene_expression$GeneSymbol
    colnames(geneExpression_dup) <- colnames(gene_expression)[1:(length(gene_expression) - 1)]
    
    return(list(gene_expression = geneExpression_dup , 
                pheno = phenoData(gset)))
    
  }   }

# DEA & RRA function ====
run_limma <- function(ge, de){
  fit <- lmFit(ge,de)
  cont <- makeContrasts(TP-NT,levels=de) # 데이터에 맞추어 manual로 설정해야 됨
  fit.cont <- contrasts.fit(fit,cont)
  fit.cont <- eBayes(fit.cont)
  res <- topTable(fit.cont,number=Inf) 
  target <- res[res$P.Value < 0.05,] %>% 
    rownames_to_column() %>% as_tibble()
  
  return(target)
}
rra_extract <- function(ml, logfc = 0.0, fdr = 0.05){
  # combine deg
  combine_degs <- names(ml) %>% 
    lapply(X = ., FUN = function(list_name){
      tmp <- ml[[list_name]] %>% 
        dplyr::filter(adj.P.Val < fdr & (logFC > logfc | logFC < -(logfc))) %>% 
        arrange(desc(logFC)) %>%
        dplyr::select(rowname, logFC)
      colnames(tmp) <- c("GENE", list_name)
      return(tmp)
    }) %>% 
    purrr::reduce(., left_join, by = "GENE") %>% 
    bind_cols(., apply(.[,-1], 1, mean, na.rm = TRUE) %>% 
                tibble(group = .)) 
  
  # up-regulated
  updown_degs <- names(ml) %>% 
    lapply(X = ., FUN = function(list_name){
      ml[[list_name]] %>% 
        dplyr::filter(adj.P.Val < fdr & (logFC > logfc | logFC < -(logfc))) %>% 
        arrange(adj.P.Val) %>% 
        dplyr::pull(rowname) %>%
        return()
    }) 
  
  # Aggregate the inputs
  # run RRA
  updown_deg_rra <- aggregateRanks(glist = updown_degs, method = "RRA") %>%
    as_tibble() %>% 
    dplyr::filter(Score < 0.05) %>% 
    arrange(Score)
  
  # # 1 - combine deg, 2 - up_down-regulated RRA
  list(combine_degs = combine_degs, updown_rra = updown_deg_rra) %>% return()
}
rra_analysis <- function(m_list, logfc = 0, fdr = 0.05, save_path =  "RData/GEO_RobustDEGs_norm.RData"){
  rra_result <- rra_extract(ml = m_list, logfc = logfc, fdr = fdr)
  combine_degs_rra <- rra_result[[1]] %>% 
    dplyr::filter(GENE %in% rra_result[[2]]$Name) %>% 
    arrange(desc(group))
  
  up_down_rra_gene <- bind_rows(head(combine_degs_rra, 20),
                                tail(combine_degs_rra, 20)) %>% 
    dplyr::select(-group)
  
  # heatmap
  combine_degs_m <- up_down_rra_gene[,-1] %>% as.matrix()
  rownames(combine_degs_m) <- up_down_rra_gene$GENE
  colnames(combine_degs_m) <- colnames(up_down_rra_gene)[2:length(colnames(up_down_rra_gene))]
  pheatmap(combine_degs_m, 
           display_numbers = TRUE,
           number_color = "black",
           fontsize_number = 10,
           border_color = "black",
           cluster_rows = F,
           cluster_cols = F,
           cellwidth = 35,
           cellheight = 10
  )
  # save
  save(up_down_rra_gene, file = save_path)
  
  return(combine_degs_rra %>% 
           dplyr::pull(1))
}

# WGCNA function ==== 
mergecutheight_test <- function(pr_name, robustdegs){
  network_list <- list()
  max_module_cnt <- 0
  max_module <- NULL
  
  for(mch_value in seq(from = 0.1, to = 0.4, by = 0.02)){
    print(paste0("mch_value : ", mch_value))
    network <- network_preprocessing(pr_name = "LIHC", robustdegs = robustdegs, mch = mch_value, time_stamp = "test")
    module_cnt <- length(network[[2]]$colors %>% unique()) - 1
    
    if(max_module_cnt <= module_cnt){
      max_module <- mch_value
      max_module_cnt <- module_cnt
    }
    # network_list[[as.character(mch_value)]] <- list(network, module_cnt)
  }
  
  return(list(max_mch = max_module, max_module_cnt = max_module_cnt))
}
network_preprocessing <- function(pr_name, robustdegs, mch = 0.25, time_stamp){
  allowWGCNAThreads(nThreads = 15)
  
  suppressMessages({
    cancer_fpkm <- load_tcga_dataset(pkl_path = "PKL/", raw_path = "RAW_DATA/", cancer_type = pr_name) %>% 
      rownames_to_column(var = "id")
    
    gene_probe_mapping <- read_delim(file = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap",
                                     delim = "\t", show_col_types = FALSE, progress = FALSE) %>% 
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
    sft <- pickSoftThreshold(robustdeg_ge, powerVector = powers, verbose = 0) %>% 
      .$powerEstimate
    
    # network construct
    # net construct
    # time_stamp <- Sys.time()
    log_save <- paste("WGCNA_LOG", pr_name, time_stamp, sep = "/")
    dir.create(paste("WGCNA_LOG", pr_name, time_stamp, sep = "/"), showWarnings = FALSE, recursive = TRUE)
    net <- blockwiseModules(datExpr = robustdeg_ge, 
                            power = sft,
                            corType = "pearson",
                            TOMType = "unsigned", 
                            minModuleSize = 50,
                            # maxBlockSize = 500,
                            reassignThreshold = 0, 
                            mergeCutHeight = mch,
                            numericLabels = TRUE, 
                            pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            saveTOMFileBase = paste0(log_save, "/TCGA-", pr_name),
                            verbose = 0)
    
  })
    
  
  return(list(deg = robustdeg_ge, network = net))
}
find_key_modulegene <- function(pr_name, network, MEs, select_clinical=NULL, mm=0.85, gs=0.25, time_stamp, binarytocategory=FALSE){
  
  # variable
  expression_sample <- rownames(network[[1]])
  nGenes <- ncol(network[[1]])
  nSamples <- nrow(network[[1]])
  log_save <- paste("WGCNA_LOG", pr_name, time_stamp, sep = "/")
  dir.create(log_save, recursive = T, showWarnings = FALSE)
  
  # UCSCXena clinical
  clinical_trait <- read_delim("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FLIHC_clinicalMatrix",
                               delim = "\t", show_col_types = FALSE, progress = FALSE)
  survival_trait <- read_delim("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FLIHC_survival.txt",
                               delim = "\t", show_col_types = FALSE, progress = FALSE) %>% 
    dplyr::select(sample, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time)
  
  immune_trait <- read_delim("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Subtype_Immune_Model_Based.txt.gz",
                             delim = "\t", show_col_types = FALSE, progress = FALSE)
  
  molecular_trait <- read_delim("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGASubtype.20170308.tsv.gz",
                                delim = "\t", show_col_types = FALSE, progress = FALSE) %>% 
    select(sampleID, Subtype_Integrative)
  
  clinical_trait <- left_join(x = clinical_trait, y = survival_trait, by = c("sampleID" = "sample")) %>% 
    left_join(x = ., y = immune_trait, by = c("sampleID" = "sample")) %>% 
    left_join(x = ., y = molecular_trait, by = "sampleID")
  
  default_clinical <- c('sample_type', 'pathologic_stage', 'Subtype_Immune_Model_Based',
                        'Subtype_Integrative')
  use_clinical <- c(default_clinical, select_clinical)
  
  # clinical trait preprocessing
  traitRows <- match(expression_sample, clinical_trait$sampleID)
  data_trait <- clinical_trait[traitRows, ] %>% 
    column_to_rownames(var = "sampleID") %>% 
    dplyr::select(all_of(use_clinical)) %>% 
    mutate(sample_type = ifelse(sample_type == "Primary Tumor", 1, 0),
           pathologic_stage = str_remove_all(string = pathologic_stage, pattern = "A|B|C"),
           # pathologic_T = str_remove_all(string = pathologic_T, pattern = "a|b")
           ) %>%  # sample type 한정, pathogenic stage
    replace(is.na(.), 0) %>% 
    mutate_all(as.factor) %>% 
    mutate_all(as.numeric) %>% 
    mutate_all(function(value){value - 1})

  # binary to category
  if(binarytocategory){
    # category to binary
    dummy <- dummyVars(" ~ .", data=data_trait)
    data_trait <- data.frame(predict(dummy, newdata=data_trait))  
  }
  
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
  GS <- lapply(X =  data_trait %>% colnames(), FUN = function(s_type){
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
    
    module_MM_GS_filtered <- lapply(X = data_trait %>% colnames(), FUN = function(t){
      GS %>% 
        select(gene, GS = t) %>% 
        filter(gene %in% module_MM_gene, abs(GS) > gs)})
    names(module_MM_GS_filtered) <- data_trait %>% colnames()
    
    return(module_MM_GS_filtered)
  })
  
  total_keyhub <- list()
  for(index in data_trait %>% colnames()){
    tmp <- lapply(X = intra_module, FUN = function(trait){
      trait[[index]]
    }) %>% bind_rows() %>% dplyr::pull(gene)
    
    if(length(tmp) != 0)
      total_keyhub[[index]] <- tmp
  }
  
  total_keyhub_merge <- names(total_keyhub) %>% lapply(X = ., FUN = function(lname){
    col_name <- lname %>% str_split(pattern = "\\.") %>% unlist %>% .[1]
    
    tibble(gene = total_keyhub[[lname]], col_name) %>% return()
  }) %>% bind_rows() %>% 
    split(x = ., f = .$col_name) %>% 
    lapply(X = ., FUN = function(df){df %>% dplyr::pull(gene)})
  
  
  # plot save
  sample_cluster_plot(network = network[[1]], clinical_trait = data_trait, save_path = log_save)
  module_cluster_plot(network = network, save_path = log_save)
  module_trait_plot(moduleTraitCor = moduleTraitCor, moduleTraitPvalue = moduleTraitPvalue,
                    data_trait = data_trait, MEs = MEs, save_path = log_save)
  gene_module_size_plot(gene_module_key_groph = gene_module_key %>% group_by(module) %>% 
                          summarise(module_size = n()),
                        save_path = log_save)
  key_hub_intersection_plot(total_keyhub = total_keyhub_merge, save_path = log_save)

  total_keyhub <- total_keyhub_merge %>% unlist() %>% unname() %>% unique()
  
  return(list(total_keyhub = total_keyhub, clinical_trait = data_trait, total_keyhub_merge = total_keyhub_merge, network = network))
  
  
}

sample_cluster_plot <- function(network, clinical_trait, save_path){
  sampleTree <-  hclust(dist(network[[1]]), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  # plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
  #      cex.axis = 1.5, cex.main = 2)
  
  traitColors <- numbers2colors(clinical_trait, signed = FALSE);
  
  png(paste0(save_path, "/sample_cluster.png"), width = 1000, height = 1000)
  plotDendroAndColors(sampleTree, traitColors,
                      groupLabels = names(clinical_trait),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
  
  
}
module_cluster_plot <- function(network, save_path){
  sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors <- labels2colors(network[[2]]$colors)
  # Plot the dendrogram and the module colors underneath
  
  png(paste0(save_path, "/module_cluster.png"), width = 1000, height = 1000)
  plotDendroAndColors(network[[2]]$dendrograms[[1]], mergedColors[network[[2]]$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)  
  dev.off()
}
module_trait_plot <- function(moduleTraitCor, moduleTraitPvalue, data_trait, MEs, save_path){
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  png(paste0(save_path, "/module_trait.png"), width = 1000, height = 1000)
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
  dev.off()
}
gene_module_size_plot <- function(gene_module_key_groph, save_path){
  p <- ggpubr::ggbarplot(gene_module_key_groph, x = "module", y = "module_size",
                         color = "module",
                         fill = "module",
                         sort.val = "asc",
                         sort.by.groups = FALSE,
                         x.text.angle = 90,          # Rotate vertically x axis texts
                         title = "Top3 Module Size",
                         ylab = "Module Size",
                         legend.title = "Module",
                         rotate = TRUE,
                         ggtheme = theme_minimal())
  ggsave(p, filename = paste0(save_path, "/top3_modulesize.png"))
  return(p)
}
key_hub_intersection_plot <- function(total_keyhub, save_path){
  
  p <- ggVennDiagram::ggVennDiagram(x = total_keyhub, label_size = 7) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none")
  
  ggsave(p, filename = paste0(save_path, "/keyhubgene_intersection.png"), dpi = 200, width = 30, height = 10)
  return(p)
}

# # plot - GS & MM correlation 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "GS for sample type",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# Machine learning function ====
# STRING function ====
retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=5, sleep=0) {
  attempts = 0
  retval = try(eval(expr))
  while (isError(retval)) {
    attempts = attempts + 1
    if (attempts >= maxErrors) {
      msg = sprintf("retry: too many retries [[%s]]", capture.output(str(retval)))
      flog.fatal(msg)
      stop(msg)
    } else {
      msg = sprintf("retry: error in attempt %i/%i [[%s]]", attempts, maxErrors, 
                    capture.output(str(retval)))
      flog.error(msg)
      warning(msg)
    }
    if (sleep > 0) Sys.sleep(sleep)
    retval = try(eval(expr))
  }
  return(retval)
}
#' Function that returns STRING network
#' @param hub_gene input character vector
#' @return network dataframe
#' @examples
#' aaa <- string_network(hub_gene = my_gene)
string_network <- function(hub_gene){
  # URLs
  string_api_url <- "https://version-11-5.string-db.org/api"
  output_format <- "tsv-no-header"
  method <- "network"
  request_url <- paste(string_api_url, output_format, method, sep = "/")
  
  # post payload
  params <- list(
    identifiers = paste0(hub_gene, collapse = "%0d"),
    species = "9606",
    caller_identity = "www.hallym.ac.kr"
  )
  
  # output
  network_colname <- c('stringId_A','stringId_B','preferredName_A','preferredName_B','ncbiTaxonId',
                       'score','nscore','fscore','pscore','ascore','escore','dscore','tscore')
  
  ppi_network <- POST(request_url, body = params, encode = "form") %>%  
    httr::content(encoding = "UTF-8") 
  colnames(ppi_network) <- network_colname
  ppi_network <- ppi_network %>% arrange(preferredName_A) %>% distinct_all()
  
  return(ppi_network)
}
