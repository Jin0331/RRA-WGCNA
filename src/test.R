# install.packages("reticulate")
library(reticulate)

# load py-function
use_condaenv(condaenv = "geo-py")
source_python("src/py-function.py")

# intra-extra analysis for top hub gene
# top_hub_gene <- extra_analysis_hub %>% pull(preferredName_A) %>% unique()
top_hub_gene <- intra_analysis_hub %>% unique()

# data-preprocessing for X, Y making
Y_col_name <- "sample_type"
DF <- data_trait %>% 
  select(starts_with(Y_col_name)) %>% 
  rownames_to_column(var = "sample") %>% 
  inner_join(x = ., y = robustdeg_ge %>% rownames_to_column(var = "sample"),
             by = "sample") %>% 
  column_to_rownames(var = "sample") %>% 
  select(1, all_of(top_hub_gene))

DF %>% rownames_to_column(var = "sample_barcode") %>% write_delim(file = "test.txt", delim = "\t", )

y_df <- DF %>% select_at(1)
x_df <- DF %>% select_at(-1)

lasso_coef <- feature_selection_LASSO(x_df, y_df)
svm_rfe_binary <- feature_selection_svm_rfecv(x_df, y_df)

lasso_validation_hub_gene <- x_df[ ,which(lasso_coef > 0)] %>% colnames()

svm_rfe_validation_hub_gene <- x_df[ ,which(svm_rfe_binary == TRUE)] %>% colnames()
intersect_gene_selected <- intersect(lasso_validation_hub_gene, svm_rfe_validation_hub_gene)

# intersection DF
DF %>% rownames_to_column(var = "sample_barcode") %>% 
  select(1:2, all_of(intersect_gene_selected)) %>% write_delim(file = "test.txt", delim = "\t")


# intersect venn
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

ggVennDiagram(x = list(LASS = lasso_validation_hub_gene, `SVM-RFE` = svm_rfe_validation_hub_gene)) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")


# DGIdb ----


for(index in 1:length(id_chunk)){
  parameters <- list(method="db2db", 
                     inputValues=id_chunk[[index]],
                     input="genbanknucleotideaccession",
                     outputs=c("genesymbol"),
                     taxonId="9606",
                     format="row")
  
  mapping_result[[index]] <- rcurl_request(json_url, parameters)
}




dgidb_interaction <- function(gene_name){
  base_url <- "https://dgidb.org"
  request_url <- paste0(base_url, "/api/v2/interactions.json?")
  result_list <- list()
  
  # chunk id
  id_chunk <- split(gene_name, ceiling(seq_along(top_hub_gene)/100))
  
  for(index in 1:length(id_chunk)){
    print(index)
    payload <-  list(genes = paste0(id_chunk[[index]], collapse = ","),
                     fda_approved_drug="true")
    
    # output
    dgidb_result <- POST(request_url, body = payload, encode = "form") %>%  
      httr::content(encoding = "UTF-8") 
  
    result_list[[index]] <- lapply(X = dgidb_result$matchedTerms, FUN = function(dgidb_element){
      gene_category <- dgidb_element$geneCategories %>% 
        sapply(X = ., FUN = function(value) {value$name}) %>% 
        paste0(collapse = ",")
      
      interaction <- dgidb_element$interactions %>% 
        sapply(X = ., FUN = function(value){
          drug_name <- value$drugName
          score <- value$score
          types <- value$interactionTypes %>% unlist() %>% paste0(collapse = "&")
          
          paste0(c(drug_name, score, types), collapse = ";") %>% 
            as_tibble() %>% 
            return()
          
          # return(drug_name)  
        }) %>% unlist() %>% 
        paste0(., collapse = "&")
      
      tibble(
        GENE_NAME = dgidb_element$geneName,
        DGI_GENE_CATEGORY = gene_category, 
        DGI_COUNT = length(dgidb_element$interactions),
        `DGI(DRUG_NAME;SCORE;TYPE)` = interaction
      )  %>% return()
      
    }) %>% bind_rows() 
  }
  
  result_list %>% bind_rows() %>% return()
}

temp <- dgidb_interaction(gene_name = top_hub_gene)




