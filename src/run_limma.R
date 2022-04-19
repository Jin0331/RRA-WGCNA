# load library
library(GEOquery)
library(limma)
library(tidyverse)

# function
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
multiple_limma <- list()


# geoquery
gse_name <- "GSE60502"
gse <- getGEO(gse_name,GSEMatrix=TRUE)
gse <- gse[[1]]

eset <- exprs(gse)
fset <- fData(gse)

# if not in hugogene and in genbank
{
  # fset_sep <- fset %>% separate(col = GB_LIST, sep = ",", into = LETTERS[1:4]) %>% as_tibble()
  # 
  # A_mapping <- biodbnet_db2db(id = fset_sep$A)
  # B_mapping <- biodbnet_db2db(id = fset_sep$B)
  # C_mapping <- biodbnet_db2db(id = fset_sep$C)
  # D_mapping <- biodbnet_db2db(id = fset_sep$D)
  # 
  # 
  # fset_sep_map <- fset_sep %>% 
  #   left_join(x = ., y = A_mapping, by = c("A" = "InputValue")) %>% 
  #   left_join(x = ., y = B_mapping, by = c("B" = "InputValue")) %>% 
  #   left_join(x = ., y = C_mapping, by = c("C" = "InputValue")) %>% 
  #   left_join(x = ., y = D_mapping, by = c("D" = "InputValue")) %>% 
  #   unite(gene_symbol, c(`Gene Symbol.x`,`Gene Symbol.y`,
  #                        `Gene Symbol.x.x`,`Gene Symbol.y.y`), 
  #         sep = ";", remove=TRUE, na.rm = TRUE) %>% 
  #   unite(GB_LIST, c(A,B,C,D), sep = ",", remove = TRUE, na.rm = TRUE)

}

colnames(fset) #"Symbol" 또는 "Gene Symbol"가 있는지 확인
rownames(eset) <- fset[,"Gene Symbol"]

# sample selection 
pset <- phenoData(gse)
View(pset@data)

# tissue:ch1, source_name_ch1
sample_select <- pset@data %>% 
  # filter(str_detect(characteristics_ch1, "HCC")) %>% 
  pull(`source_name_ch1`)
grp <- pset@data %>% 
  # filter(str_detect(characteristics_ch1, "HCC")) %>% 
  pull(`source_name_ch1`) %>% 
  lapply(X = ., FUN = tolower) %>% 
  unlist() %>% 
  as.factor()

# eset <- eset[, sample_select]


# duplicate remove probe / gene
probe_MAD <- apply(eset,1,mad) %>% 
  tibble(id = seq(length(.)), gene_name = names(.), MAD = .) %>% 
  arrange(desc(MAD))
probe_MAD_dup <- probe_MAD[which(!duplicated(probe_MAD$gene_name)),] %>% 
  filter(gene_name != "") %>% 
  arrange(id)
eset <- eset[probe_MAD_dup$id, ]

# Limma
grp %>% unique()
design <- model.matrix(~0 + grp)
head(grp)
colnames(design) <- c("NT","TP") # 데이터에 맞추어 manual로 설정해야 됨

fit <- lmFit(eset,design)
cont <- makeContrasts(TP-NT,levels=design) # 데이터에 맞추어 manual로 설정해야 됨
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont)
res <- topTable(fit.cont,number=Inf) 
target <- res[res$P.Value <0.05,] %>% 
  rownames_to_column() %>% as_tibble()

# save
multiple_limma[[gse_name]] <- target
save(multiple_limma, file = "GEO_integrated.RData")
