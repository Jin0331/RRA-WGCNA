# function
source("src/r-function.R")

# LOCAL VARIABLE
# rdata_path <- "RData/"
# load(file = "RData/HCC_GEO_RobustDEGs_norm.RData")
pr_name <- "LIHC"
time_stamp <- Sys.time()

# weighted gene co-expression network costruction
mch_test <- purrr::quietly(mergecutheight_test)(pr_name = "LIHC", robustdegs = robustdegs)
network <- network_preprocessing(pr_name = "LIHC", robustdegs = robustdegs, mch = mch_test$result$max_mch,time_stamp = time_stamp)

# network variable
moduleLabels <- network[[2]]$colors
moduleColors <-  labels2colors(network[[2]]$colors)
MEs <- network[[2]]$MEs;
geneTree <- network[[2]]$dendrograms[[1]]
MEs0 <- moduleEigengenes(network[[1]], moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
    
# find key hub gene
key_hub_gene <- find_key_modulegene(pr_name = "LIHC", network = network, MEs = MEs, time_stamp = time_stamp)
