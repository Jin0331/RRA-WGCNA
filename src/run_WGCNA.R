library(WGCNA)
library(GEOquery)
library(tidyverse)

# function
source("src/function.R")

# base data
gse_name <- "GSE14520"
gse <- getGEO(gse_name,GSEMatrix=TRUE)
gse <- gse[[1]]