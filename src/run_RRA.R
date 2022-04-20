library(RobustRankAggreg)
library(pheatmap)
library(tidyverse)

# function
source("src/function.R")

# load
load("~/WORK/GEO/RData/GEO_integrated.RData")
rra_result <- rra_extract(m_list = multiple_limma, 
                          fdr = 0.0
                          )

combine_degs_rra <- rra_result[[1]] %>% 
  filter(GENE %in% c(rra_result[[2]] %>% pull(1) %>% .[1:20],
                     rra_result[[3]] %>% pull(1) %>% .[1:20]))
up_down_rra_gene <- bind_rows(rra_result[[2]], rra_result[[3]])

# heatmap
combine_degs_m <- combine_degs_rra[,-1] %>% as.matrix()
rownames(combine_degs_m) <- combine_degs_rra$GENE
colnames(combine_degs_m) <- colnames(combine_degs_rra)[2:length(colnames(combine_degs_rra))]

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
save(up_down_rra_gene, file = "RobustDEGs.RData")
