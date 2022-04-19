library(RobustRankAggreg)
library(pheatmap)
library(tidyverse)

# function
rra_extract <- function(m_list, logfc = 0.0, fdr = 0.05){
  # combine deg
  combine_degs <- names(m_list) %>% 
    lapply(X = ., FUN = function(list_name){
      tmp <- multiple_limma[[list_name]] %>% 
        filter(adj.P.Val < fdr & (logFC > logfc | logFC < -(logfc))) %>% 
                 arrange(desc(logFC)) %>% 
                 select(rowname, logFC)
               colnames(tmp) <- c("GENE", list_name)
               return(tmp)
    }) %>% reduce(., left_join, by = "GENE")
  
  # up-regulated
  up_degs <- names(m_list) %>% 
    lapply(X = ., FUN = function(list_name){
      m_list[[list_name]] %>% 
        filter(adj.P.Val < fdr & logFC > logfc) %>%
        pull(rowname) %>%
        return()
    }) 
      
  # down-regulated
  down_degs <- names(m_list) %>% 
    lapply(X = ., FUN = function(list_name){
      m_list[[list_name]] %>% 
        filter(adj.P.Val < fdr & logFC < -(logfc)) %>% 
        pull(rowname) %>% 
        return()
    }) 
      
  # Aggregate the inputs
  # run RRA
  up_deg_rra <- aggregateRanks(glist = up_degs, method = "RRA") %>%
    as_tibble() %>% 
    filter(Score < 0.05) %>% 
    arrange(Score)
  
  down_deg_rra <- aggregateRanks(glist = down_degs, method = "RRA") %>%
    as_tibble() %>% 
    filter(Score < 0.05) %>% 
    arrange(Score)
  
  # 1 - combine deg, 2 - up-regulated RRA 3 - down-regulated RRA
  list(combine_deg = combine_degs, up_deg_rra = up_deg_rra, down_deg_rra = down_deg_rra) %>% 
    return()
}


# load
load("~/WORK/GEO/GEO_integrated.RData")
rra_result <- rra_extract(m_list = multiple_limma)

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
