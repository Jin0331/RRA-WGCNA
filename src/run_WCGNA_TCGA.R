# function
source("src/r-function.R")

# load RData
if(!exists("robustdegs")){
  load(paste0(base_dir, "/Step2_RRA_DEG.RData"))
}

# weighted gene co-expression network costruction
mch_test <- purrr::quietly(mergecutheight_test)(pr_name = pr_name, robustdegs = robustdegs) 
mch_test_value <- mch_test$result$max_mch

network <- network_preprocessing(pr_name = pr_name, robustdegs = robustdegs, 
                                 mch = mch_test_value,
                                 time_stamp = time_stamp)

# network variable
moduleLabels <- network[[2]]$colors
moduleColors <-  labels2colors(network[[2]]$colors)
MEs <- network[[2]]$MEs;
geneTree <- network[[2]]$dendrograms[[1]]
MEs0 <- moduleEigengenes(network[[1]], moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
    
# find key hub gene
key_hub_gene <- find_key_modulegene(pr_name = pr_name, base_dir = base_dir,
                                    network = network, MEs = MEs, 
                                    select_clinical = c("hist_hepato_carc_fact"),
                                    mm = 0.85, gs = 0.2)

save(key_hub_gene, file = paste0(base_dir, "/Step3_KEY_HUB_GENE.RData"))

# Gene Onthlogy

library(clusterProfiler)
library(org.Hs.eg.db)

geneList <- key_hub_gene$intra_module_gene$turquoise$GS
names(geneList) <- bitr(key_hub_gene$intra_module_gene$turquoise$gene, 
                        fromType="ALIAS", 
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db, 
                        drop = FALSE) %>% 
  distinct(ALIAS, .keep_all = TRUE) %>% 
  pull(2)

geneList <- geneList %>% sort(decreasing = T)

# ORA
## BP
ego_ora_bp <- enrichGO(names(geneList), 
                       OrgDb=org.Hs.eg.db, 
                       ont="BP", 
                       pAdjustMethod = "bonferroni",
                       pvalueCutoff=0.01,
                       qvalueCutoff=0.05,
                       readable=T)

## CC
ego_ora_cc <- enrichGO(names(geneList), 
                        OrgDb=org.Hs.eg.db, 
                        ont="CC", 
                        pAdjustMethod = "bonferroni",
                        pvalueCutoff=0.01,
                        qvalueCutoff=0.05,
                        readable=T)

# MF
ego_ora_mf <- enrichGO(names(geneList), 
                       OrgDb=org.Hs.eg.db, 
                       ont="MF", 
                       pAdjustMethod = "bonferroni",
                       pvalueCutoff=0.01,
                       qvalueCutoff=0.05,
                       readable=T)

kk <- enrichKEGG(gene = names(geneList),
                 organism = 'hsa',
                 keyType = "SYMBOL",
                 pvalueCutoff = 0.05)

library(enrichplot)
barplot(ego_ora_mf, showCategory=30)
