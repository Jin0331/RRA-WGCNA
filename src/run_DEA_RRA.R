# function
source("src/r-function.R")

# ======= RUN DEA & RRA ======= 
# save
multiple_limma <- GSE_manual()
save(multiple_limma, file = paste0(base_dir, "/Step1_RAW.RData"))

# rra & robust DEGs
if(!exists("multiple_limma"))
  load(paste0(base_dir, "/Step2_RAW.RData"))

robustdegs <- rra_analysis(m_list = multiple_limma, logfc = 0, fdr = 0.05, save_path = base_dir)
save(robustdegs, file = paste0(base_dir, "/Step2_RRA_DEG.RData"))