# BiocManager::install("geometadbfile <- getSQLiteFile()")
# BiocManager::install("GEOmetadb")

library(GEOmetadb)
library(tidyverse)

# geometadb connect
# geometadbfile <- getSQLiteFile() # geometadb download
geometabdbfile <- "GEOmetadb.sqlite"

con <- dbConnect(SQLite(), geometabdbfile)
geo_tables <- dbListTables(con)

gse <- tbl(con, "gse") %>% collect()
gsm <- tbl(con, "gsm") %>% collect()



#
gse_list <- c('GSE3500','GSE4108','GSE6222','GSE7474','GSE27150','GSE40873','GSE45114','GSE117361','GSE166163')
gsm %>% filter(series_id == gse_list[9]) %>% View()
# gse %>% filter(gse == gse_list[4]) %>% View()

gsm %>% filter(series_id == gse_list[9]) %>% 
  write_delim(file = paste0("gsm/", gse_list[9], "_sample.tsv"), delim = "\t")
