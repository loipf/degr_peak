#!/usr/bin/env Rscript
library(data.table)
library(magrittr)
library(glue)

#Load mapping and peaks
CES_data = fread("peak_data.txt", header = TRUE)
gene_info = fread("mapping_data.txt") %>%
  .[top_probe_indicator == 1,] %>% #only jetset top probes
  .[order(BP_Mapping), .(ENTREZID, PROBESET, CHR_Mapping, BP_Mapping)] #required columns

#generate peak file with mapping
CES_data_with_mapping = gene_info[CES_data, on = "ENTREZID==rn", nomatch = 0]

for (chromosome in CES_data_with_mapping[, unique(CHR_Mapping)]){
  fwrite(CES_data_with_mapping[CHR_Mapping == chromosome, 1:4], glue("mapping_{chromosome}.txt"), sep = "\t")
}

fwrite(CES_data_with_mapping, glue("peaks_data.txt"), sep = "\t")
