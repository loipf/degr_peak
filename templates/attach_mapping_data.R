#!/usr/bin/env Rscript
library(data.table)
library(magrittr)
library(glue)
library(assertthat)

#Load mapping and peaks
CES_data = fread("peak_data.txt", header = TRUE)
gene_info = fread("mapping_data.txt") %>% complete.cases %>%
  .[order(BP_Mapping), .($params.gene_id, CHR_Mapping, BP_Mapping)] #required columns

#chack that each gene has only one mapping
assert_that(gene_info[, anyDuplicated($params.gene_id)] == 0)

CES_data[, $params.gene_id := as.character($params.gene_id)]
gene_info[, $params.gene_id := as.character($params.gene_id)]

#generate peak file with mapping
CES_data_with_mapping = gene_info[CES_data, on = "$params.gene_id==$params.gene_id", nomatch = 0]
 #check that each there is some matches between peak data and mapping data
  assert_that(nrow(CES_data_with_mapping) > 0, msg="Zero matches between mapping data and peak data... using column: $params.gene_id")

for (chromosome in CES_data_with_mapping[, unique(CHR_Mapping)]){
  fwrite(CES_data_with_mapping[CHR_Mapping == chromosome, 1:3], glue("mapping_{chromosome}.txt"), sep = "\t")
}

fwrite(CES_data_with_mapping, glue("peaks_data.txt"), sep = "\t")
