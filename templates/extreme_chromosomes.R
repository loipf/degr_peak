#!/usr/bin/env Rscript
library(data.table)
library(magrittr)
library(glue)
library(parallel)
library(assertthat)
###########################
#Permutation Test for obtaining Amplified or Deleted chromosomes
#stategy is to generate permutations and generate a p value for perms that got higher than real data
###########################

#Fread mappings with source and densty amtrix
n_permutations = $params.permutations

all_mappings_with_source = fread("peak_data.txt")
setkey(all_mappings_with_source, CHR_Mapping)

clust <- makeCluster($params.cores, type="FORK")
#### MASTER FOR LOOP to loop through all sources
regions_list <- parLapply(clust, 4:ncol(all_mappings_with_source), function(source_colnumber)
{
  source_name = colnames(all_mappings_with_source)[source_colnumber]
  source_values = all_mappings_with_source[,source_colnumber, with=F][[1]]
  
  #generate matrix with original data and permutations
  set_seed = 123456 + source_colnumber
  set.seed(set_seed)

  original_and_permutations = cbind(source_values, replicate(n_permutations, sample(source_values))) %>% as.data.table
  original_and_permutations[, CHR_Mapping := all_mappings_with_source[,CHR_Mapping]]
  
  #calculate medians per chromosome for real and permutations
  chromosome_level_medians <- original_and_permutations[, lapply(.SD, median), by = CHR_Mapping]

  source_chrom_values <- chromosome_level_medians[, 1:2]
  setnames(source_chrom_values, c("CHR_Mapping", "original_chr_values"))

  #Sorting both separately with dummy to be explicit that same method must be applied to not break code
  chromosome_level_medians_mtx <- chromosome_level_medians %>% .[,-1] %>% as.matrix
  dummy_abs <- function(x) abs(x)
  #sorting separately real data and permutations to call p values
  ordering_of_source = order(dummy_abs(chromosome_level_medians_mtx[,1]), decreasing = TRUE)
  sorted_permutations = apply(dummy_abs(chromosome_level_medians_mtx), 2, function(x) sort(x, decreasing = TRUE))

  # call p value when with the num of permutation higher than real data
  p_values = array(0, dim(sorted_permutations)[1])
  for(i in 1:length(p_values))
  {
      current_value = sorted_permutations[i, 1]
      p_values[i] = sum(sorted_permutations[i,-1] > current_value)/n_permutations
  }
  

  chromosome_status = source_chrom_values[ordering_of_source, original_chr_values/abs(original_chr_values)]
  chromosome_name = source_chrom_values[ordering_of_source, CHR_Mapping]
  chromosome_pvalues = p_values
  chromosome_magnitude = sorted_permutations[, 1]
  
  return(data.frame(chromosome = chromosome_name,
                    chrom_pvalue = chromosome_pvalues,
                    extreme_value_magnitude = chromosome_magnitude,
                    extreme_value_status = chromosome_status,
                    name = source_name
        ))
})
stopCluster(clust)

extreme_valued_region_chromosomes <- do.call("rbind", regions_list) %>% as.data.table

fwrite(extreme_valued_region_chromosomes, "extreme_valued_chromosomes.txt", sep = "\t")
