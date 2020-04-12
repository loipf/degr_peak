#!/usr/bin/env Rscript
library(data.table)
library(magrittr)
library(glue)
###########################
#Permutation Test for obtaining Amplified or Deleted genomic positions
###########################

#Fread mappings with source and densty amtrix
n_permutations = $params.permutations
FDR = $params.FDR

density_matrix <- readRDS("density_matrix.rds")
#read source and select genes int he order found on the density matrix
mappings_with_source_unordered = fread("peak_data.txt")
setkey(mappings_with_source_unordered, ENTREZID)
all_mappings_with_source = mappings_with_source_unordered[J(rownames(density_matrix) %>% as.numeric), nomatch = 0]

regions_list = list()
#### MASTER FOR LOOP to loop through all sources
for (source in 5:ncol(all_mappings_with_source))
{
  mappings_with_source = all_mappings_with_source[, c(1:4, source), with = FALSE]
  
  source_to_evaluate <- mappings_with_source[[5]] # last column is the source
  
  source_name = colnames(mappings_with_source)[5] 
  source_number = substr(source_name, 2, length(source_name)) %>% as.numeric()
  set_seed = 123456 + source_number
  
  permutations <- cbind(source_to_evaluate, replicate(n_permutations, sample(source_to_evaluate)))
  
  smoothened_permutations = density_matrix %*% permutations #smothening the permuted values
  
  sorted_permutations = apply(abs(smoothened_permutations), 2, function(x) sort(x, decreasing = TRUE))
  
  expected_permutation_index = 1 #initialized to minimum location
  candidate_permutation_index = 0 #initialized to 0 just to fail the while condition at start
  
  while(expected_permutation_index != candidate_permutation_index) # While FDR expected index is not equal to candidate index
  {
    #Update expectations
    candidate_permutation_index = expected_permutation_index
    
    #What is the biggest value at this candidate row (candidate threshold)
    candidate_threshold = sorted_permutations[candidate_permutation_index, -1] %>% max
    #Given this value where is the first index of source where candidate threshold is bigger
    source_index = which.max(sorted_permutations[, 1] < candidate_threshold)
    #Given this index where is the minimum row of permutations to surpass FDR
    expected_permutation_index =  ceiling(source_index * FDR)
  }
  #Source index is now the minimun index where at least one permutation surpasses FDR
  
  
  cutoffs = array(0, n_permutations)
  for(j in 2:(n_permutations + 1))
  {
    for(i in source_index:dim(sorted_permutations)[1])
    {
      current_value = sorted_permutations[i, 1]
      if(length(which(sorted_permutations[,j] > current_value))/i >= FDR)
      {
        cutoffs[j-1] = current_value
        break
      }
    }
  }
  
  #### Calling extreme valued regions
  CL = $params.FDR_CI
  state_deciding_cutoff = $params.state_deciding_cutoff
  first_cutoff = quantile(cutoffs, CL)
  
  call_region <- function(values, cutoff){
    ifelse(values > cutoff, 1, ifelse(values < -cutoff, -1, 0))
  }
  
  candidate_regions = call_region(smoothened_permutations[,1], first_cutoff)
  re_smoothened = density_matrix %*% candidate_regions
  final_regions = call_region(re_smoothened, state_deciding_cutoff)
  
  mappings_with_source[, extreme_valued_region := final_regions]
  #mappings_with_source[2:4. extreme_valued_region := 1]
  
  #print mappings with source
  
  ### Transforming extreme valued regions to bed file
  mappings_with_source[, contiguous_region_id := rleid(extreme_valued_region)]
  
  region_lenghts = mappings_with_source[, .N, by = contiguous_region_id][, N]
  region_start = mappings_with_source[, .SD[1, BP_Mapping], by = contiguous_region_id][, V1]
  region_end = mappings_with_source[, .SD[.N, BP_Mapping], by = contiguous_region_id][, V1]
  region_status = mappings_with_source[, .SD[1, extreme_valued_region], by = contiguous_region_id][, V1] # taking first status just for convenience
  region_chromosome = mappings_with_source[, .SD[1, as.character(CHR_Mapping)], by = contiguous_region_id][, V1] # taking first mapping just for convenience
  
  regions_list[[source - 4]] = data.frame(chrom = region_chromosome,
                                              chromStart = region_start,
                                              chromEnd = region_end,
                                              name = source_name,
                                              extreme_value_region_status = region_status,
                                              hyperparameter = glue("interval_coverage:{$params.interval_coverage} interval_coverage_CI:{$params.interval_coverage_CI} interval_coverage_sd_ratio:{$params.interval_coverage_sd_ratio} permutations:{$params.permutations} FDR:{$params.FDR} FDR_CI:{$params.FDR_CI} state_deciding_cutoff:{$params.state_deciding_cutoff}")
  )
}

extreme_valued_region_segments <- do.call("rbind", regions_list)

fwrite(extreme_valued_region_segments, "extreme_valued_region_segments.txt", sep = "\t")