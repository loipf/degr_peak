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
mappings_with_source <- fread("source.txt")
density_matrix <- readRDS("density_matrix.rds")


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

region_lenghts = mappings_with_source[, .N, by = contiguous_region_id]
region_start = mappings_with_source[, .SD[1, BP_Mapping], by = contiguous_region_id]
region_end = mappings_with_source[, .SD[.N, BP_Mapping], by = contiguous_region_id]
region_status = mappings_with_source[, .SD[1, extreme_valued_region], by = contiguous_region_id] # taking first status just for convenience
region_chromosome = mappings_with_source[, .SD[1, CHR_Mapping], by = contiguous_region_id] # taking first mapping just for convenience

extreme_valued_region_segments = data.table(chrom = region_chromosome,
                                            chromStart = region_start,
                                            chromEnd = region_end,
                                            name = source_name,
                                            extreme_value_region_status = region_status,
                                            hyperparameter = glue("interval_coverage:{$params.interval_coverage} interval_coverage_CI:{$params.interval_coverage_CI} $params.interval_coverage_sd_ratio:{$params.interval_coverage_sd_ratio} FDR:{$params.FDR} FDR_CI:{$params.FDR_CI} state_deciding_cutoff:{$params.state_deciding_cutoff}")
)

fwrite(extreme_valued_region_segments, "extreme_valued_region_segments.txt", sep = "\t")