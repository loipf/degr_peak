#!/usr/bin/env Rscript
library(data.table)
library(magrittr)
library(glue)

peaks_with_mapping = fread("peaks_with_mapping.txt") %>%
  .[order(BP_Mapping)]

########################
#Finding proper interval for sliding Gaussian Kernel
#Choosing that interval where 10 or more no. of probesets 
#corresponding to that chromosome are there in +/- 3*interval for 95% of the cases
########################

mappings <- peaks_with_mapping[, BP_Mapping]
num_of_mappings = length(mappings)

interval_coverage = $params.interval_coverage
interval_coverage_CI = $params.interval_coverage_CI

linear_nearest_n <- function(interval_coverage, mappings, order_of_mappings) {
  size_mappings <- length(mappings)
  nearest_neighbor_distances <- matrix(NA, interval_coverage, size_mappings)
  
  distance_buffer <- rep(NA, interval_coverage)
  round_counter = 0
  previous_mapping = NA
  #Prefill variable first mapping. Exceptional case without previous mapping
  previous_mapping = mappings[1]
  
  if(order_of_mappings == "forward"){
    sequence = 1:size_mappings
  } else if(order_of_mappings == "reverse"){
    sequence = rev(1:size_mappings)
  }
  
  #First pass left to right
  for (n_mapping in sequence)
  {
    current_mapping = mappings[n_mapping]
    
    distance_buffer = distance_buffer + abs(current_mapping - previous_mapping) # add new covered distance
    distance_buffer[round_counter + 1] = 0 #reset mapping that is to far away already
    
    #save to matrix
    nearest_neighbor_distances[,n_mapping] <- distance_buffer
    
    previous_mapping = current_mapping
    round_counter = (round_counter + 1) %% interval_coverage
  }
  
  return(nearest_neighbor_distances)
}

selected_quantile <- rbind(linear_nearest_n(interval_coverage, mappings, "forward"),
                           linear_nearest_n(interval_coverage, mappings, "reverse")) %>%
  apply(2, sort) %>% # order distances low to high columnwise
  lapply(function(x) x[interval_coverage + 1]) %>% unlist %>% # nth smallest distance guarantes k nearest neighbors (+1 bcs current gene is coutned twice)
  quantile(interval_coverage_CI) # coverage percentage guarantee

###########################
#Creating the Density Matrix to get sliding gaussian Kernel
###########################
library(truncnorm)

interval_coverage_sd_ratio = $params.interval_coverage_sd_ratio

density_matrix = matrix(0, num_of_mappings, num_of_mappings)

for(j in 1:num_of_mappings)
{
  #generate truncated normal distribution at each mapping
  den_tnorm = dtruncnorm(mappings, #use all ordered bp mappings as template
                         a = mappings[1], #left truncate at first mapping
                         b = mappings[num_of_mappings], #right truncate at last mapping
                         mean = mappings[j], # centered on current point 
                         sd = selected_quantile / interval_coverage_sd_ratio) # with standard deviation equal to the selected interval / 3
  
  density_matrix[j,] = den_tnorm/sum(den_tnorm)
}

for (n_source in 5:ncol(peaks_with_mapping))
{
  mappings_with_source <- peaks_with_mapping[, c(1:4, n_source), with = FALSE]
  save(density_matrix, mappings_with_source, file = glue("peak_call_input{n_source - 4}.RData"))
}