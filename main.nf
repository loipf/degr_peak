params.outdir = "results_degr/"
params.verbose = null
params.cores = 4
params.id = ""


process attach_mapping_data {
  label 'rscript'
  label 'regular'
  cpus = params.cores
  echo {params.verbose != null ? true : false}
  
  input:
    file "peak_data.txt" from Channel.fromPath("$params.peak_data")
    file "mapping_data.txt" from Channel.fromPath("$params.mapping_data")
  
  output:
    file "mapping_*.txt" into mappings
    file "peaks_data.txt" into peaks_data

  
  script:
    template 'attach_mapping_data.R'
  
}

peaks_data.
into {peak_data_evc; peak_data_evr}

process smoothing_matrix {
  label 'rscript'
  label 'regular'
  cpus = params.cores
  echo {params.verbose != null ? true : false}
  
  input:
    file "peaks_with_mapping.txt" from mappings.flatten()
  
  output:
    file "density_matrix.rds" into density_matrix
  
  script:
    template 'smoothing_matrix.R'
  
}

process degr_peak_calling {
  label 'rscript'
  label 'regular'
  cpus = params.cores
  echo {params.verbose != null ? true : false}
  
  input:
    set "density_matrix.rds", "peak_data.txt" from density_matrix.combine(peak_data_evr)
  
  output:
    file "extreme_valued_region_segments.txt" into extreme_valued_regions
  
  script:
    template 'degr_peak_calling.R'
    
}

process extreme_valued_chromosomes {
  label 'rscript'
  label 'regular'
  cpus = params.cores
  echo {params.verbose != null ? true : false}
  
  input:
    file "peak_data.txt" from peak_data_evc
  
  output:
    file "extreme_valued_chromosomes.txt" into extreme_valued_chromosomes

  when:
  params.chromosome_level != "NOT_PROVIDED"
  
  script:
    template 'extreme_chromosomes.R'
  
}


extreme_valued_regions
.collectFile(name: 'extreme_valued_regions_all_chromosomes.txt', skip: 1, keepHeader: true)
.set{ extreme_valued_regions_all_chromosomes }

extreme_valued_regions_all_chromosomes
.mix(extreme_valued_chromosomes)
.set{ all_results }


process print_results {
  label 'rscript'
  label 'regular'
  echo {params.verbose != null ? true : false}
  
  publishDir "${params.outdir}/", mode: 'copy', saveAs: { filename -> "${params.id}_$filename" }
      
  input:
    file result from all_results
  
  output:
    file result
  
  """
  echo Printing
  """
  }
