params.outdir = "results_degr/"
params.verbose = null


process attach_mapping_data {
  label 'rscript'
  label 'regular'
  cpus 8
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

process smoothing_matrix {
  label 'rscript'
  label 'regular'
  cpus 8
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
  echo {params.verbose != null ? true : false}
  
  input:
    set "density_matrix.rds", "peak_data.txt" from density_matrix.combine(peaks_data)
  
  output:
    file "extreme_valued_region_segments.txt" into extreme_valued_regions
  
  script:
    template 'degr_peak_calling.R'
  
}

extreme_valued_regions
.collectFile(name: 'extreme_valued_regions_all_chromosomes.txt', skip: 1, keepHeader: true)
.set{ extreme_valued_regions_all_chromosomes }

process print_results {
  label 'rscript'
  label 'regular'
  echo {params.verbose != null ? true : false}
  
  publishDir "${params.outdir}/", mode: 'copy', saveAs: { filename -> "${params.id}_$filename" }
      
  input:
    file result from extreme_valued_regions_all_chromosomes
  
  output:
    file result
  
  """
  echo Printing
  """
  }