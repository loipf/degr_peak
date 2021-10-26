# degr_peak

General parameters: 
```
--mapping_data	A separated text file with genomic mapping columns (gene_id, CHR_Mapping, BP_Mapping) 
--peak_data 	A separated text file with gene_ids on the first column (does not need to be named that way) 
--gene_id	Name of matching gene_id column in input files (Default gene_id) 
--id 		Prefix to attach to output file (Default No prefix) 
```

OPTIONAL:

General parameters:
```
--outdir		Directory path to print results (Default degr_results) 
--cores			Number of cores used for degr (Default 4) 
--chromosome_level	To output whole chromosome outlier calling
```

Window Parameters:
```
--interval_coverage 		Default 10 
--interval_coverage_CI		Default 0.95 
--interval_coverage_sd_ratio	Default 3 
```

Permutation test parameters:
```
--permutations              Default 1000
--permutation_strategy      Use *same* or *other* chromosomes for permutation reference (Default same)
--FDR 			    Default 0.05
--FDR_CI		    Default 0.5
--state_deciding_cutoff	    Default 0.85 
```



---
### required R-packages:
```R
install.packages("truncnorm")
install.packages("data.table")
```
other packages used but should come with R-installation:
```R
install.packages("magrittr")
install.packages("parallel")
install.packages("glue")
install.packages("assertthat")
```



---
### TACNA pipeline
```
wget https://raw.githubusercontent.com/loipf/degr_peak/master/big_machine.config
nextflow loipf/degr_peak --mapping_data MAPPING_FILE.tsv --peak_data /PATH/TO/ica_independent_components_consensus.tsv --outdir OUTDIR -c big_machine.config --cores 10 --permutations 1000 --permutation_strategy other
```



