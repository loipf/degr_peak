# degr_peak

General paramenters: 

--mapping_data	A separated text file with genomic mapping columns (gene_id, CHR_Mapping, BP_Mapping) \
--peak_data 	A separated text file with gene_ids on the first column. \
--gene_id		Name of matching gene_id column in input files (Default gene_id) \
--id 			Prefix to attach to output file (Default No prefix) 

OPTIONAL:

General paramenters:

--outdir		Directory path to print results (Default degr_results) 
--cores			Number of cores used for degr (Default 4) 

Window Parameters:

--interval_coverage 			Default 10 \
--interval_coverage_CI			Default 0.95 \
--interval_coverage_sd_ratio	Default 3 

Permutation test parameters:

--permutations				Default 1000 \
--FDR 						Default 0.05 \
--FDR_CI = 0.5 \
--state_deciding_cutoff		Default 0.85 

--chromosome_level		To output whole chromosome outlier calling
