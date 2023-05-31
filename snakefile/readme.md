# Efficient motifs
1. emotif_ILP.snakefile
	This script splits all vntr loci into `numOfSplits` files and runs ILP solver on each file with 9 threads. 

# Comparison with SV-based method (Jasmine)
1. process_before_sv-vntr.snakefile
	This script filters SNPs and indels under specified length.

2. separate_Jasmine_merge_by_chrom_indel20.snakefile
	This script separates SVs by chromosomes and Jasmine to merge SVs across samples. 

3. sv-vntr_combine.snakefile
	This script compares SV alleles overlapping with VNTRs with vamos alleles. 

# Plot the heatmap-plot for special loci 
1. plot_special_loci.snakefile
	This script extracts VCF for specific VNTR loci. 

# Comparison with Greedy method
1. vamos_greedy_compare.snakefile
	This script compares the results from greedy and vamos methods. 
