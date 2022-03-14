# efficient motifs
1. emotif_ILP.snakefile
	This script splits all vntr loci into `numOfSplits` files and runs ILP solver on each file with 5 threads

# assembly study
1. asm_Anno.snakefile
	This script annotates asm with original & efficient motifs
2. asm_plot.snakefile & Assembly_alleles.ipynb
	This script plots histogram of # of alleles for asm vcf
3. visualAllele.py
	This script generates allele heatmaps for vntr loci based on annotation vcf using 32 threads
	`python3 visualAllele.py lra anno-delta-0.5`
	`python3 visualAllele.py lra anno-original-motifs`

# simulation study
1. ReadAnno_ConsensusInside.snakefile 
	This script aligns reads, uses `vamos --liftover` to take out the subseq of each vntr, uses `vamos --conseq_anno` to annotate.

2. ReadAnno_nonConsensus.snakefile
	This script aligns simulated reads, use `vamos --raw_anno` to annotate. 

3. CompAnno.snakefile
	This script compares the `vamos annotation` with `replacement annotation`.
	NOTE: This needs to adjust manually for `nonconsensus` and `consensus-inside`, due to different read types setting.

4. Comp_plot.snakefile
	This script generates violin plot to evaluate the `vamos annotation` VS `replacement annotation`

