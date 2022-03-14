# efficient motifs
1. emotif_ILP.snakefile
	split all vntr loci into `numOfSplits` files and run ILP solver on each file with 5 threads

# assembly study
1. asm_Anno.snakefile
	annotate asm with original & efficient motifs
2. asm_plot.snakefile & Assembly_alleles.ipynb
	plot histogram of # of alleles
3. visualAllele.py
	allele heatmap for vntr loci using 32 threads
	`python3 visualAllele.py lra anno-delta-0.5`
	`python3 visualAllele.py lra anno-original-motifs`

# simulation study
1. ReadAnno_ConsensusInside.snakefile 
	align reads + use `vamos --liftover` to take out the subseq of each vntr + use `vamos --conseq_anno` to annotate 

2. ReadAnno_nonConsensus.snakefile
	align reads + use `vamos --raw_anno` to annotate 

3. CompAnno.snakefile
	compare the vamos annotation with swap annotation
	NOTE: This needs to adjust manually for `nonconsensus` and `consensus-inside`, due to different read types setting.

4. Comp_plot.snakefile
	generate violin plot