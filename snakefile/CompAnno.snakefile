import os
import tempfile
import subprocess
import os.path
import re

"""
This snakefile compare the swap-anno of the ground truth VS the vamos-anno of the read bam
"""

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
# Config
configfile: "/project/mchaisso_100/cmb-16/jingwenr/trfCall/script/CompAnno.json"

input_path = config["input_path"] 
sim_path = config["sim_path"] 
emotif_path = config["emotif_path"]

genomes_index = config["genomes_index"]
delta_threshold = config["delta_threshold"]
gt_sim_mode = config["gt_sim_mode"]
read_sim_mode = config["read_sim_mode"]
mode = config["mode"]
numOfVntrs = int(config["numOfVntrs"])
vntr_idx = [str(i) for i in range(1, numOfVntrs + 1)]

# read_type_dict = {"CCS" : "ccs_0001.fastq", "CLR" : "clr_0001.fastq", "ONT" : "fa", "PerfectRead" : "combined.fasta", }

"""
NOTE: (manual) for consensus-inside
"""
read_type_dict = {"CCS" : "ccs_0001.fastq", "ONT" : "fa"}
read_type = ["CCS", "ONT"]
consensus = ["consensus-inside"]

"""
NOTE: (manual) for nonconsensus
"""
read_type_dict = {"CCS" : "ccs_0001.fastq", "ONT" : "fa", "PerfectRead" : "combined.fasta"}
read_type = ["CCS", "ONT", "GroundTruth"]
consensus = ["nonconsensus"]


rule all:
    input:
        read_combine_anno_vcf = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/comp/read-{read_type}/comp-delta-{delta_threshold}.{consensus}/vamos.anno.{read_sim_mode}.delta-{delta_threshold}.{consensus}.vcf", \
        			consensus = consensus, mode = mode, read_type = read_type, read_sim_mode = read_sim_mode, genomes_index = genomes_index, delta_threshold = delta_threshold, input_path = input_path, gt_sim_mode = gt_sim_mode),
        attach = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/gt/vntr.gt.{vntr_idx}.{read_sim_mode}.bed", \
        			mode = mode, read_sim_mode = read_sim_mode, genomes_index = genomes_index, gt_sim_mode = gt_sim_mode, input_path = input_path, vntr_idx = vntr_idx),
        gt_combined_anno = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/gt/vntr.gt.anno.{read_sim_mode}.vcf", \
        			mode = mode, genomes_index = genomes_index, gt_sim_mode = gt_sim_mode, input_path = input_path, read_sim_mode = read_sim_mode),
        swap_gt = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/gt/vntr.gt.swap.anno.{read_sim_mode}.delta-{delta_threshold}.bed", \
        			mode = mode, read_type = read_type, genomes_index = genomes_index, gt_sim_mode = gt_sim_mode, delta_threshold = delta_threshold, input_path = input_path, read_sim_mode = read_sim_mode),
        log = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/comp/read-{read_type}/comp-delta-{delta_threshold}.{consensus}/log.{read_sim_mode}.txt", \
        			consensus = consensus, mode = mode, read_type = read_type, genomes_index = genomes_index, delta_threshold = delta_threshold, input_path = input_path, gt_sim_mode = gt_sim_mode, read_sim_mode = read_sim_mode),

"""
reads vamos anno vcf
"""
rule CombineReadVCFs:
	input:
		vcfheader = config["vcfheader"],
		vcf = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/vamos/read-{read_type}/locus{vntr_idx}/tmp/vamos.anno.{read_sim_mode}.delta-{delta_threshold}.noheader.{consensus}.vcf", \
					vntr_idx = vntr_idx, allow_missing = True)
	output:
		"{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/comp/read-{read_type}/comp-delta-{delta_threshold}.{consensus}/vamos.anno.{read_sim_mode}.delta-{delta_threshold}.{consensus}.vcf"
	params:
		combineReadVCFsh = config["combineReadVCFsh"]
	shell:"""
		cat {input.vcfheader} {input.vcf} > {output}
	"""

# VNTR loc hg38, original motif anno, sim seq name
rule AttachToVNTRHg38:
	input:
		aggregate =  emotif_path + "/{mode}/aggregate.{genomes_index}.tsv",
		gt = sim_path + "/{gt_sim_mode}/set{genomes_index}/locus{vntr_idx}/vntr.locus{vntr_idx}{read_sim_mode}.truth.bed",
	output:
		"{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/gt/vntr.gt.{vntr_idx}.{read_sim_mode}.bed"
	shell:"""
		mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/gt"
		awk -v vntr={wildcards.vntr_idx} 'BEGIN{{OFS="\\t"}}{{ if (NR==FNR) {{coord[FNR]=$1;}} else {{print coord[vntr+1], $2, $1;}} }}' {input.aggregate} {input.gt}  > {output}
	"""

rule CombineGTVCFs:
	input:
		expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/gt/vntr.gt.{vntr_idx}.{read_sim_mode}.bed", vntr_idx = vntr_idx, allow_missing = True)
	output:
		"{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/gt/vntr.gt.anno.{read_sim_mode}.vcf"
	shell:"""
		cat {input} > {output}
	"""

rule SwapAnno:
	input:
		emotifs = emotif_path + "/{mode}/out-{genomes_index}-delta-{delta_threshold}/emotifs.tsv",
		gt_bed = "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/gt/vntr.gt.anno.{read_sim_mode}.vcf"
	output:
		"{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/gt/vntr.gt.swap.anno.{read_sim_mode}.delta-{delta_threshold}.bed"
	params:
		swapAnno = config["swapAnno_py"]
	shell: """
		python3 {params.swapAnno} {input.emotifs} {input.gt_bed} {output}
	"""

rule AnnoComp:
	input:
		gt_bed = "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/gt/vntr.gt.swap.anno.{read_sim_mode}.delta-{delta_threshold}.bed",
		read_vcf = "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/comp/read-{read_type}/comp-delta-{delta_threshold}.{consensus}/vamos.anno.{read_sim_mode}.delta-{delta_threshold}.{consensus}.vcf"
	params:
		anno_comp = config["AnnoComp_py"]
	log: 
		"{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/comp/read-{read_type}/comp-delta-{delta_threshold}.{consensus}/log.{read_sim_mode}.txt"
	shell:"""
		mkdir -p "{input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/comp/read-{wildcards.read_type}/comp-delta-{wildcards.delta_threshold}.{wildcards.consensus}/"
		python3 {params.anno_comp} {input.gt_bed} {input.read_vcf} {wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/comp/read-{wildcards.read_type}/comp-delta-{wildcards.delta_threshold}.{wildcards.consensus}/gt_reads.swap.comp.anno {wildcards.delta_threshold} {wildcards.read_type}
		echo "done!" > {log}
	"""




