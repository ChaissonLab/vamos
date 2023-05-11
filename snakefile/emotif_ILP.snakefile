import os
import tempfile
import subprocess
import os.path
import re

"""
This snakefile generate emotifs for vntr loci

Needs update: 

"numOfSplits" : "10000",
"delta_threshold" : ["0.1", "0.2", "0.3"],
"genomes_index" : ["64"],
"mode" : ["summary_2023-03-31"],
"""

SD = os.path.dirname(workflow.snakefile)

# Config
configfile: "/project/mchaisso_100/cmb-16/jingwenr/trfCall/vamos/snakefile/configs/snakefile.json"

input_path = config["input_path"]
numOfSplits = int(config["numOfSplits"])
splits_index = range(numOfSplits)
genomes_index = config["genomes_index"]
delta_threshold = config["delta_threshold"]
mode = config["mode"]

rule all:
	input:
		splitfa = expand("{input_path}/emotifs/{mode}/split_{genomes_index}/seq/vntr_motif.{splits_index}.tsv", \
		         mode = mode, input_path = input_path, genomes_index = genomes_index, splits_index = splits_index),
		splitILP = expand("{input_path}/emotifs/{mode}/initial-{genomes_index}-delta-{delta_threshold}/emotifs.{splits_index}.tsv", mode = mode, input_path = input_path, genomes_index = genomes_index, delta_threshold = delta_threshold, splits_index = splits_index),
		emo=expand("{input_path}/emotifs/{mode}/out-{genomes_index}-delta-{delta_threshold}/emotifs.tsv", mode = mode, input_path = input_path, genomes_index = genomes_index, delta_threshold = delta_threshold),
		vm=expand("{input_path}/emotifs/{mode}/out-{genomes_index}-delta-{delta_threshold}/vntrs_motifs.bed", mode = mode, input_path = input_path, genomes_index = genomes_index, delta_threshold = delta_threshold)

rule splitSeq:
	input:
		"{input_path}/emotifs/{mode}/aggregate.{genomes_index}.tsv"
	output:
		"{input_path}/emotifs/{mode}/split_{genomes_index}/seq/vntr_motif.{splits_index}.tsv"
	params:
		numOfSplits = config["numOfSplits"],
		grid_opts = config["grid_split"]
	shell:"""
		mkdir -p "{wildcards.input_path}/emotifs/{wildcards.mode}/split_{wildcards.genomes_index}/seq"
		awk '{{split($2, a, ","); if ((NR >= 2) && ((NR - 1) % {params.numOfSplits} == {wildcards.splits_index}) && (length(a) <= 1000)) {{print $0;}} }}' {input} > {output}
	"""

rule runILP_solver:
	input:
		"{input_path}/emotifs/{mode}/split_{genomes_index}/seq/vntr_motif.{splits_index}.tsv"
	output:
		"{input_path}/emotifs/{mode}/initial-{genomes_index}-delta-{delta_threshold}/emotifs.{splits_index}.tsv"
	params:
		ILP_Solver_py = config["ILP_Solver_py"],
		grid_opts = config["grid_ILP"]
	shell:"""
		mkdir -p "{wildcards.input_path}/emotifs/{wildcards.mode}/initial-{wildcards.genomes_index}-delta-{wildcards.delta_threshold}"
		time python3 {params.ILP_Solver_py} {input} {output} 9 {wildcards.delta_threshold}
	"""

rule combineResult:
	input:
		split_file = expand("{input_path}/emotifs/{mode}/initial-{genomes_index}-delta-{delta_threshold}/emotifs.{splits_index}.tsv", splits_index=splits_index, allow_missing=True),
		header = config["header"]
	output:
		"{input_path}/emotifs/{mode}/out-{genomes_index}-delta-{delta_threshold}/emotifs.tsv"
	shell:"""
		mkdir -p "{wildcards.input_path}/emotifs/{wildcards.mode}/out-{wildcards.genomes_index}-delta-{wildcards.delta_threshold}"
		cat {input.header} {input.split_file} > {output}
	"""

rule getVNTRBed:
	input:
		"{input_path}/emotifs/{mode}/out-{genomes_index}-delta-{delta_threshold}/emotifs.tsv"
	output:
		vm = "{input_path}/emotifs/{mode}/out-{genomes_index}-delta-{delta_threshold}/vntrs_motifs.bed",
		v = "{input_path}/emotifs/{mode}/out-{genomes_index}-delta-{delta_threshold}/vntrs.bed",
		m = "{input_path}/emotifs/{mode}/out-{genomes_index}-delta-{delta_threshold}/motifs.csv"
	shell:"""
		awk '{{split($1, a, ":|-"); if (NR > 1) {{print a[1]"\\t"a[2]"\\t"a[3]"\\t"$7;}} }}' {input} > {output.vm}
		awk '{{split($1, a, ":|-"); if (NR > 1) {{print a[1]"\\t"a[2]"\\t"a[3];}} }}' {input} > {output.v}
		awk '{{split($1, a, ":|-"); if (NR > 1) {{print$7;}} }}' {input} > {output.m}
	""" 
