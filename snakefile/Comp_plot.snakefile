import os
import tempfile
import subprocess
import os.path
import re


# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
# Config
configfile: "/project/mchaisso_100/cmb-16/jingwenr/trfCall/script/CompAnno.json"

input_path = config["input_path"] 
sim_path = config["sim_path"] 
emotif_path = config["emotif_path"]

# vntr_idx = ["{0:03}".format(i) for i in range(1, 88)]
genomes_index = config["genomes_index"]
delta_threshold = config["delta_threshold"]
gt_sim_mode = config["gt_sim_mode"]
read_sim_mode = config["read_sim_mode"]
mode = config["mode"]

numOfVntrs = int(config["numOfVntrs"])
vntr_idx = [str(i) for i in range(1, numOfVntrs + 1)]

# read_type_dict = {"CCS" : "ccs_0001.fastq", "CLR" : "clr_0001.fastq", "ONT" : "fa", "PerfectRead" : "combined.fasta", }
read_type_dict = {"CCS" : "ccs_0001.fastq", "ONT" : "fa"}
read_type = ["CCS", "ONT", "GroundTruth"]

consensus = [ "consensus-inside", "nonconsensus"]
software = ["comp"]


rule all:
    input:
        # edist_plot = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/plot/editSwapVamos.delta-{delta_threshold}-{software}.{consensus}.png", \
              # input_path=input_path, mode=mode, genomes_index=genomes_index, gt_sim_mode=gt_sim_mode, delta_threshold=delta_threshold, consensus=consensus, software=software),
        p1_png = expand("{input_path}/{mode}/result-{genomes_index}/plot/editSwapVamos.delta-{delta_threshold}-{software}.consensus-inside.group.png", \
              input_path=input_path, mode=mode, genomes_index=genomes_index, delta_threshold=delta_threshold, software=software),
        p2_png = expand("{input_path}/{mode}/result-{genomes_index}/plot/editSwapVamos.delta-{delta_threshold}-{software}.nonconsensus.group.png", \
                    input_path=input_path, mode=mode, genomes_index=genomes_index, delta_threshold=delta_threshold, software=software),
        p1_eps = expand("{input_path}/{mode}/result-{genomes_index}/plot/editSwapVamos.delta-{delta_threshold}-{software}.consensus-inside.group.eps", \
              input_path=input_path, mode=mode, genomes_index=genomes_index, delta_threshold=delta_threshold, software=software),
        p2_eps = expand("{input_path}/{mode}/result-{genomes_index}/plot/editSwapVamos.delta-{delta_threshold}-{software}.nonconsensus.group.eps", \
                    input_path=input_path, mode=mode, genomes_index=genomes_index, delta_threshold=delta_threshold, software=software)

rule edist_plot:
  output:
    "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/plot/editSwapVamos.delta-{delta_threshold}-{software}.{consensus}.png"
  params:
    editSwapVamos_py = config["editSwapVamos_py"]
  shell:"""
    mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/plot/"
    python3 {params.editSwapVamos_py} {wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/comp {output} {wildcards.consensus} {wildcards.delta_threshold} {wildcards.gt_sim_mode} {wildcards.software}
  """

rule group_edist_plot:
  params:
    together_editSwapVamos_py = config["together_editSwapVamos_py"]
  output:
    p1_png = "{input_path}/{mode}/result-{genomes_index}/plot/editSwapVamos.delta-{delta_threshold}-{software}.consensus-inside.group.png",
    p2_png = "{input_path}/{mode}/result-{genomes_index}/plot/editSwapVamos.delta-{delta_threshold}-{software}.nonconsensus.group.png",
    p1_eps = "{input_path}/{mode}/result-{genomes_index}/plot/editSwapVamos.delta-{delta_threshold}-{software}.consensus-inside.group.eps",
    p2_eps = "{input_path}/{mode}/result-{genomes_index}/plot/editSwapVamos.delta-{delta_threshold}-{software}.nonconsensus.group.eps"
  shell:"""
    mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/plot"
    python3 {params.together_editSwapVamos_py} {wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index} {output.p1_png} {output.p2_png} {output.p1_eps} {output.p2_eps} {wildcards.delta_threshold} {wildcards.software}
  """ 

