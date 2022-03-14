import os
import tempfile
import subprocess
import os.path
import re

"""
This snakefile generates vamos annotations for read bam (inside-consensus, non-consensus)
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
# consensus = config["consensus"]

# read_type_dict = {"CCS" : "ccs_0001.fastq", "CLR" : "clr_0001.fastq", "ONT" : "fa", "PerfectRead" : "combined.fasta"}
read_type_dict = {"CCS" : "ccs_0001.fastq", "ONT" : "fa", "PerfectRead" : "combined.fasta"}
read_type = ["CCS", "ONT", "GroundTruth"]
consensus = ["nonconsensus"]

rule all:
  input:
    read_bam = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/align/read-{read_type}/locus{vntr_idx}/reads-{read_sim_mode}.nonconsensus.bam", \
      mode = mode, read_type = read_type, read_sim_mode = read_sim_mode, genomes_index = genomes_index, input_path = input_path, vntr_idx = vntr_idx, gt_sim_mode = gt_sim_mode),
    emotif = expand("{input_path}/{mode}/result-{genomes_index}/emotif_{delta_threshold}/motif-{vntr_idx}.csv", \
      mode = mode, genomes_index = genomes_index, delta_threshold = delta_threshold, input_path = input_path, vntr_idx = vntr_idx), 
    read_anno_vcf = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/vamos/read-{read_type}/locus{vntr_idx}/vamos.anno.{read_sim_mode}.delta-{delta_threshold}.{consensus}.vcf",
      mode = mode, read_type = read_type, read_sim_mode = read_sim_mode, genomes_index = genomes_index, delta_threshold = delta_threshold, input_path = input_path, vntr_idx = vntr_idx, gt_sim_mode = gt_sim_mode, consensus=consensus),                        
    tmp = expand("{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/vamos/read-{read_type}/locus{vntr_idx}/tmp/vamos.anno.{read_sim_mode}.delta-{delta_threshold}.noheader.{consensus}.vcf",
      mode = mode, read_type = read_type, read_sim_mode = read_sim_mode, genomes_index = genomes_index, delta_threshold = delta_threshold, input_path = input_path, vntr_idx = vntr_idx, gt_sim_mode = gt_sim_mode, consensus=consensus), 

###############################
########### ALIGN #############
###############################
rule AlignPerfectReads:
  input:
    read = lambda wc: sim_path + f"/{wc.gt_sim_mode}/set{wc.genomes_index}/locus{wc.vntr_idx}/reads.locus{wc.vntr_idx}{wc.read_sim_mode}.combined.fasta",
    ref = sim_path + "/{gt_sim_mode}/set{genomes_index}/locus{vntr_idx}/vntr.locus{vntr_idx}{read_sim_mode}.fa"
  output:
    "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/align/read-PerfectRead/locus{vntr_idx}/reads-{read_sim_mode}.nonconsensus.bam"
  params:
    minimap2 = config["minimap2"],
    samtools = config["samtools"]
  shell: """
    mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/align/read-PerfectRead/locus{wildcards.vntr_idx}/"
    {params.minimap2} -ax map-pb {input.ref} {input.read} | {params.samtools} sort > {output}
    {params.samtools} index {output} 
  """

rule AlignGroundTruth:
  input:
    read = lambda wc: sim_path + f"/{wc.gt_sim_mode}/set{wc.genomes_index}/locus{wc.vntr_idx}/vntr.locus{wc.vntr_idx}{wc.read_sim_mode}.fa",
    ref = sim_path + "/{gt_sim_mode}/set{genomes_index}/locus{vntr_idx}/vntr.locus{vntr_idx}{read_sim_mode}.fa"
  output:
    "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/align/read-GroundTruth/locus{vntr_idx}/reads-{read_sim_mode}.nonconsensus.bam"
  params:
    minimap2 = config["minimap2"],
    samtools = config["samtools"]
  shell: """
    mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/align/read-GroundTruth/locus{wildcards.vntr_idx}/"
    {params.minimap2} -ax map-pb {input.ref} {input.read} | {params.samtools} sort > {output}
    {params.samtools} index {output} 
  """


rule AlignONTReads:
  input:
    read = lambda wc: sim_path + f"/{wc.gt_sim_mode}/set{wc.genomes_index}/locus{wc.vntr_idx}/reads.locus{wc.vntr_idx}{wc.read_sim_mode}.fa",
    ref = sim_path + "/{gt_sim_mode}/set{genomes_index}/locus{vntr_idx}/vntr.locus{vntr_idx}{read_sim_mode}.fa"
  output:
    "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/align/read-ONT/locus{vntr_idx}/reads-{read_sim_mode}.nonconsensus.bam"
  params:
    minimap2 = config["minimap2"],
    samtools = config["samtools"]
  shell: """
    mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/align/read-ONT/locus{wildcards.vntr_idx}/"
    {params.minimap2} -ax map-ont {input.ref} {input.read} | {params.samtools} sort > {output}
    {params.samtools} index {output} 
  """

rule AlignCCSReads:
  input:
    read =lambda wc: sim_path + f"/{wc.gt_sim_mode}/set{wc.genomes_index}/locus{wc.vntr_idx}/reads.locus{wc.vntr_idx}{wc.read_sim_mode}.ccs_0001.fastq",
    ref = sim_path + "/{gt_sim_mode}/set{genomes_index}/locus{vntr_idx}/vntr.locus{vntr_idx}{read_sim_mode}.fa"
  output:
    "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/align/read-CCS/locus{vntr_idx}/reads-{read_sim_mode}.nonconsensus.bam"
  params:
    minimap2 = config["minimap2"],
    samtools = config["samtools"]
  shell: """
    mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/align/read-CCS/locus{wildcards.vntr_idx}/"
    {params.minimap2} -ax map-pb {input.ref} {input.read} | {params.samtools} sort > {output}
    {params.samtools} index {output} 
  """

rule AlignCLRReads:
  input:
    read = lambda wc: sim_path + \
        f"/{wc.gt_sim_mode}/set{wc.genomes_index}/locus{wc.vntr_idx}/reads.locus{wc.vntr_idx}{wc.read_sim_mode}.clr_0001.fastq",
    ref = sim_path + "/{gt_sim_mode}/set{genomes_index}/locus{vntr_idx}/vntr.locus{vntr_idx}{read_sim_mode}.fa"
  output:
    "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/align/read-CLR/locus{vntr_idx}/reads-{read_sim_mode}.nonconsensus.bam"
  params:
    minimap2 = config["minimap2"],
    samtools = config["samtools"]
  shell: """
    mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/align/read-CLR/locus{wildcards.vntr_idx}/"
    {params.minimap2} -ax map-pb {input.ref} {input.read} | {params.samtools} sort > {output}
    {params.samtools} index {output} 
  """

rule AssignMotifsToVNTR:
  input:
    emotifs = emotif_path + "/{mode}/out-{genomes_index}-delta-{delta_threshold}/emotifs.tsv",
    aggregate =  emotif_path + "/{mode}/aggregate.{genomes_index}.tsv"
  output:
    "{input_path}/{mode}/result-{genomes_index}/emotif_{delta_threshold}/motif-{vntr_idx}.csv"
  shell:"""
    mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/emotif_{wildcards.delta_threshold}/"
    awk -v v={wildcards.vntr_idx} 'BEGIN{{OFS="\\t"}}{{ if (NR==FNR) {{ if (FNR==(v+1)) {{coord[FNR] = $1;}} }} else {{ if ($1==coord[(v+1)]) {{print$7; exit;}} }} }}' {input.aggregate} {input.emotifs} > {output}
  """ 

###############################
########### VAMOS #############
###############################

rule RunVamos_nonconsensus:
  input:
    bam = "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/align/read-{read_type}/locus{vntr_idx}/reads-{read_sim_mode}.nonconsensus.bam",
    emotif = "{input_path}/{mode}/result-{genomes_index}/emotif_{delta_threshold}/motif-{vntr_idx}.csv",
    vntr = sim_path + "/{gt_sim_mode}/set{genomes_index}/locus{vntr_idx}/vntr.locus{vntr_idx}{read_sim_mode}.bed"
  output:
    vcf = "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/vamos/read-{read_type}/locus{vntr_idx}/vamos.anno.{read_sim_mode}.delta-{delta_threshold}.nonconsensus.vcf",
  params:
    vamos = config["vamos"],
    grid_opts = config["grid_vamos"]
  shell:"""
    {params.vamos} --raw_anno -i {input.bam} -v {input.vntr} -m {input.emotif} -s sample -o {output.vcf} 
  """

rule RemoveVCFHeader:
  input:
    "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/vamos/read-{read_type}/locus{vntr_idx}/vamos.anno.{read_sim_mode}.delta-{delta_threshold}.{consensus}.vcf"
  output:
    tmp = "{input_path}/{mode}/result-{genomes_index}/{gt_sim_mode}/vamos/read-{read_type}/locus{vntr_idx}/tmp/vamos.anno.{read_sim_mode}.delta-{delta_threshold}.noheader.{consensus}.vcf"
  shell:"""
    mkdir -p "{wildcards.input_path}/{wildcards.mode}/result-{wildcards.genomes_index}/{wildcards.gt_sim_mode}/vamos/read-{wildcards.read_type}/locus{wildcards.vntr_idx}/tmp/"
    awk 'BEGIN{{OFS="\\t"}}{{ if ($1!~/^#/) {{print$0;}} }}' {input} > {output.tmp}
  """
