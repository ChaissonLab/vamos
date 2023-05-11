import os
import tempfile
import subprocess
import os.path
import re
from collections import defaultdict

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
configfile: SD + "/configs/novel_vntr.config"

path = config["path"]
datasets = ["hprc"]

asms = defaultdict(list)
asms_paths = defaultdict(str)
asms_fa = defaultdict(lambda: defaultdict(str))
vntr_loci = defaultdict(list)

for dataset in datasets:
  asms_paths[dataset] = f"{path}/{dataset}/trf"
  asms[dataset] = [asm.rstrip(".trf.bed") for asm in os.listdir(asms_paths[dataset]) if asm.endswith(".trf.bed")]

for dataset in datasets:
  for asm in asms[dataset]:
    sample_id, hap = asm.split('_')
    if hap == "paternal":
      asms_fa[dataset][asm] = f'/project/mchaisso_100/cmb-16/bidagu/databases/{dataset}-v2/fasta/{sample_id}.{hap}.1.fa'
    else:
      asms_fa[dataset][asm] = f'/project/mchaisso_100/cmb-16/bidagu/databases/{dataset}-v2/fasta/{sample_id}.{hap}.0.fa'

# for dataset in datasets:
  # vntr_loci[dataset] = [vntr.rstrip('.bed') for vntr in os.listdir(f'{path}/{dataset}/asm/vntr_bed/bed') if vntr.endswith(".bed")]

rule all:
    input:
      vntr_fa = expand(path + "/{dataset}/asm/asm_vntr_fa/{asm}.fa", dataset=datasets, asm=asms[dataset])

rule extract_vntr_seq:
  input:
    fa = lambda wc: asms_fa[wc.dataset][wc.asm],
    bed = path + "/{dataset}/asm/asm_vntr_bed/bed/{asm}.bed"
  output:
    vntr_fa = path + "/{dataset}/asm/asm_vntr_fa/{asm}.fa"
  params:
    grid_opts = config["grid_small"],
    samtools = config["samtools"]
  shell:"""
    mkdir -p "{path}/{dataset}/asm/asm_vntr_fa"
    {params.samtools} faidx {input.fa} -r {input.bed} > {output}
  """

