import os
import tempfile
import subprocess
import os.path
import re
from collections import defaultdict

"""
Filter SV calls by SNP and indels, need to specify indel length
"""

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
configfile: SD + "/sv-vntr.config"

hgsvc_path = "/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/sv-vntr/sv/hgsvc/individual"
hgsvc_asms = [file.rstrip("_dipcall.dip.vcf") for file in os.listdir(hgsvc_path) if file.endswith('_dipcall.dip.vcf')]

hprc_path = "/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/sv-vntr/sv/hprc/individual"
hprc_asms = [file.rstrip("_hprc.original.vcf") for file in os.listdir(hprc_path) if file.endswith('_hprc.original.vcf')]

rule all:
    input:
        sv_nosnp_hgsvc = [f"{hgsvc_path}/{hgsvc_asm}_hgsvc.vcf" for hgsvc_asm in hgsvc_asms],
        sv_nosnp_hprc = [f"{hprc_path}/{hprc_asm}_hprc.vcf" for hprc_asm in hprc_asms],

rule filterSNP_Indels_hgsvc:
  input:
    hgsvc_path + "/{hgsvc_asm}_dipcall.dip.vcf"
  output:
    hgsvc_path + "/{hgsvc_asm}_hgsvc.vcf"
  params:
    filter_SNP_Indels_py = config["filter_SNP_Indels_py"],
    python3 = config["python3"],
    grid_opts = config["grid_small"]
  shell: """
  {params.python3} {params.filter_SNP_Indels_py} -i {input} -o {output} -a "{wildcards.hgsvc_asm}_hgsvc" -l 30 
  """


rule filterSNP_Indels_hprc:
  input:
    hprc_path + "/{hprc_asm}_hprc.original.vcf"
  output:
    hprc_path + "/{hprc_asm}_hprc.vcf"
  params:
    filter_SNP_Indels_py = config["filter_SNP_Indels_py"],
    python3 = config["python3"],
    grid_opts = config["grid_small"]
  shell: """
  {params.python3} {params.filter_SNP_Indels_py} -i {input} -o {output} -a "{wildcards.hprc_asm}_hprc" -l 30
  """


