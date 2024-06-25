import os
import tempfile
import subprocess
import os.path
import re
from collections import defaultdict

"""
This script separates SVs by chromosomes
And thread through chromosomes
"""


# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
configfile: SD + "/sv-vntr.config"

path = config["path"]

asm_path = f"/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/sv-vntr/sv/combineAll"
vntr_loci_path = config["vntr_loci_path"]

chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6",\
          "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", \
          "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", \
          "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]

asms = [file for file in os.listdir("/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/sv-vntr/sv/combineAll/individual/vcf/indel20") \
      if file.endswith('_hprc.vcf') or file.endswith('_hgsvc.vcf')]

rule all:
    input:
        asm_sv_vcf_chr = [f"{asm_path}/individual/vcf/indel20/{chrom}/{asm}" for asm in asms for chrom in chroms],
        jasmine = [f"{asm_path}/individual/vcf/indel20/{chrom}/combineAll.jasmine.merged.vcf" for chrom in chroms],
        noheaderVCF = [f"{asm_path}/individual/vcf/indel20/{chrom}/combineAll.jasmine.merged.without_header.vcf" for chrom in chroms],
        combinedVCF = f"{asm_path}/individual/vcf/indel20/combineAll.jasmine.merged.vcf",
        finalVCF = f"{asm_path}/individual/vcf/indel20/combineAll.final.vcf"


rule separateVCFByChr:
  input: 
    asm_sv_vcf = "{asm_path}/individual/vcf/indel20/{asm}"
  output:
    asm_sv_vcf_chr = "{asm_path}/individual/vcf/indel20/{chrom}/{asm}"
  shell:"""
    mkdir -p "{wildcards.asm_path}/individual/vcf/indel20/{wildcards.chrom}"
    awk -v chr={wildcards.chrom} 'BEGIN{{OFS="\\t";}}{{if ($1~/^#/) {{print $0;}} else {{if ($1==chr) {{print$0;}}}} }}' {input} > {output}
  """

rule JasmineMergeSV_indel20:
  input:
    "{asm_path}/individual/vcf/indel20/{chrom}/filelist.txt"
  output:
    "{asm_path}/individual/vcf/indel20/{chrom}/combineAll.jasmine.merged.vcf"
  params:
    jasmine = "~/.conda/envs/vamos/bin/jasmine",
    grid_opts = "sbatch -c 2 --mem=20G --time=24:00:00 -p qcb --output=/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/sv-vntr/sv/combineAll/individual/vcf/indel20/submit/log/%j.out",
  shell:"""
    cd "{wildcards.asm_path}/individual/vcf/indel20/{wildcards.chrom}"
    time {params.jasmine}  file_list={input} out_file={output} threads=2 --output_genotypes
    echo "finish step 1"

    time {params.jasmine}  --dup_to_ins --postprocess_only out_file={output} threads=2
    echo "finish step 2"
  """


rule removeHeader:
  input:
    "{asm_path}/individual/vcf/indel20/{chrom}/combineAll.jasmine.merged.vcf"
  output:
    "{asm_path}/individual/vcf/indel20/{chrom}/combineAll.jasmine.merged.without_header.vcf"
  shell:"""
    awk 'BEGIN{{OFS="\\t";}}{{if ($1!~/^#/) {{print $0;}} }}' {input} > {output}
  """

rule combineVCF:
  input:
    vcf = [f"{asm_path}/individual/vcf/indel20/{chrom}/combineAll.jasmine.merged.without_header.vcf" for chrom in chroms],
    header = "/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/sv-vntr/sv/combineAll/individual/vcf/header.txt"
  output:
    "{asm_path}/individual/vcf/indel20/combineAll.jasmine.merged.vcf"
  shell:"""
    cat {input.header} {input.vcf} > {output}
  """

rule postprocess_jasmine:
  input:
    vcf = "{asm_path}/individual/vcf/indel20/combineAll.jasmine.merged.vcf",
    filelist = "{asm_path}/individual/vcf/indel20/filelist.txt",
  output:
    "{asm_path}/individual/vcf/indel20/combineAll.final.vcf"
  params:
    python3 = "~/.conda/envs/lra/bin/python3",
    process_jasmine_out_combineAll_py = config["process_jasmine_out_combineAll_py"]
  shell:"""
    {params.python3} {params.process_jasmine_out_combineAll_py} \
    -v {input.vcf} -f {input.filelist} -o {output}
  """