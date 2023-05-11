import os
import tempfile
import subprocess
import os.path
import re
from collections import defaultdict

"""
This script count SV allele intersecting with VNTRs
Write comparison summary
"""

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
configfile: SD + "configs/sv-vntr.config"

path = config["path"]
compressions = ["original", "q-0.1", "q-0.2", "q-0.3"]
datasets = ["combineAll"]
modes = ["individual"]
# svtypes = ["indel30", "indel20", "indel10"]
svtypes = ["indel20"]

asm_path = defaultdict(str)
for dataset in datasets:
  asm_path[dataset] = f"/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/sv-vntr/sv/{dataset}"

vntr_loci_path = config["vntr_loci_path"]

rule all:
    input:
        vntr_bed = f"{path}/vntr/vntr.loci.bed",
        vntr_ovp = [f"{path}/ovp/{dataset}/{mode}/vcf/{svtype}/{dataset}.vntr.intersect.bed" \
                    for dataset in datasets for mode in modes for svtype in svtypes], 
        allalleles = [f"{path}/vamos/{dataset}/{mode}/{compression}/{svtype}/{compression}.all.bed" \
                    for compression in compressions for dataset in datasets for mode in modes for svtype in svtypes],
        highvamosanno_lowsvalleles = [f"{path}/vamos/{dataset}/{mode}/{compression}/{svtype}/{compression}.highvamosanno_lowsvalleles.bed" \
                    for compression in compressions for dataset in datasets for mode in modes for svtype in svtypes],
        lowvamosanno_highsvalleles = [f"{path}/vamos/{dataset}/{mode}/{compression}/{svtype}/{compression}.lowvamosanno_highsvalleles.bed" \
                    for compression in compressions for dataset in datasets for mode in modes for svtype in svtypes],
        sv_vntr_allele_summary = [f"{path}/summary/{dataset}/{svtype}/sv_vntr_allele_summary.txt" \
                    for dataset in datasets for svtype in svtypes],




        # sv_allele = [f"{path}/sv/{dataset}/{mode}/{dataset}.sv_alleles.bed" for compression in compressions for dataset in datasets for mode in modes], 
        # vamos_allele = [f"{path}/vamos/{dataset}/{mode}/{compression}/{compression}.vamos_alleles.bed" for compression in compressions for dataset in datasets for mode in modes],         
        # png = [f"{path}/vamos/{dataset}/{mode}/{compression}/{compression}.highvamosanno_lowsvalleles.dotplot.png" for compression in compressions for mode in modes]


#######################################
## Select SVs
## 1. Remove SNP and small indels
## 2. find SVs intersect with VNTRs
#######################################

rule extractVNTRloci:
  input:
    vntr_loci_path
  output:
    path + "/vntr/vntr.loci.bed"
  shell:"""
    mkdir -p f"{path}/vntr"
    awk '$1!~/^#/' {input} | awk 'BEGIN{{OFS="\\t";}}{{print $1, $2, $3;}}' > {output}
  """

rule extractSVcoordinates:
  input:
    path + "/sv/{dataset}/{mode}/vcf/{svtype}/{dataset}.final.vcf"
  output:
    path + "/sv/{dataset}/{mode}/vcf/{svtype}/{dataset}.sv.nosnv.bed"
  params:
    postprocessing_py = config["postprocessing_py"],
    python3 = config["python3"]
  shell: """
  {params.python3} {params.filter_SNP_py} -i {input} -o {output}
  """

"""
output 
VNTR loci, SV loci, SV_type
"""
rule overlap:
  input:
    vntr = path + "/vntr/vntr.loci.bed",
    sv = path + "/sv/{dataset}/{mode}/vcf/{svtype}/{dataset}.sv.nosnv.bed"
  output:
    vntr_ovp = path + "/ovp/{dataset}/{mode}/vcf/{svtype}/{dataset}.vntr.intersect.bed",
    sv_ovp = path + "/ovp/{dataset}/{mode}/vcf/{svtype}/{dataset}.sv.intersect.bed"
  params:
    bedtools = config["bedtools"]
  shell: """
    mkdir -p f"{path}/ovp/{wildcards.dataset}/{wildcards.mode}/vcf/{wildcards.svtype}"
    {params.bedtools} intersect -wo -a {input.vntr} -b {input.sv} > {output.vntr_ovp}
    {params.bedtools} intersect -wo -a {input.sv} -b {input.vntr} > {output.sv_ovp}
  """


rule summarizeVNTR_allele_all_allele:
  input:
    vntr = path + "/vntr/vntr.loci.bed",
    ovp = path + "/ovp/{dataset}/{mode}/vcf/{svtype}/{dataset}.vntr.intersect.bed",
    asm = lambda wc: f"{asm_path[wc.dataset]}/{wc.mode}/annotation/{wc.compression}"
  output:
   path + "/vamos/{dataset}/{mode}/{compression}/{svtype}/{compression}.all.bed"
  params:
    compare_sv_vamos_allele_combineAll_py = config["compare_sv_vamos_allele_combineAll_py"],
    python3 = config["python3"]
  shell:"""
    mkdir -p  f"{path}/vamos/{wildcards.dataset}/{wildcards.mode}/{wildcards.compression}/{wildcards.svtype}"
    {params.python3} {params.compare_sv_vamos_allele_combineAll_py} -v {input.vntr} -a {input.asm} -s {input.ovp} -o {output}
  """


rule summarizeVNTR_allele_high_vamos_allele:
  input:
   path + "/vamos/{dataset}/{mode}/{compression}/{svtype}/{compression}.all.bed"
  output:
   path + "/vamos/{dataset}/{mode}/{compression}/{svtype}/{compression}.highvamosanno_lowsvalleles.bed"
  shell:"""
    awk 'BEGIN{{OFS="\\t";}}{{if ($2 > 0 && $4 > 0 && $4 > $2 && $4 - $2 >= 5) {{print $0;}} }}' {input} > {output}
  """

rule summarizeVNTR_allele_low_vamos_allele:
  input:
   path + "/vamos/{dataset}/{mode}/{compression}/{svtype}/{compression}.all.bed"
  output:
   path + "/vamos/{dataset}/{mode}/{compression}/{svtype}/{compression}.lowvamosanno_highsvalleles.bed"
  shell:"""
    awk 'BEGIN{{OFS="\\t";}}{{if ($2 > 0 && $4 > 0 && $2 > $4 && $2 - $4 >= 5) {{print $0;}} }}' {input} > {output}
  """

"""
Edit $2, $3, $4, $5
 str(numSValleles), str(numSVLenalleles), \
 str(numLenAnno), str(numStringAnno), \
"""
rule writeSummary:
    input:
      sv_nosnp = path + "/sv/{dataset}/{mode}/vcf/{svtype}/{dataset}.final.vcf",
      sv_ovp = path + "/ovp/{dataset}/{mode}/vcf/{svtype}/{dataset}.sv.intersect.bed",
      all_alleles = path + "/vamos/{dataset}/{mode}/{compression}/{svtype}/{compression}.all.bed"
    output: 
      path + "/summary/tmp/{svtype}/{dataset}_{mode}_{compression}.sv_vntr_allele_summary.txt"
    shell:"""
      mkdir -p f"{path}/summary/tmp/{wildcards.svtype}"

      echo "==============================================" >> {output}
      echo "{wildcards.mode}" >> {output}
      echo "==============================================" >> {output}

      echo "{wildcards.dataset} {wildcards.compression}" >> {output}
      awk 'END{{print "Total SVs: ", NR}}' {input.sv_nosnp} >> {output}
      awk 'END{{print "VNTR SVs: ", NR}}' {input.sv_ovp} >> {output}
      awk 'BEGIN{{SUM=0; Nlines=0;}}{{ if ($2 > 0 && $4 > 0) {{SUM+=$2; Nlines+=1;}}  }}END{{print "VNTR SV alleles: ", SUM/Nlines, Nlines;}}' {input.all_alleles} >> {output}
      awk 'BEGIN{{SUM=0; Nlines=0;}}{{ if ($2 > 0 && $4 > 0) {{SUM+=$3; Nlines+=1;}}  }}END{{print "VNTR SV Len alleles: ", SUM/Nlines, Nlines;}}' {input.all_alleles} >> {output}

      echo "{wildcards.dataset} {wildcards.compression}" >> {output}
      awk 'BEGIN{{SUM=0;}}{{if ($2 > 0 && $4 > 0) {{SUM+=$4; Nlines+=1; }} }}END{{print "vamos length alleles: ", SUM/Nlines, Nlines;}}'  {input.all_alleles} >> {output}
      awk 'BEGIN{{SUM=0;}}{{if ($2 > 0 && $4 > 0) {{SUM+=$5; Nlines+=1; }}  }}END{{print "vamos vec alleles: ", SUM/Nlines, Nlines;}}'  {input.all_alleles} >> {output}
      echo "\n\n" >> {output}
    """

rule concatSummary:
  input:
    sorted([f"{path}/summary/tmp/{svtype}/{dataset}_{mode}_{compression}.sv_vntr_allele_summary.txt" for svtype in svtypes for compression in compressions for dataset in datasets for mode in modes])
  output:
    path + "/summary/{dataset}/{svtype}/sv_vntr_allele_summary.txt"
  shell:"""
    mkdir -p f"{path}/summary/{wildcards.dataset}/{wildcards.svtype}"
    cat {input} > {output}
  """
