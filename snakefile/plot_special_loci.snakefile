import os
import tempfile
import subprocess
import os.path
import re
from collections import defaultdict

path = "/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/special_loci"
compressions = ["original", "q-0.1", "q-0.2", "q-0.3"]


annotation_vcf_path = defaultdict(str)
for compression in compressions:
  annotation_vcf_path[compression] = f"/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/sv-vntr/sv/combineAll/individual/annotation/{compression}"


vntr_motif_file = defaultdict(str)
for compression in compressions:
  vntr_motif_file[compression] = f"/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/vntr_motif_database/out-148/{compression}_adjusted/vntr_motifs.bed"


assembly_files = [file.rstrip(".anno.vcf") for file in os.listdir(annotation_vcf_path["original"]) if file.endswith(".anno.vcf")]


asm_path = "/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/sv-vntr/sv/combineAll/individual/annotation"


rule all:
    input:
        ACAN_asm_anno = [f"{path}/ACAN/{compression}/{asm}.anno.vcf" \
                              for compression in compressions for asm in assembly_files],
        WDR7_asm_anno = [f"{path}/WDR7/{compression}/{asm}.anno.vcf" \
                              for compression in compressions for asm in assembly_files],
        CACNA1C_asm_anno = [f"{path}/CACNA1C/{compression}/{asm}.anno.vcf" \
                              for compression in compressions for asm in assembly_files],
        ACAN_bed =  [f"{path}/ACAN/{compression}/vntr_motif.bed" for compression in compressions],
        WDR7_bed = [f"{path}/WDR7/{compression}/vntr_motif.bed" for compression in compressions],
        CACNA1C_bed = [f"{path}/CACNA1C/{compression}/vntr_motif.bed" for compression in compressions],


#######################################
## Select ACAN
#######################################


rule selectACAN:
  input:
    original_asm_anno = asm_path + "/{compression}/{assembly_file}.anno.vcf",
  output:
    asm_anno = path + "/ACAN/{compression}/{assembly_file}.anno.vcf",
  shell:"""
    mkdir -p f"{path}/ACAN/{wildcards.compression}"
    awk 'BEGIN{{OFS="\\t";}}{{if ($1=="chr15" && $2=="88855422" && $8~/88857301/) {{print $0;}} }}' {input.original_asm_anno} > {output.asm_anno}
  """

rule selectACAN_bed:
  input:
    lambda wc: f"{vntr_motif_file[wc.compression]}",
  output:
    path + "/ACAN/{compression}/vntr_motif.bed"
  shell:"""
    mkdir -p f"{path}/ACAN/{wildcards.compression}"
    grep -w 88855422 {input}  | grep -w 88857301 | grep -w chr15 > {output}
  """

#######################################
## Select WDR7
#######################################

rule selectWDR7:
  input:
    original_asm_anno = asm_path + "/{compression}/{assembly_file}.anno.vcf",
  output:
    asm_anno = path + "/WDR7/{compression}/{assembly_file}.anno.vcf",
  shell:"""
    mkdir -p f"{path}/WDR7/{wildcards.compression}"
    awk 'BEGIN{{OFS="\\t";}}{{if ($1=="chr18" && $2=="57024494" && $8~/57024955/) {{print $0;}} }}' {input.original_asm_anno} > {output.asm_anno}
  """

rule selectWDR7_bed:
  input:
    lambda wc: f"{vntr_motif_file[wc.compression]}",
  output:
    path + "/WDR7/{compression}/vntr_motif.bed"
  shell:"""
    mkdir -p f"{path}/WDR7/{wildcards.compression}"
    grep -w 57024494 {input}  | grep -w 57024955 | grep -w chr18 > {output}
  """


#######################################
## Select CACNA1C
#######################################

rule selectCACNA1C:
  input:
    original_asm_anno = asm_path + "/{compression}/{assembly_file}.anno.vcf",
  output:
    asm_anno = path + "/CACNA1C/{compression}/{assembly_file}.anno.vcf",
  shell:"""
    mkdir -p f"{path}/CACNA1C/{wildcards.compression}"
    awk 'BEGIN{{OFS="\\t";}}{{if ($1=="chr12" && $2=="2255790" && $8~/2256090/) {{print $0;}} }}' {input.original_asm_anno} > {output.asm_anno}
  """

rule selectCACNA1C_bed:
  input:
    lambda wc: f"{vntr_motif_file[wc.compression]}",
  output:
    path + "/CACNA1C/{compression}/vntr_motif.bed"
  shell:"""
    mkdir -p f"{path}/CACNA1C/{wildcards.compression}"
    grep -w 2255790 {input}  | grep -w 2256090 | grep -w chr12 > {output}
  """






