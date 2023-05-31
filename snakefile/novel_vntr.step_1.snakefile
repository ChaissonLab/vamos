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
datasets = ["hgsvc", "hprc"]

asms = defaultdict(list)
asms_paths = defaultdict(str)
asms_fa = defaultdict(lambda: defaultdict(str))

for dataset in datasets:
  asms_paths[dataset] = f"{path}/{dataset}/trf"
  asms[dataset] = [asm.rstrip(".trf.bed") for asm in os.listdir(asms_paths[dataset]) if asm.endswith(".trf.bed")]

for dataset in datasets:
  for asm in asms[dataset]:
    sample_id, hap = asm.split('_')
    if hap == "paternal":
      asms_fa[dataset][asm] = f'{path}/{dataset}/fasta/{sample_id}.{hap}.1.fa'
    else:
      asms_fa[dataset][asm] = f'{path}/{dataset}/fasta/{sample_id}.{hap}.0.fa'

rule all:
    input:
      filtered_trf_intv = [f"{path}/{dataset}/filter_trf/{asm}.trf.filter.bed" for dataset in asms for asm in asms[dataset]],
      merged_trf_intv = [f"{path}/{dataset}/merged_trf/{asm}.trf.filter.bed" for dataset in asms for asm in asms[dataset]],
      flank_coordinate = [f"{path}/{dataset}/flanked_bed/{asm}.flank.bed" for dataset in asms for asm in asms[dataset]],
      flank_seq = [f"{path}/{dataset}/flanked_seq/{asm}.flank.fa" for dataset in asms for asm in asms[dataset]],
      align_bam = [f"{path}/{dataset}/align_bam/{asm}.flank.bam" for dataset in asms for asm in asms[dataset]],
      align_bam_filter = [f"{path}/{dataset}/align_bam_convert_coordinate/{asm}.flank.bed" for dataset in asms for asm in asms[dataset]],
      asm_mapped_region = [f"{path}/{dataset}/mapped_bed/{asm}.map.bed" for dataset in asms for asm in asms[dataset]],
      asm_filed_bed = [f"{path}/{dataset}/mapped_filtered_bed/{asm}.map.filtered.bed" for dataset in asms for asm in asms[dataset]],
      combined_bed = [f"{path}/{dataset}/mapped_filtered_bed/combine.map.filtered.bed" for dataset in asms],
      ref_merged_bed = [f"{path}/{dataset}/ref_bed/ref.merged.bed" for dataset in asms],
      ref_merged_filter_bed = [f"{path}/{dataset}/ref_bed/ref.merged.repeat_filter.center_filter.bed" for dataset in asms],
      ref_merged_filter_bed_combined = f"{path}/ref_bed/ref.merged.repeat_filter.center_filter.merged_dataset.bed",
      ref_merged_filter_bed_combined_filter_len = f"{path}/ref_bed/ref.merged.filtered.repeat_filter.center_filter.merged_dataset.len_filter.bed",
      asm_bed = [f"{path}/{dataset}/asm/asm_vntr_bed/bed/{asm}.bed" for dataset in asms for asm in asms[dataset]],
      asm_vntr_fa = [f"{path}/{dataset}/asm/asm_vntr_fa/{asm}.fa" for dataset in asms for asm in asms[dataset]],
      asm_vntr_bida_fa = [f"{path}/{dataset}/asm/asm_vntr_bida_fa/{asm}.fa" for dataset in asms for asm in asms[dataset]]
      # output = expand(path + "/{dataset}/ref_seq/log.txt", dataset = datasets)


rule filter_trf_intv:
  input:
    lambda wc: asms_paths[wc.dataset] + "/{asm}.trf.bed"
  output:
    path + "/{dataset}/filter_trf/{asm}.trf.filter.bed"
  params:
    grid_opts = config["grid_small"]
  shell:"""
    awk 'BEGIN{{OFS="\\t";}}{{if (length($10) >= 6) {{print $0;}} }}' {input} | sort -k1,1 -k2,2n -k3,3n > {output}
  """

rule merge_trf_intv:
  input:
    path + "/{dataset}/filter_trf/{asm}.trf.filter.bed"
  output:
    path + "/{dataset}/merged_trf/{asm}.trf.filter.bed",
  params:
    grid_opts = config["grid_small"],
    python3 = config["python3"],
    merge_intv_py = config["merge_intv_py"]
  shell:"""
    mkdir -p "{path}/{wildcards.dataset}/merged_trf/"
    {params.python3} {params.merge_intv_py} -i {input} -o {output} 
  """

rule extract_flanking_coordinate:
  input:
    intv = path + "/{dataset}/merged_trf/{asm}.trf.filter.bed",
    fa = lambda wc: asms_fa[wc.dataset][wc.asm]
  output:
    path + "/{dataset}/flanked_bed/{asm}.flank.bed"
  params:
    grid_opts = config["grid_small"],
    python3 = config["python3"],
    flank_coordinate_py = config["flank_coordinate_py"]
  shell:"""
  {params.python3} {params.flank_coordinate_py} -i {input.intv} -o {output} -f {input.fa}
  """

rule extract_flanking_seq:
  input:
    region = path + "/{dataset}/flanked_bed/{asm}.flank.bed",
    fa = lambda wc: asms_fa[wc.dataset][wc.asm]
  output:
    path + "/{dataset}/flanked_seq/{asm}.flank.fa"
  params:
    grid_opts = config["grid_small"],
    samtools = config["samtools"]
  shell:"""
  {params.samtools} faidx {input.fa} -r {input.region} > {output}
  """

rule align_flank_seq:
  input:
    path + "/{dataset}/flanked_seq/{asm}.flank.fa"
  output:
    path + "/{dataset}/align_bam/{asm}.flank.bam"
  params:
    mm2 = config["mm2"],
    ref_mm2 = config["ref_mm2"],
    samtools = config["samtools"],
    grid_opts = config["grid_align"]
  shell:"""
    mkdir -p "{path}/{dataset}/align_bam"
    {params.mm2} -a {params.ref_mm2} {input} -t 16 | {params.samtools} sort -@4 > {output}
    {params.samtools} index -@16 {output} 
  """

# filter unmapped, not primary, duplicate, low mapq < 20
rule filter_align_bam:
  input:
    path + "/{dataset}/align_bam/{asm}.flank.bam"
  output:
    path + "/{dataset}/align_bam_convert_coordinate/{asm}.flank.bed"
  params:
    samtools = config["samtools"],
    grid_opts = config["grid_median"]
  shell:"""
    mkdir -p "{path}/{wildcards.dataset}/align_bam_convert_coordinate"
    {params.samtools} view -F 1284 -q 20 -@12 {input} | awk 'BEGIN{{OFS="\\t";}}{{print $1, $2, $3, $4, $4+5000, $5, $6;}}' > {output}
  """

rule validate_mapped_trf_intv:
  input:
    ref_region = path + "/{dataset}/align_bam_convert_coordinate/{asm}.flank.bed",
    asm_region = path + "/{dataset}/flanked_bed/{asm}.flank.bed"
  output:
    path + "/{dataset}/mapped_bed/{asm}.map.bed"
  params:
    asm_trf_map_bed_py = config["asm_trf_map_bed_py"],
    python3 = config["python3"]
  shell:"""
    mkdir -p "{path}/{wildcards.dataset}/mapped_bed/"
    {params.python3} {params.asm_trf_map_bed_py} -r {input.ref_region} -a {input.asm_region} -o {output} -t {wildcards.dataset} -p {wildcards.asm}
  """

rule filter_exisiting_vntr:
  input:
    asm_bed = path + "/{dataset}/mapped_bed/{asm}.map.bed", 
    vntr_bed = "/project/mchaisso_100/mchaisso/projects/vamos/HGSVC_SV_Callsets/tandem.merged.bed"
  output:
    asm_bed = path + "/{dataset}/mapped_filtered_bed/{asm}.map.filtered.bed"
  params:
    bedtools = config["bedtools"]
  shell:"""
    mkdir -p "{path}/{wildcards.dataset}/mapped_filtered_bed"
    {params.bedtools} intersect -a {input.asm_bed} -b {input.vntr_bed} -wao | \
    awk 'BEGIN{{OFS="\\t"}}{{ if ($6 ==".") {{print $1, $2, $3, $4, $5;}} }}' - > {output}
  """

rule combineBed:
  input:
    lambda wc: expand(path + "/{dataset}/mapped_filtered_bed/{asm}.map.filtered.bed", asm=asms[wc.dataset], allow_missing=True)
  output:
    bed = path + "/{dataset}/mapped_filtered_bed/combine.map.filtered.bed",
    srt_bed = path + "/{dataset}/mapped_filtered_bed/combine.map.filtered.srt.bed"
  shell:"""
    [ -e file ] && rm {output}
    cat {input} >> {output.bed}
    sort -k1,1 -k2,2n -k3,3n {output.bed} > {output.srt_bed}
  """

rule merge_ref_intv:
  input:
    path + "/{dataset}/mapped_filtered_bed/combine.map.filtered.srt.bed"
  output:
    path + "/{dataset}/ref_bed/ref.merged.bed"
  params:
    grid_opts = config["grid_small"],
    python3 = config["python3"],
    merge_ref_intv_py = config["merge_ref_intv_py"]
  shell:"""
    mkdir -p "{path}/{wildcards.dataset}/ref_bed/"
    {params.python3} {params.merge_ref_intv_py} -i {input} -o {output} 
  """

rule filter_repeat_element:
  input:
    before = path + "/{dataset}/ref_bed/ref.merged.bed",
    repeat = "/project/mchaisso_100/shared/references/hg38/regions/RepeatMasker/RepeatMasker.bed"
  params:
    bedtools = config["bedtools"]
  output:
     path + "/{dataset}/ref_bed/ref.merged.repeat_filter.bed",
  shell:"""
    {params.bedtools} intersect -a {input.before} -b {input.repeat} -wao | \
    awk 'BEGIN{{OFS="\\t"}}{{ if ($10 == ".") {{print $1, $2, $3, $4, $5, $6, $7, $8, $9;}} }}' - > {output}
  """

rule filter_centermere:
  input:
    before = path + "/{dataset}/ref_bed/ref.merged.repeat_filter.bed",
    centermere = "/project/mchaisso_100/cmb-16/bidagu/databases/masks/para_centro_grch38_from_quentin.bed"
  params:
    bedtools = config["bedtools"]
  output:
     path + "/{dataset}/ref_bed/ref.merged.repeat_filter.center_filter.bed",
  shell:"""
    {params.bedtools} intersect -a {input.before} -b {input.centermere} -wao | \
    awk -v d={wildcards.dataset} 'BEGIN{{OFS="\\t"}}{{ if ($10 == ".") {{print $1, $2, $3, $4, $5, $6, $7, $8, $9, d;}} }}' - > {output}
  """

rule combine_ref_bed_dataset:
  input:
    [f"{path}/{dataset}/ref_bed/ref.merged.repeat_filter.center_filter.bed" for dataset in datasets]
  output:
    f"{path}/ref_bed/ref.merged.repeat_filter.center_filter.bed"
  shell:"""
    mkdir -p "{path}/ref_bed"
    cat {input} | sort -k1,1 -k2,2n -k3,3n - >> {output}
  """

rule merge_intv_dataset:
  input:
    f"{path}/ref_bed/ref.merged.repeat_filter.center_filter.bed"
  output:
    f"{path}/ref_bed/ref.merged.repeat_filter.center_filter.merged_dataset.bed"
  params:
    python3 = config["python3"],
    grid_opts = config["grid_small"],
    merge_ref_intv_dataset_py = config["merge_ref_intv_dataset_py"]  
  shell:"""
    {params.python3} {params.merge_ref_intv_dataset_py} -i {input} -o {output} 
  """

rule filter_by_len_comparison:
  input:
    f"{path}/ref_bed/ref.merged.repeat_filter.center_filter.merged_dataset.bed"
  output:
    f"{path}/ref_bed/ref.merged.filtered.repeat_filter.center_filter.merged_dataset.len_filter.bed"
  params:
    python3 = config["python3"],
    grid_opts = config["grid_small"],
    intv_length_filter_py = config["intv_length_filter_py"]  
  shell:"""
    {params.python3} {params.intv_length_filter_py} -i {input} -o {output} 
  """

rule extract_intv_asm_and_ref:
  input: 
    bed = f"{path}/ref_bed/ref.merged.filtered.repeat_filter.center_filter.merged_dataset.len_filter.bed",
    fa = lambda wc: asms_fa[wc.dataset][wc.asm]
  output:
    bed = path + "/{dataset}/asm/asm_vntr_bed/bed/{asm}.bed",
    full_bed = path + "/{dataset}/asm/asm_vntr_bed/full_bed/{asm}.full.bed"
  params:
    python3 = config["python3"],
    grid_opts = config["grid_small"],
    extract_intv_asm_and_ref_py = config["extract_intv_asm_and_ref_py"]
  shell:"""
    mkdir -p "{path}/{wildcards.dataset}/asm/asm_vntr_bed"
    mkdir -p "{path}/{wildcards.dataset}/asm/asm_vntr_bed/full_bed"
    mkdir -p "{path}/{wildcards.dataset}/asm/asm_vntr_bed/bed"

    [ -e file ] && rm {output}

    {params.python3} {params.extract_intv_asm_and_ref_py} -i {input.bed} -o {output.bed} -f {output.full_bed} -a {wildcards.asm} -d {wildcards.dataset}
  """

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

    if [ -s {input.bed} ]; then
        {params.samtools} faidx {input.fa} -r {input.bed} > {output}
    else
        touch {output}
    fi
  """

rule arrange_format_for_Bida:
  input:
    input_fa = path + "/{dataset}/asm/asm_vntr_fa/{asm}.fa",
    input_bed = path + "/{dataset}/asm/asm_vntr_bed/full_bed/{asm}.full.bed"
  output:
    output = path + "/{dataset}/asm/asm_vntr_bida_fa/{asm}.fa"
  params:
    python3 = config["python3"],
    grid_opts = config["grid_small"],
    arrange_format_for_Bida_py = config["arrange_format_for_Bida_py"]
  shell:"""
    mkdir -p "{path}/{wildcards.dataset}/asm/asm_vntr_bida_fa"
    if [ -s {input.input_bed} ]; then
        {params.python3} {params.arrange_format_for_Bida_py} -f {input.input_fa} -b {input.input_bed} -o {output}
    else
        touch {output}
    fi
  """

# rule arrange_ref_asm_fa:
#   input:
#     input_fa_dir = path + "/{dataset}/asm/asm_vntr_fa",
#     input_bed_dir = path + "/{dataset}/asm/asm_vntr_bed/full_bed"
#   output:
#     output = path + "/{dataset}/ref_seq/log.txt"
#   params:
#     python3 = config["python3"],
#     grid_opts = config["grid_small"],
#     arrange_ref_asm_fa_py = config["arrange_ref_asm_fa_py"]
#   shell:"""
#     {params.python3} {params.arrange_ref_asm_fa_py} -i {input.input_fa_dir} -b {input.input_bed_dir} -o {path}/{wildcards.dataset}/ref_seq > {output}
#   """
