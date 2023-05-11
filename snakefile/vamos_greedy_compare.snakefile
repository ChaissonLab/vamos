import os
import tempfile
import subprocess
import os.path
import re

configfile: "/project/mchaisso_100/cmb-16/jingwenr/trfCall/vamos/snakefile/configs/snakefile.json"

modes = ["q-0.1"]
input_path = "/project/mchaisso_100/cmb-16/jingwenr/trfCall/revision/greedy_comparison"
asms = [asm.rstrip(".anno.vcf") for asm in os.listdir(f"{input_path}/anno/vamos/original") if asm.endswith(".anno.vcf")]
filter_orignal_motifs_threshold = [0, 5, 10, 20, 30, 50]

rule all:
    input:
        summary = expand( "{input_path}/compare_result/{asm}_{mode}_{filter}_compare_stats.bed", \
                        input_path = input_path, asm = asms, mode=modes, filter = filter_orignal_motifs_threshold),
        summary_final = expand("{input_path}/summary/{mode}/{mode}_{filter}_compare_stats_summary.bed", \
            input_path = input_path, mode=modes, filter = filter_orignal_motifs_threshold)


rule RunVamos_Liftover:
  input:
    bam = "/project/mchaisso_100/cmb-16/bidagu/working/vamos/hgsvc_2022-09-01/{mode}/analysis/{asm}/aln/{asm}.aln.bam",
    vntr = "/project/mchaisso_100/cmb-16/jingwenr/trfCall/emotifs/hgsvc_2022-06-25_delta_3/out-64-delta-0.2/{mode}.emotifs.bed"
  output:
    "{input_path}/analysis/{mode}/{asm}/vcf/{asm}.vcf"
  params:
    vamos = config["vamos"],
    grid_opts = "sbatch -c 64 --mem=50G --time=24:00:00 -p qcb"
  shell:"""
    mkdir -p "{wildcards.input_path}/analysis/{wildcards.mode}/{wildcards.asm}/vcf/"
    {params.vamos} --readwise -b {input.bam} -r {input.vntr} -s {wildcards.asm} -o {output} -t 64
  """

rule CompareTwoAnnotation: 
    input: 
        vamos_anno_vcf = "{input_path}/anno/vamos/{mode}/{asm}.anno.vcf",
        greedy_anno_vcf = "{input_path}/anno/greedy/{asm}.anno.vcf",
        liftover_seq = "{input_path}/liftover_asm_seq/{asm}.fa",
        original_vntr_motif = "/project/mchaisso_100/cmb-16/bidagu/pipelines/vamos/pipeline/databases/original_adjusted/original_motifs.set148.tsv"
    output:
        "{input_path}/compare_result/{asm}_{mode}_{filter}_compare_stats.bed"
    params:
        compare_with_greedy_py = config["compare_with_greedy_py"],
        grid_opts = "sbatch -c 1 --mem=10G --time=24:00:00 -p qcb",
        python3 = "~/.conda/envs/lra/bin/python3"
    shell:"""
        mkdir -p "{wildcards.input_path}/compare_result"
        {params.python3} {params.compare_with_greedy_py} \
        -i {wildcards.asm} \
        -o {output} \
        -m {wildcards.mode} \
        -v {input.vamos_anno_vcf} \
        -g {input.greedy_anno_vcf} \
        -l {input.liftover_seq} \
        -f {input.original_vntr_motif} \
        -p {wildcards.filter}
    """

rule combineCompare:
  input:
    expand("{input_path}/compare_result/{asm}_{mode}_{filter}_compare_stats.bed", asm = asms, allow_missing = True)
  output:
    "{input_path}/summary/{mode}/{mode}_{filter}_compare_stats_summary.bed"
  params:
    grid_opts = "sbatch -c 1 --mem=5G --time=24:00:00 -p qcb",
  shell:"""
    mkdir -p "{wildcards.input_path}/summary/{wildcards.mode}"
    cat {input} > {output}
  """