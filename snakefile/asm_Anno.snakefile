import os
import tempfile
import subprocess
import os.path
import re

"""
This snakefile annotate assemblies using original & efficient motifs
"""

SD = os.path.dirname(workflow.snakefile)

# Config
configfile: "/project/mchaisso_100/cmb-16/jingwenr/trfCall/script/snakefile.json"

input_path = config["input_path"]
assembly = config["assembly"]
aligner = ["lra"]
date = config["date"]
mode = config["mode"]
genomes_index = ["64"]
delta_threshold = ["0.1", "0.2", "0.3", "0.4", "0.5"]

rule all:
    input:
        bam = expand("{input_path}/assembly/align/{assembly}.{aligner}.bam", input_path = input_path, assembly = assembly, aligner = aligner),
        omotif = expand("{input_path}/emotifs/{mode}/out-{genomes_index}-original_motifs/motifs.csv",\
                input_path = input_path, mode = mode, genomes_index = genomes_index),
        vntr = expand("{input_path}/emotifs/{mode}/out-{genomes_index}-original_motifs/vntr.bed",\
                input_path = input_path, mode = mode, genomes_index = genomes_index),
        original_vcf = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/{assembly}.{aligner}.vcf", \
                input_path = input_path, assembly = assembly, aligner = aligner, date = date, mode = mode, genomes_index = genomes_index),
        efficient_vcf = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/{assembly}.{aligner}.vcf", \
                input_path = input_path, assembly = assembly, aligner = aligner, date = date, mode = mode, genomes_index = genomes_index, delta_threshold = delta_threshold),

rule align_lra:
    input:
        "{input_path}/hgsvc/{assembly}"
    output:
        "{input_path}/assembly/align/{assembly}.lra.bam"
    params:
        lra = config["lra"],
        ref = config["ref"],
        samtools = config["samtools"],
        grid_opts = config["grid_align"]
    shell:"""
        mkdir -p "{wildcards.input_path}/assembly/align/"
        {params.lra} align -CONTIG {params.ref} {input} -t 16 -p s -H | {params.samtools} sort -@4 > {output}
        {params.samtools} index -@16 {output} 
    """

rule align_mm2:
    input:
        "{input_path}/hgsvc/{assembly}"
    output:
        "{input_path}/assembly/align/{assembly}.mm2.bam"
    params:
        mm2 = config["mm2"],
        ref_mm2 = config["ref_mm2"],
        samtools = config["samtools"],
        grid_opts = config["grid_align"]
    shell:"""
        mkdir -p "{wildcards.input_path}/assembly/align/"
        {params.mm2} -a {params.ref_mm2} {input} -t 16 | {params.samtools} sort -@4 > {output}
        {params.samtools} index -@16 {output} 
    """

"""
based on original motifs
"""
rule processOmotif:
    input:
        "{input_path}/emotifs/{mode}/aggregate.{genomes_index}.tsv"
    output:
        omotif = "{input_path}/emotifs/{mode}/out-{genomes_index}-original_motifs/motifs.csv",
        vntr = "{input_path}/emotifs/{mode}/out-{genomes_index}-original_motifs/vntr.bed"
    shell:"""
        awk '{{if (NR > 1) {{print $2;}} }}' {input} | awk '{{split($1, a, ","); if (length(a) <= 256) {{print$1;}} }}' > {output.omotif}
        awk '{{if (NR > 1) {{print $1, $2;}} }}'  {input} | awk 'BEGIN{{OFS="\\t";}}{{split($1, b,":|-"); split($2,a,","); if (length(a) <= 256) {{print b[1], b[2], b[3];}} }}' > {output.vntr}
    """

rule RunVamos_original:
    input:
        bam = "{input_path}/assembly/align/{assembly}.{aligner}.bam",
        omotif = "{input_path}/emotifs/{mode}/out-{genomes_index}-original_motifs/motifs.csv",
        vntr = "{input_path}/emotifs/{mode}/out-{genomes_index}-original_motifs/vntr.bed"
    output:
        vcf = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/{assembly}.{aligner}.vcf"
    params:
        vamos = config["vamos"],
        grid_opts = config["grid_vamos"]
    shell:"""
        mkdir -p "{wildcards.input_path}/{wildcards.date}/{wildcards.mode}/result-{wildcards.genomes_index}/assembly/anno-original-motifs/"
        {params.vamos}  -i {input.bam} -v {input.vntr} -m {input.omotif} -s {wildcards.assembly} -o {output.vcf} -t 32
    """

rule RunVamos_efficient:
    input:
        bam = "{input_path}/assembly/align/{assembly}.{aligner}.bam",
        emotif = "{input_path}/emotifs/{mode}/out-{genomes_index}-delta-{delta_threshold}/motifs.csv",
        vntr = "{input_path}/emotifs/{mode}/out-{genomes_index}-delta-{delta_threshold}/vntrs.bed"
    output:
        vcf = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/{assembly}.{aligner}.vcf"
    params:
        vamos = config["vamos"],
        grid_opts = config["grid_vamos"]
    shell:"""
        mkdir -p "{wildcards.input_path}/{wildcards.date}/{wildcards.mode}/result-{wildcards.genomes_index}/assembly/anno-delta-{wildcards.delta_threshold}/"
        {params.vamos} -i {input.bam} -v {input.vntr} -m {input.emotif} -s {wildcards.assembly} -o {output.vcf} -t 32 
    """

