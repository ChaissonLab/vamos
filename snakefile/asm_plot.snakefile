import os
import tempfile
import subprocess
import os.path
import re

SD = os.path.dirname(workflow.snakefile)

# Config
configfile: "/project/mchaisso_100/cmb-16/jingwenr/trfCall/script/snakefile.json"

input_path = config["input_path"]
assembly = config["assembly"]
aligner = ["lra"]
date = config["date"]
mode = config["mode"]
genomes_index = ["64"]
delta_threshold=["0.1", "0.2", "0.3", "0.4", "0.5"]

rule all:
    input:
      efficient_plot_png = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/plot/{aligner}.numAlleles.png", \
                    input_path = input_path, date = date, mode = mode, genomes_index = genomes_index, delta_threshold = delta_threshold, aligner = aligner),
      efficient_hist_png = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/plot/{aligner}.numAlleles_hist.png", \
                    input_path = input_path, date = date, mode = mode, genomes_index = genomes_index, delta_threshold = delta_threshold, aligner = aligner),
      efficient_plot_eps = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/plot/{aligner}.numAlleles.eps", \
                    input_path = input_path, date = date, mode = mode, genomes_index = genomes_index, delta_threshold = delta_threshold, aligner = aligner),
      efficient_hist_eps = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/plot/{aligner}.numAlleles_hist.eps", \
                    input_path = input_path, date = date, mode = mode, genomes_index = genomes_index, delta_threshold = delta_threshold, aligner = aligner),
      original_scatter_png = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/plot/{aligner}.numAlleles.png",\
                    input_path = input_path, date = date, mode = mode, genomes_index = genomes_index, aligner = aligner),
      original_hist_png = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/plot/{aligner}.numAlleles_hist.png",\
                    input_path = input_path, date = date, mode = mode, genomes_index = genomes_index, aligner = aligner),
      original_scatter_eps = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/plot/{aligner}.numAlleles.eps",\
                    input_path = input_path, date = date, mode = mode, genomes_index = genomes_index, aligner = aligner),
      original_hist_eps = expand("{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/plot/{aligner}.numAlleles_hist.eps", \
                    input_path = input_path, date = date, mode = mode, genomes_index = genomes_index, aligner = aligner)

rule plot_efficient:
  output:
    scatter_png = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/plot/{aligner}.numAlleles.png",
    hist_png = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/plot/{aligner}.numAlleles_hist.png",
    scatter_eps = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/plot/{aligner}.numAlleles.eps",
    hist_eps = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-delta-{delta_threshold}/plot/{aligner}.numAlleles_hist.eps"
  params:
    Assembly_alleles_py = config["Assembly_alleles_py"]
  shell:"""
    mkdir -p "{wildcards.input_path}/{wildcards.date}/{wildcards.mode}/result-{wildcards.genomes_index}/assembly/anno-delta-{wildcards.delta_threshold}/plot"
    python3 {params.Assembly_alleles_py} {wildcards.input_path}/{wildcards.date}/{wildcards.mode}/result-{wildcards.genomes_index}/assembly/ anno-delta-{wildcards.delta_threshold} {output.scatter_png} {output.hist_png} {output.scatter_eps} {output.hist_eps} {wildcards.aligner}
  """


rule plot_original:
  output:
    scatter_png = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/plot/{aligner}.numAlleles.png",
    hist_png = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/plot/{aligner}.numAlleles_hist.png",
    scatter_eps = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/plot/{aligner}.numAlleles.eps",
    hist_eps = "{input_path}/{date}/{mode}/result-{genomes_index}/assembly/anno-original-motifs/plot/{aligner}.numAlleles_hist.eps"
  params:
    Assembly_alleles_py = config["Assembly_alleles_py"]
  shell:"""
    mkdir -p "{wildcards.input_path}/{wildcards.date}/{wildcards.mode}/result-{wildcards.genomes_index}/assembly/anno-original-motifs/plot"
    python3 {params.Assembly_alleles_py} {wildcards.input_path}/{wildcards.date}/{wildcards.mode}/result-{wildcards.genomes_index}/assembly/ anno-original-motifs {output.scatter_png} {output.hist_png} {output.scatter_eps} {output.hist_eps} {wildcards.aligner}
  """