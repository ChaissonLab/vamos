#!~/.conda/envs/lra/bin/python3
import os
import sys
import re
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt
from collections import defaultdict 
import numpy as np
import UnionFind
import edlib

os.chdir(r"/project/mchaisso_100/cmb-16/jingwenr/trfCall")

if len(sys.argv) == 8:
    input_path = sys.argv[1]
    branch = sys.argv[2]
    out_scatter_png = sys.argv[3]
    out_hist_png = sys.argv[4]
    out_scatter_eps = sys.argv[5]
    out_hist_eps = sys.argv[6]
    aligner = sys.argv[7]
else:
    print("require input_path, branch_assembly,out_scatter_png, out_hist_png, out_scatter_eps, out_hist_eps, aligner")
    exit(1)

assembly = ["v12_HG00512_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA20509_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03065_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG02011_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA19650_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG03486_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_HG00732_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA18534_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG00096_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_NA19983_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA18939_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG00731_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA19238_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG02818_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG01596_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG02587_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00513_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_HG01114_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03371_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG00864_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03683_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03009_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_NA20847_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA12329_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03732_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG01505_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00171_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03125_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA19239_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_HG02492_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03009_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03683_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG00864_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA20847_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG01596_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03371_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG01114_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00513_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG02587_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03125_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG00171_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG02492_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA19239_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA12329_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG01505_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03732_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG02011_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03486_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA19650_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03065_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_NA20509_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00512_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA12878_giab_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA19238_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG02818_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA18534_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00732_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_HG00731_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA18939_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_NA19983_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00096_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta"]

"""
read the assembly annotations
"""
vntr2anno = defaultdict(lambda: defaultdict(int)) # vntr -> annos -> # of samples supporting

for asm in assembly:
    file = input_path + branch + "/" + asm + "." + aligner + ".vcf"
    with open(file, 'r') as f:
        lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.startswith("#"):
                continue
            fields = re.split('\t', line.strip('\n'))
            chrm = fields[0]
            start = fields[1]
            infos = re.split(";", fields[7])
            end = re.split("=", infos[0])[1]
            anno = re.split("=", infos[3])[1]
            p_anno = ','.join([motif[6:] for motif in re.split(",", anno)])
            vntr = chrm + ":" + start + "-" + end
            vntr2anno[vntr][p_anno] += 1
                 
                 
"""
convert to dataframe (vntr: # of different alleles)
"""
def to_dataframe (_vntr2anno, ed):
    vntr_list = []
    allele_num = []
    vntr2alleles = defaultdict(int)
    for vntr, alleles in _vntr2anno.items():
        vntr_list.append(vntr)
        allele_num.append(len(alleles))
        vntr2alleles[vntr] = len(alleles)  

    allele_num_srt = [a for _, a in sorted(zip(vntr_list, allele_num), key=lambda pair: pair[0])]
    vntr_list_srt = sorted(vntr_list)
    data = {"vntr" : vntr_list_srt, "allele_num": allele_num_srt, "index" : list(range(0 + ed * len(vntr2alleles) , len(vntr2alleles) + ed * len(vntr2alleles))), "edit distance": [ed] * len(vntr2alleles)}
    out_df = pd.DataFrame(data) 
    return out_df

"""
allow `ed` difference in the annotations in order to be counted as same allele
"""
def count_asm_alleles (_vntr2anno, vntr2anno_out, ed):
    for vntr, alleles in _vntr2anno.items():
        n = len(alleles)
        annos = list(alleles.keys())

        sz = [alleles[anno] for anno in annos]
        dsu = UnionFind.DSU(n, sz)
        
        # check pairwise distance 
        for i, x in enumerate(annos):
            for j, y in enumerate(annos):
                if i >= j:
                    continue
                x_list = re.split(",", x)
                y_list = re.split(",", y)
                if edlib.align(x_list, y_list)['editDistance'] <= ed:
                    dsu.union(i, j)
        
        for i, x in enumerate(annos):
            if annos[dsu.find(i)] in vntr2anno_out[vntr]:
                continue
            vntr2anno_out[vntr][annos[dsu.find(i)]] = dsu.sz(i)

# dataframe 
df_vntr2anno_0 = to_dataframe(vntr2anno, 0)
df_vntr2anno_0['allele_num'].mean()
df_vntr2anno_0 


# vntr2anno_1: # of different alleles (allowing 1 difference in annotation)
vntr2anno_1 = defaultdict(lambda: defaultdict(int))
count_asm_alleles (vntr2anno, vntr2anno_1, 1)
df_vntr2anno_1 = to_dataframe(vntr2anno_1, 1)
df_vntr2anno_1['allele_num'].mean()

vntr2anno_2 = defaultdict(lambda: defaultdict(int))
count_asm_alleles (vntr2anno, vntr2anno_2, 2)
df_vntr2anno_2 = to_dataframe(vntr2anno_2, 2)
df_vntr2anno_2['allele_num'].mean()

vntr2anno_3 = defaultdict(lambda: defaultdict(int))
count_asm_alleles (vntr2anno, vntr2anno_3, 3)
df_vntr2anno_3 = to_dataframe(vntr2anno_3, 3)
df_vntr2anno_3['allele_num'].mean()

frames = [df_vntr2anno_0, df_vntr2anno_1, df_vntr2anno_2, df_vntr2anno_3]
result = pd.concat(frames)

"""
hist
"""
plt.figure(figsize=(15, 12))
ax = sns.histplot(result, x="allele_num", hue="edit distance", multiple="dodge", discrete=True, palette="muted", shrink=.8, legend=True)
ax.set(xlim=(0, 30))
ax.tick_params(axis="x", labelsize=15)
ax.tick_params(axis="y", labelsize=15)
plt.xlabel("number of alleles", fontsize=20)
plt.ylabel("frequency", fontsize=20)
plt.setp(ax.get_legend().get_texts(), fontsize='15') 
plt.setp(ax.get_legend().get_title(), fontsize='20') 
plt.savefig(out_hist_png, dpi=300, format="png")
plt.savefig(out_hist_eps, dpi=300, format="eps")
plt.show()
plt.close()


"""
scatter plot
"""
fig, ax = plt.subplots(figsize=(15, 12))
p1 = sns.scatterplot(data=df_vntr2anno_0, x="index", y="allele_num", hue="allele_num", s=3, palette="muted")
ax.set_xlabel("vntr",fontsize=20)
ax.set_ylabel("number of alleles",fontsize=20)
ax.tick_params(labelsize=15)
ax.get_legend().remove()
plt.savefig(out_scatter_png, dpi=300, format="png")
plt.savefig(out_scatter_eps, dpi=300, format="eps")
plt.show() 
plt.close()



































