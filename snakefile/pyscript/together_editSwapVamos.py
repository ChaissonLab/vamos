#!~/.conda/envs/lra/bin/python3
import os
import sys
import re
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt
from collections import defaultdict 
from matplotlib.collections import PathCollection
import numpy as np
import edlib

os.chdir(r"/project/mchaisso_100/cmb-16/jingwenr/trfCall")

if len(sys.argv) == 8:
    input_path = sys.argv[1]
    out1_png = sys.argv[2]
    out2_png = sys.argv[3]
    out1_eps = sys.argv[4]
    out2_eps = sys.argv[5]
    delta = sys.argv[6]
    software = sys.argv[7]
else:
    print("require input_file, out_consensus_png, out_nonconsensus_png, out_consensus_eps, out_nonconsensus_eps delta software")
    exit(1)

consensus = ["consensus-inside", "nonconsensus"]
# gt_sim_mode = ["r1-cov30", "r2-cov30", "r1-g1-l2to4-cov30", "r2-g1-l2to4-cov30"]
read_type = ["GroundTruth", "CCS", "ONT"]
# gt_sim_dict = {"r1-cov30" : "1% error", "r2-cov30" : "2% error", "r1-g1-l2to4-cov30" : "1% error + 1% gap", "r2-g1-l2to4-cov30" : "2% error + 1% gap"}
gt_sim_dict = {"r1-cov30" : "1% error", "r2-cov30" : "2% error"}
gt_sim_mode = ["r1-cov30", "r2-cov30"]


read_type_list = []
delta_list = []
vntr_list = []
edit_list = []
anno_len = []
gt_sim_mode_list = []
consensus_list = []

for c in consensus:
    if software == 'sdcomp' and c == "nonconsensus":
        continue
    for r in read_type:
        if r == "GroundTruth" and c == "consensus-inside":
            continue
        for g in gt_sim_mode:
            input_file = input_path + "/" + g + "/comp/read-" + r + "/" + software + "-delta-" + delta + "." + c + "/gt_reads.swap.comp.anno.summary.bed"
            with open(input_file, 'r') as f:
                lines = f.readlines()
                for idx, line in enumerate(lines):
                    fields = re.split('\t', line.strip('\n'))
                    vntr_list.append(fields[0])
                    edit_list.append(int(fields[1]))
                    anno_len.append(int(fields[2]))
                    read_type_list.append(fields[3])
                    delta_list.append(float(fields[4]))
                    gt_sim_mode_list.append(gt_sim_dict[g])
                    consensus_list.append(c)

data_dict = {"read_type" : read_type_list, 
             "delta" : delta_list,
             "edit_dist" : edit_list,
             "anno_len" : anno_len,
             "vntr" : vntr_list,
             "relative_error" : [edit_list[i] / anno_len[i] for i in range(len(anno_len))],
             "gt_sim_mode" : gt_sim_mode_list,
             "consensus" : consensus_list}

data = pd.DataFrame(data_dict)


plt.figure(figsize=(18, 12))
sub_data2 = data.loc[(data['consensus'] == 'nonconsensus') & (data['read_type'] == 'GroundTruth')]
sub_data1 = data.loc[data['consensus'] == 'consensus-inside']
frames = [sub_data2, sub_data1]
result = pd.concat(frames)
sns.swarmplot(x ='gt_sim_mode', y ='relative_error', hue="read_type", data = result, \
						s=1.5, color="white", edgecolor="gray", dodge=True)
ax = sns.violinplot(x="gt_sim_mode", y="relative_error", hue="read_type", data = result, \
                        palette="muted", scale="width", inner="quartile", linewidth=2)
ax.set_ylim(-0.1, 1.0)
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[0:3], ["Assembly (-raw_anno)", "HiFi (-liftover + -conseq_anno)", "ONT (-liftover + -conseq_anno)"], fontsize=20)
plt.xlabel("Simulation setting of ground truth sequence", fontsize=30)
plt.ylabel("Relative difference", fontsize=30)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
plt.savefig(out1_png, dpi=300, format="png")
plt.savefig(out1_eps, dpi=300, format="eps")
plt.close() 


plt.figure(figsize=(18, 12))
sns.swarmplot(x ='gt_sim_mode', y ='relative_error', hue="read_type", data = data.loc[data['consensus'] == 'nonconsensus'], \
                      s=1.5, color="white", edgecolor="gray", dodge=True)
ax = sns.violinplot(x="gt_sim_mode", y="relative_error", hue="read_type", data = data.loc[data['consensus'] == 'nonconsensus'], \
                        palette="muted", scale="width",inner="quartile", linewidth=2)
ax.set_ylim(-0.1, 1.0)
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[0:3], ["Assembly (-raw_anno)", "HiFi (-raw_anno)", "ONT (-raw_anno)"], fontsize=20)
plt.xlabel("Simulation setting of ground truth sequence", fontsize=30)
plt.ylabel("Relative difference", fontsize=30)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
plt.savefig(out2_png, dpi=300, format="png")
plt.savefig(out2_eps, dpi=300, format="eps")
plt.close() 