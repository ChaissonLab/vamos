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

if len(sys.argv) == 7:
    input_path = sys.argv[1]
    out_plot = sys.argv[2]
    consensus = sys.argv[3]
    delta = sys.argv[4]
    gt_sim_mode = sys.argv[5]
    software = sys.argv[6]
else:
    print("require input_file out_plot consensus-mode delta gt_sim_mode software")
    exit(1)

read_type = ["CCS", "ONT", "GroundTruth"]

read_type_list = []
delta_list = []
vntr_list = []
edit_list = []
anno_len = []

for r in read_type:
    input_file = input_path + "/read-" + r + "/" + software + "-delta-" + delta + "." + consensus + "/gt_reads.swap.comp.anno.summary.bed"
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for idx, line in enumerate(lines):
            fields = re.split('\t', line.strip('\n'))
            vntr_list.append(fields[0])
            edit_list.append(int(fields[1]))
            anno_len.append(int(fields[2]))
            read_type_list.append(fields[3])
            delta_list.append(float(fields[4]))
                    
data_dict = {"read_type" : read_type_list, 
             "delta" : delta_list,
             "edit_dist" : edit_list,
             "anno_len" : anno_len,
             "vntr" : vntr_list,
             "relative_error" : [edit_list[i] / anno_len[i] for i in range(len(anno_len))] }

data = pd.DataFrame(data_dict)
# subdata = data.loc[data['delta'] == 0.1]

# ax = sns.violinplot(x="read_type", y="edit_dist", data=data)
# sns.stripplot(x ='read_type', y ='edit_dist', data = data, color= "white", s = 1.5)
# ax.set(xlabel="read type", ylabel = "edit distance between swap and vamos annotation")
# # sns.swarmplot(x ='read_type', y ='edit_dist', data = data, color= "white", s = 4)
# plt.savefig(out_plot, dpi=300, format="png")
# plt.show() 

ax = sns.violinplot(x="read_type", y="relative_error", data=data)
sns.stripplot(x ='read_type', y ='relative_error', data = data, color= "white", s = 1.5)
ax.set(xlabel="read type", ylabel = "relative error between swap and vamos annotation", title=gt_sim_mode)
ax.set(ylim=(-0.2, 1.0))
# sns.swarmplot(x ='read_type', y ='edit_dist', data = data, color= "white", s = 4)
plt.savefig(out_plot, dpi=300, format="png")
plt.show() 