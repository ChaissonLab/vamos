#!~/.conda/envs/lra/bin/python3
import os
import sys
import re
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt

os.chdir(r"/project/mchaisso_100/cmb-16/jingwenr/trfCall")

if len(sys.argv) == 6:
    input_file = sys.argv[1]
    out_plot_png = sys.argv[2]
    out_plot_eps = sys.argv[3]
    out_box_png = sys.argv[4]
    out_box_eps = sys.argv[5]
else:
    print("require input_file out_plot_png out_plot_eps, out_box_png, out_box_eps")
    exit(1)

#original_motif_size, efficient_motif_size, delta
omotif_size = []
emotif_size = []
delta = []

with open(input_file, 'r') as f:
    lines = f.readlines()
    for idx, line in enumerate(lines):
        fields = re.split('\t', line.strip('\n'))
        omotif_size.append(int(fields[0]))
        emotif_size.append(int(fields[1]))
        delta.append(str(100 * float(fields[2])) + "% quantile" )
        
data = {"omotif_size" : omotif_size, "emotif_size" : emotif_size, "ratio" : [emotif_size[i]/omotif_size[i] for i in range(len(omotif_size))], "delta" : delta}
df = pd.DataFrame(data) 
max_original = df["omotif_size"].max()
size = len(omotif_size)

fig, ax = plt.subplots(figsize=(6, 4))

p1 = sns.scatterplot(data=df, x="omotif_size", y="emotif_size", hue="delta", s=15, ax=ax, palette=sns.color_palette('rocket', n_colors=3))
ax.set(ylim=(0, max_original))
ax.set(xlim=(0, max_original))
ax.set(xlabel="original motif set size", ylabel = "efficient motif set size")
p2 = sns.lineplot(x=range(0, max_original), y=range(0, max_original), linestyle='--', ax=ax)
ax.legend(title='delta', loc='upper left')
plt.savefig(out_plot_png, dpi=300, format="png")
plt.savefig(out_plot_eps, dpi=300, format="eps")
plt.close()

plt.figure(figsize=(15, 12))
ax = sns.boxplot(x="delta", y="ratio", hue="delta", palette="Set2", data=df, width=0.5)
ax.set_xlabel("delta",fontsize=20)
ax.set_ylabel("efficient motifs set size / original motifs set size",fontsize=20)
ax.tick_params(labelsize=15)
ax.legend(title='delta', loc='lower left')
plt.setp(ax.get_legend().get_texts(), fontsize='15') 
plt.setp(ax.get_legend().get_title(), fontsize='15') 
plt.savefig(out_box_png, dpi=300, format="png")
plt.savefig(out_box_eps, dpi=300, format="eps")
plt.close()



