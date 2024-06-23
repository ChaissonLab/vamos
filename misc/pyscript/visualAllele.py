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
import multiprocessing 

os.chdir(r"/project/mchaisso_100/cmb-16/jingwenr/trfCall")

if len(sys.argv) == 3:
    aligner = sys.argv[1]
    branch = sys.argv[2]
else:
    print("require aligner branch-delta")
    exit(1)
    

input_path = "/project/mchaisso_100/cmb-16/jingwenr/trfCall/2022-03-04/hgsvc_2022-03-04/result-64/assembly/"
assembly = ["v12_HG00512_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA20509_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03065_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG02011_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA19650_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG03486_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_HG00732_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA18534_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG00096_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_NA19983_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA18939_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG00731_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA19238_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG02818_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG01596_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG02587_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00513_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_HG01114_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03371_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG00864_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03683_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03009_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_NA20847_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA12329_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03732_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG01505_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00171_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03125_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA19239_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_HG02492_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03009_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03683_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG00864_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA20847_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG01596_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03371_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG01114_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00513_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG02587_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03125_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG00171_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG02492_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_NA19239_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA12329_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG01505_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG03732_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta","v12_HG02011_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03486_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA19650_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG03065_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_NA20509_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00512_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA12878_giab_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA19238_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_HG02818_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_NA18534_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00732_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta","v12_HG00731_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta","v12_NA18939_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_NA19983_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta","v12_HG00096_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta"]

class Allele:
    def __init__ (self, aligner, branch):
        self.aligner = aligner
        self.branch = branch
        self.vntr2anno = defaultdict(lambda: defaultdict(list)) # vntr -> assembly -> annos
        self.vntr_space = {}
        self.read()
        self.appendLength()
        self.path = input_path + self.branch + "/plot" + self.aligner
        print(self.path)
        if not os.path.exists(self.path):
            os.makedirs(self.path)
            print("The new directory is created!")

    def read (self):
        """
        read the assembly annotations
        """
        for asm in assembly:
            file = input_path + self.branch + "/" + asm + "." + self.aligner + ".vcf"
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
                    p_anno = [int(motif[6:]) for motif in re.split(",", anno)]
                    vntr = chrm + ":" + start + "-" + end
                    self.vntr2anno[vntr][asm] = p_anno

    def appendLength(self):
        """
        make annotation the same length per VNTR
        """
        for vntr in self.vntr2anno.keys():
            max_length = 0
            for asm in self.vntr2anno[vntr].keys():
                max_length = max(max_length, len(self.vntr2anno[vntr][asm]))

            for asm in self.vntr2anno[vntr].keys():
                if len(self.vntr2anno[vntr][asm]) < max_length:
                    self.vntr2anno[vntr][asm].extend([257] * (max_length -  len(self.vntr2anno[vntr][asm])))   
                    self.vntr_space[vntr] = 1
                else:
                    self.vntr_space[vntr] = 0


    """
    heatmap for one VNTR
    """
    def heatmap(self, vntr):
        asms = list(self.vntr2anno[vntr].keys())
        data_anno = [self.vntr2anno[vntr][asm] for asm in asms]
        np_annp = np.array(data_anno)
        max_length = max([max(self.vntr2anno[vntr][asm]) for asm in asms]) + 1
        if self.vntr_space[vntr] == 1:
            cmap = sns.color_palette("deep", 256)
            cmap.append((1,1,1)) 
        else:
            cmap = sns.color_palette("deep", max_length)
        print(max_length)
        
        # the index of the position of yticks
#         num_ticks = len(asms)
#         yticks = np.linspace(0, len(asms) - 1, num_ticks, dtype=np.int)
        # the content of labels of these yticks
#         yticklabels = [asms[idx] for idx in yticks]

        plt.figure(figsize=(8, 6))
        ax = sns.heatmap(np_annp, cmap=cmap)
#         ax = sns.heatmap(np_annp, cmap=cmap, yticklabels=yticklabels)
#         ax.set_yticklabels(yticklabels, rotation=0, fontsize="10")
        plt.title(vntr)
        
        plt.show()
        plt.savefig(self.path + "/allele." + vntr + "." + self.aligner + ".png", dpi=100, format="png")
        plt.close()


allele = Allele(aligner, branch)
vntrs = list(allele.vntr2anno.keys())
vntrs.sort()
print(len(vntrs))

def mp_worker(idx_vntr):
    try:
        allele.heatmap(vntrs[idx_vntr])
        return idx_vntr
    except Exception as e:
        print('Caught exception in worker thread')  
        traceback.print_exc() 
        print()
        raise e

def mp_handler():
    pool = multiprocessing.Pool(32)
    all_vntr = range(len(vntrs))

    try:
        for idx_vntr in pool.imap_unordered(mp_worker, all_vntr):
            print("processing VNTR " + str(idx_vntr))
    except:
        print ("Outer exception caught!")
    pool.close()
    pool.join()
    print("finish allele")

if __name__=='__main__':
    mp_handler()



# if __name__=="__main__":
#   allele = Allele(aligner, branch)
#   vntrs = list(allele.vntr2anno.keys())
#   vntrs.sort()
#   print(len(vntrs))
#   for vntr in vntrs:
#       print(vntr)
#       allele.heatmap(vntr)
#   print("finish allele")



  # mm2_allele = Allele("mm2", "anno-delta-0.2")
  # vntrs = list(mm2_allele.vntr2anno.keys())
  # vntrs.sort()
  # print(len(vntrs))
  # for vntr in vntrs:
  #     print(vntr)
  #     mm2_allele.heatmap(vntr)

  # print("finish mm2_allele")

  # lra_allele = Allele("lra", "anno-delta-0.5")
  # vntrs = list(lra_allele.vntr2anno.keys())
  # vntrs.sort()
  # print(len(vntrs))
  # for vntr in vntrs:
  #     print(vntr)
  #     lra_allele.heatmap(vntr)

  # print("finish lra_allele")

  # original_mm2_allele = Allele("mm2", "anno-original-motifs")
  # vntrs = list(original_mm2_allele.vntr2anno.keys())
  # vntrs.sort()
  # print(len(vntrs))
  # for vntr in vntrs:
  #     print(vntr)
  #     original_mm2_allele.heatmap(vntr)

  # print("finish original_mm2_allele")

  # original_lra_allele = Allele("lra", "anno-original-motifs")
  # vntrs = list(original_lra_allele.vntr2anno.keys())
  # vntrs.sort()
  # print(len(vntrs))
  # for vntr in vntrs:
  #     print(vntr)
  #     original_lra_allele.heatmap(vntr) 

  # print("finish original_lra_allele")

