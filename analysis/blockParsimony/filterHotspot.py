import sys
import os
import logging
import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

mutRate = '../../relate/2025-12-15_window2Mb-100kb-expanded/s06_hotspotWindow/all.mutRate.filterHotspot.tsv'
allBlock = 'all.blockTrue.tsv'
selectedBlock = 'all.selectedBlock.reduced.sorted.tsv'

allDict = {}
with open(mutRate) as f:
    for i,line in enumerate(f):
        if i == 0: continue
        fields = line.strip().split()
        chr,start,end = fields[:3]
        allDict[(chr,start,end)] = [-1,-1]

with open(allBlock) as f:
    for line in f:
        fields = line.strip().split(':')[1].split()
        chr,start,end,annotated,block,temp1,temp2,totalAffected = fields[:8]
        patterns = fields[8:]
        if (chr,start,end) in allDict:
            temp = [ len(p.split(',')) for p in patterns if p != 'NA' ]
            temp = sum(temp) / int(annotated)
            allDict[(chr,start,end)] = [0,temp]

with open(selectedBlock) as f:
    for line in f:
        fields = line.strip().split()
        chr,start,end = fields[:3]
        if (chr,start,end) in allDict:
            allDict[(chr,start,end)][0] += 1

out = open('all.mutRate.filterHotspot.tagBlock.tsv', 'w')
#out.write('\t'.join(['chr','start','end','numPattern','aveBlock']) +'\n')
for (chr,start,end),(numPattern,aveBlock) in allDict.items():
    out.write('\t'.join([chr,start,end,str(numPattern),str(aveBlock)]) +'\n')
out.close()

