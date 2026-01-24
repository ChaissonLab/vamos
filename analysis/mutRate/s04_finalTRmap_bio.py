import sys
import os
import logging
import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

mutRate = '../s06_hotspotWindow/all.mutRateBio.filterHotspot.tsv'
hotspot = '../s06_hotspotWindow/all.sorted.hotspot.tsv'
segdupIdentity = 'grch38_superdups.max_identity.bed'
segdupCNV = 'all.segdup.sort.tsv'
outFile = 'all.mutRateBio.filterHotspot.tagSegdup.tsv'

rateDict = {}
with open(mutRate) as f:
    for i,line in enumerate(f):
        if i == 0:
            header = line.strip().split()
            continue
        chr,start,end,startTree,endTree,TR2Tree,nSNP,treeLen = line.strip().split('\t')[:8]
        rateDict[(chr,start,end)] = [line.strip().split()]

hotDict = {}
with open(hotspot) as f:
    for line in f:
        chr,start,end = line.strip().split('\t')[:3]
        hotDict[(chr,start,end)] = 'hotspot'

segdupDict = {}
with open(segdupIdentity) as f:
    for line in f:
        chr,start,end,identity = line.strip().split()
        segdupDict[(chr,start,end)] = f'{identity}__0__0__0'

with open(segdupCNV) as f:
    for line in f:
        chr,start,end,mapped,annotated,mean,median = line.strip().split()
        identity = segdupDict[(chr,start,end)].split('__')[0]
        segdupDict[(chr,start,end)] = f'{identity}__{mapped}__{mean}__{median}'

rateDict = bedlib.bedIntersectTag(rateDict,hotDict)
rateDict = bedlib.bedIntersectTag(rateDict,segdupDict)

out = open(outFile, 'w')
out.write('\t'.join(header+['hotspot','segDupIdentity','segDupCNV']) +'\n')
for (chr,start,end),(line,hotspot,segdups) in rateDict.items():
    if hotspot == 'NULL':
        hotspot = 'normal'
    else:
        hotspot = 'hotspot'
    if segdups == 'NULL':
        identity = -1
        cnv = -1
    else:
        identities = [ float(s.split('__')[0]) for s in segdups.split(',') ]
        mappeds = [ float(s.split('__')[1]) for s in segdups.split(',') ]
        means = [ float(s.split('__')[2]) for s in segdups.split(',') ]
        medians = [ float(s.split('__')[3]) for s in segdups.split(',') ]
        identity = max(identities)
        cnv = max(means)
    out.write('\t'.join(line+[hotspot,str(identity),str(cnv)]) +'\n')
out.close()

