import sys
import os
import logging
import random

import numpy as np
from scipy.stats import fisher_exact

import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

windowFinal = sys.argv[1]
regionFinal = sys.argv[2]
segdupBed = sys.argv[3]
outFile = sys.argv[4]

# read in the window file
winDict = {}
with open(windowFinal) as f:
    for line in f:
        chr,startPointCoor,endPointCoor,startPointIndex,endPointIndex,windowTag,\
            meanRateWindow,meanRateTotal,meanPairDistWindow,meanPairDistTotal,stdPairDistTotal = line.strip().split('\t')
        winDict[(chr,startPointCoor,endPointCoor)] = []

# read in the hotspot file
hotspotDict, numHotspot = {}, 0
with open(regionFinal) as f:
    for line in f:
        chr,start,end,startIndex,endIndex,hotspot,SegmentPairDist = line.strip().split('\t')
        hotspotDict[(chr,start,end)] = hotspot
        if hotspot == 'hotspot': numHotspot += 1

# read in segdup file
segdupDict = {}
with open(segdupBed) as f:
    for line in f:
        chr,start,end = line.strip().split('\t')[:3]
        segdupDict[(chr,start,end)] = 'segdup'

# intersect window with hotspot/normal beds
winDict = bedlib.bedIntersectTag(winDict, hotspotDict)
# intersect window with segdup beds
winDict = bedlib.bedIntersectTag(winDict, segdupDict)

# a: hotspot and in segdup
# b: hotspot but not in segdup
# c: normal and in segdup
# d: normal but not in segdup
a,b,c,d = 0,0,0,0
for (chr,start,end),v in winDict.items():
    if v[-2] == 'hotspot':
        if 'segdup' in v[-1]:
            a += 1
        else:
            b += 1
    elif v[-2] == 'normal':
        if 'segdup' in v[-1]:
            c += 1
        else:
            d += 1
    else:
        logging.info(f'Problem window without clear hotspot tag: {chr}:{start}-{end}, {v[-2]}')

# fisher exact test
table = np.array([[a, b],[c, d]])
OR, p = fisher_exact(table, alternative='two-sided')

# output
out = open(outFile, 'w')
out.write(f'p-value={p}\t{a+b}\t{numHotspot}\n\n')
for (chr,start,end),v in winDict.items():
    out.write('\t'.join([chr,start,end,v[-2],v[-1]]) +'\n')
out.close()

