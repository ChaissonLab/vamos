import sys
import os
import logging
import math
import numpy as np
import bedlib
#import ruptures as rpt
#import matplotlib.pyplot as plt

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

rateFile = sys.argv[1]
CHR = sys.argv[2]
outFile = sys.argv[3]


def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)


logging.info(f'handling {CHR}...')

# get total number of loci on this chromosome
total, starts, ends, rates, pairDists = 0, [], [], [], []
with open('../s04_rate/all.mutRate.tsv') as f:
    for line in f:
        chr,start,end,treeS,treeE,TR2Tree,nSNP,treeLen,winS,winE,period,type,\
                annotated,totalAllele,spanned,totalEd,totalDist,totalLen,rate,rateNormBefore,rateNormAfter = line.strip().split('\t')
        if chr == 'chr': continue # skip header
        rate = float(line.strip().split('\t')[-1])
        # skip constant loci
        if rate == 0: continue
        # skip loci with no tree
        if int(TR2Tree) > 10000: continue
        if chr == CHR:
            total += 1
            starts.append(start)
            ends.append(end)
            rates.append(rate)
            temp = ( int(annotated)*(int(annotated)-1) ) / 2
            pairDists.append(float(totalDist)/temp)
meanRateTotal = np.mean(np.array(rates))
stdRateTotal = np.std(np.array(rates))
medianRateTotal = np.median(np.array(rates))
percentile75RateTotal = np.percentile(np.array(rates), 75)
percentile80RateTotal = np.percentile(np.array(rates), 80)
percentile90RateTotal = np.percentile(np.array(rates), 90)
meanPairDistTotal = np.mean(np.array(pairDists))
stdPairDistTotal = np.std(np.array(pairDists))
percentile25PairDistTotal = np.percentile(np.array(pairDists), 25)

startPointsIndex, startPointsCoor = [1], [starts[0]]
endPointsIndex, endPointsCoor = [], []

windowRates, windowPairDists = [], []
i, windowRate, windowPairDist = 0, [], []
with open('../s04_rate/all.mutRate.tsv') as f:
    for line in f:
        chr,start,end,treeS,treeE,TR2Tree,nSNP,treeLen,winS,winE,period,type,\
                annotated,totalAllele,spanned,totalEd,totalDist,totalLen,rate,rateNormBefore,rateNormAfter = line.strip().split('\t')
        if chr == 'chr': continue # skip header
        rate = float(line.strip().split('\t')[-1])

        # skip constant loci
        if rate == 0: continue
        # skip loci with no tree
        if int(TR2Tree) > 10000: continue
        if chr != CHR: continue
        i += 1
        windowRate.append(rate)
        temp = ( int(annotated)*(int(annotated)-1) ) / 2
        pairDist = float(totalDist)/temp
        windowPairDist.append(pairDist)

        # break point at 20 TRs (one window has 20 TRs)
        if i % 20 == 0:
            # this TR as the end point
            endPointsIndex.append(i)
            endPointsCoor.append(end)
            # next TR as the next start point if this TR is not the last
            if i < total:
                startPointsIndex.append(i+1)
                startPointsCoor.append(starts[i])
            # record rates for this window
            windowRates.append(windowRate)
            windowPairDists.append(windowPairDist)
            windowRate, windowPairDist = [], []
        # a centromere break point
        elif i != 1 and int(start) - int(ends[i-2]) > 1000000:
            # this TR as the end point
            endPointsIndex.append(i-1)
            endPointsCoor.append(ends[i-2])
            # next TR as the next start point if this TR is not the last
            if i < total:
                startPointsIndex.append(i)
                startPointsCoor.append(starts[i-1])
            # record rates for this window
            windowRates.append(windowRate[:-1])
            windowPairDists.append(windowPairDist[:-1])
            windowRate, windowPairDist = [windowRate[-1]], [windowPairDist[-1]]

# add last TR on the chromosome as the last endpoint if it is not
if total not in endPointsIndex:
    endPointsIndex.append(total)
    endPointsCoor.append(ends[total-1])
    windowRates.append(windowRate)
    windowPairDists.append(windowPairDist)
logging.info(f'total number of TRs: {total}')
logging.info(f'total number of windows: {len(startPointsIndex)}, {len(endPointsIndex)}, {len(windowRates)}, {len(windowPairDists)}')


# make crude peak decision and output all windows as a single bed file
outWindowRaw = outFile.replace('.tsv','.window.tsv')
out = open(outWindowRaw, 'w')
windowTags = []
for i,windowRate in enumerate(windowRates):
    windowPairDist = windowPairDists[i]
    if not windowRate:
        # empty window, may happen at a centromere truncation
        print(i)
        meanRateWindow = 0
        stdRateWindow = 0
        medianRateWindow = 0
        percentile75RateWindow = 0
        percentile90RateWindow = 0
        meanPairDistWindow = 0
        continue
    else:
        meanRateWindow = np.mean(np.array(windowRate))
        stdRateWindow = np.std(np.array(windowRate))
        medianRateWindow = np.median(np.array(windowRate))
        percentile75RateWindow = np.percentile(np.array(windowRate), 75)
        percentile90RateWindow = np.percentile(np.array(windowRate), 90)
        meanPairDistWindow = np.mean(np.array(windowPairDist))

    if meanRateWindow > meanRateTotal + stdRateTotal:
        windowTag = 'hotspot'
        windowTags.append('hotspot')
    else:
        windowTag = 'normal'
        windowTags.append('normal')

    temp = [CHR, startPointsCoor[i], endPointsCoor[i], startPointsIndex[i], endPointsIndex[i]]
    out.write('\t'.join(map(str,temp+[windowTag, meanRateWindow, meanRateTotal, meanPairDistWindow, meanPairDistTotal, stdPairDistTotal])) +'\n')

out.close()


# 1. meanRateWindow and meanRateTotal were used to decide if a window is hotspot
# 2. hotspot windows are to be merged if they are close enough
# 3. meanPairDistWindow and meanPairDistTotal are to be used to decide if a hotspot segment should be filtered by overall small tree branches

# for permutation test in complex pipeline analysis, the step where permutation happens depends on what effect needs to be permuted.
# for the artifact that is to be removed later on, permutation should happen at a step where this artifact is not disrupted.

# So, in the permutation test, we shuffle window tags before merging and final hotspot segment filtering

