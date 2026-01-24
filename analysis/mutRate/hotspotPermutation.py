import sys
import os
import logging
import random

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

hotspotWindows = sys.argv[1]
randomSeed = sys.argv[2]
outFile = sys.argv[3]

# read in the outWindowRaw file
chrs, startPointsCoor, endPointsCoor, startPointsIndex, endPointsIndex, windowTags = [], [], [], [], [], []
meanRateWindows, meanRateTotals, meanPairDistWindows, meanPairDistTotals, stdPairDistTotals = [], [], [], [], []
idxs = []
with open(hotspotWindows) as f:
    for i,line in enumerate(f):
        chr,startPointCoor,endPointCoor,startPointIndex,endPointIndex,\
            windowTag,meanRateWindow,meanRateTotal,meanPairDistWindow,meanPairDistTotal,stdPairDistTotal = line.strip().split('\t')

        chrs.append(chr)
        startPointsCoor.append(startPointCoor)
        endPointsCoor.append(endPointCoor)
        startPointsIndex.append(startPointIndex)
        endPointsIndex.append(endPointIndex)

        windowTags.append(windowTag)
        meanRateWindows.append(meanRateWindow)
        meanRateTotals.append(meanRateTotal)
        meanPairDistWindows.append(meanPairDistWindow)
        meanPairDistTotals.append(meanPairDistTotal)
        stdPairDistTotals.append(stdPairDistTotal)

        idxs.append(i)

random.seed(randomSeed)
random.shuffle(idxs)

out = open(outFile, 'w')
for i,idx in enumerate(idxs):
    oldItems = [chrs[i],startPointsCoor[i],endPointsCoor[i],startPointsIndex[i],endPointsIndex[i]]
    newItems = [windowTags[idx],meanRateWindows[idx],meanRateTotals[idx],meanPairDistWindows[idx],meanPairDistTotals[idx],stdPairDistTotals[idx]]
    out.write('\t'.join(oldItems+newItems) +'\n')
out.close()

