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

hotspotWindows = sys.argv[1]
CHR = sys.argv[2]
outFile = sys.argv[3]


def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)


logging.info(f'handling {CHR}...')

# 1. meanRateWindow and meanRateTotal were used to decide if a window is hotspot
# 2. hotspot windows are to be merged if they are close enough
# 3. meanPairDistWindow and meanPairDistTotal are to be used to decide if a hotspot segment should be filtered by overall small tree branches

# for permutation test in complex pipeline analysis, the step where permutation happens depends on what effect needs to be permuted.
# for the artifact that is to be removed later on, permutation should happen at a step where this artifact is not disrupted.

# So, in the permutation test, we shuffle window tags before merging and final hotspot segment filtering


##### all remaining code re-uses data from the "outWindowRaw" output file from hotspotWindow.py
# read in the outWindowRaw file
startPointsCoor, endPointsCoor, startPointsIndex, endPointsIndex, windowTags = [], [], [], [], []
meanPairDistsWindow = []
with open(hotspotWindows) as f:
    for line in f:
        chr,startPointCoor,endPointCoor,startPointIndex,endPointIndex,windowTag,\
            meanRateWindow,meanRateTotal,meanPairDistWindow,meanPairDistTotal,stdPairDistTotal = line.strip().split('\t')

        startPointsCoor.append(startPointCoor)
        endPointsCoor.append(endPointCoor)
        startPointsIndex.append(startPointIndex)
        endPointsIndex.append(endPointIndex)
        windowTags.append(windowTag)

        meanPairDistsWindow.append(meanPairDistWindow)

# merge windows to produce more smooth result
k, currentTag = 0, windowTags[0]
startPointsIndexMerge, endPointsIndexMerge = [], []
startPointsCoorMerge, endPointsCoorMerge = [], []
numWindows, segmentTags = [], []
# initialize
startPointIndexMerge, endPointIndexMerge = startPointsIndex[0], endPointsIndex[0]
startPointCoorMerge, endPointCoorMerge = startPointsCoor[0], endPointsCoor[0]
for i,windowTag in enumerate(windowTags):
    # same tag and no centromere/bin break, keep merging
    if windowTag == currentTag and int(startPointsCoor[i])-int(endPointCoorMerge) < 1000000:
        k += 1
        endPointIndexMerge = endPointsIndex[i]
        endPointCoorMerge = endPointsCoor[i]
    # change tag or centromere break, record this segment and start new segment if this is not the very last window
    else:
        segmentTags.append(currentTag)
        numWindows.append(k)
        startPointsIndexMerge.append(startPointIndexMerge)
        startPointsCoorMerge.append(startPointCoorMerge)
        endPointsIndexMerge.append(endPointIndexMerge)
        endPointsCoorMerge.append(endPointCoorMerge)

        currentTag = windowTag
        k = 1
        startPointIndexMerge = startPointsIndex[i]
        startPointCoorMerge = startPointsCoor[i]
        endPointIndexMerge = endPointsIndex[i]
        endPointCoorMerge = endPointsCoor[i]
# record the very last segment
segmentTags.append(currentTag)
numWindows.append(k)
startPointsIndexMerge.append(startPointIndexMerge)
startPointsCoorMerge.append(startPointCoorMerge)
endPointsIndexMerge.append(endPointIndexMerge)
endPointsCoorMerge.append(endPointCoorMerge)

# discard normal segment with less than certain number of windows (except for edge cases)
cut = 4
segmentTagsFilter = []
numWindowsFilter = []
startPointsIndexMergeFilter = []
startPointsCoorMergeFilter = []
endPointsIndexMergeFilter = []
endPointsCoorMergeFilter = []
for i,k in enumerate(numWindows):
    if k < cut and segmentTags[i] == 'normal':
        if i != 0 and i != len(numWindows)-1: # not chromosome ends
            leftGap = int(startPointsCoorMerge[i]) - int(endPointsCoorMerge[i-1])
            rightGap = int(startPointsCoorMerge[i+1]) - int(endPointsCoorMerge[i])
            if leftGap < 1000000 and rightGap < 1000000: # not right/left of centromere/bin
                continue
    segmentTagsFilter.append(segmentTags[i])
    numWindowsFilter.append(k)
    startPointsIndexMergeFilter.append(startPointsIndexMerge[i])
    startPointsCoorMergeFilter.append(startPointsCoorMerge[i])
    endPointsIndexMergeFilter.append(endPointsIndexMerge[i])
    endPointsCoorMergeFilter.append(endPointsCoorMerge[i])

# now finally merge adjacent segments with the same tag
currentTag = segmentTagsFilter[0]
startPointsIndexFinal = []
startPointsCoorFinal = []
endPointsIndexFinal = []
endPointsCoorFinal = []
segmentTagsFinal = []

startPointIndexFinal, endPointIndexFinal = startPointsIndexMergeFilter[0], endPointsIndexMergeFilter[0]
startPointCoorFinal, endPointCoorFinal = startPointsCoorMergeFilter[0], endPointsCoorMergeFilter[0]
for i,segmentTagFilter in enumerate(segmentTagsFilter):
    # same tag and no centromere break, keep merging
    if segmentTagFilter == currentTag and int(startPointsCoorMergeFilter[i])-int(endPointCoorFinal) < 1000000:
        endPointIndexFinal = endPointsIndexMergeFilter[i]
        endPointCoorFinal = endPointsCoorMergeFilter[i]
    # change tag or centromere break, record this segment and start new segment if this is not the very last window
    else:
        segmentTagsFinal.append(currentTag)
        startPointsIndexFinal.append(startPointIndexFinal)
        startPointsCoorFinal.append(startPointCoorFinal)
        endPointsIndexFinal.append(endPointIndexFinal)
        endPointsCoorFinal.append(endPointCoorFinal)

        currentTag = segmentTagFilter
        startPointIndexFinal = startPointsIndexMergeFilter[i]
        startPointCoorFinal = startPointsCoorMergeFilter[i]
        endPointIndexFinal = endPointsIndexMergeFilter[i]
        endPointCoorFinal = endPointsCoorMergeFilter[i]
# record the very last segment
segmentTagsFinal.append(currentTag)
startPointsIndexFinal.append(startPointIndexFinal)
startPointsCoorFinal.append(startPointCoorFinal)
endPointsIndexFinal.append(endPointIndexFinal)
endPointsCoorFinal.append(endPointCoorFinal)

outRaw = outFile.replace('.tsv', '.merge.raw.tsv')
out = open(outRaw, 'w')
# output final segments
for i,segmentTagFinal in enumerate(segmentTagsFinal):
    # treat hotspot of only one window as noise
    if segmentTagFinal == 'hotspot' and int(endPointsIndexFinal[i]) - int(startPointsIndexFinal[i]) < 30:
        segmentTagFinal = 'normal'
    temp = [ startPointsCoorFinal[i], endPointsCoorFinal[i], startPointsIndexFinal[i], endPointsIndexFinal[i] ]
    out.write('\t'.join(map(str, [CHR]+temp+[segmentTagFinal])) +'\n')
out.close()


# remove hotspot of constantly small tree branch
regionDict = {}
with open(outRaw) as f:
    for line in f:
        chr,startS,endS,startSI,endSI,hotspot = line.strip().split('\t')
        regionDict[(chr,startS,endS)] = [startSI,endSI,hotspot]
windowDict = {}
with open(hotspotWindows) as f:
    for line in f:
        chr,startPointCoor,endPointCoor,startPointIndex,endPointIndex,windowTag,\
            meanRateWindow,meanRateTotal,meanPairDistWindow,meanPairDistTotal,stdPairDistTotal = line.strip().split('\t')
        windowDict[(chr,startPointCoor,endPointCoor)] = f'{meanPairDistWindow}'
meanPairDistTotal = float(meanPairDistTotal)
stdPairDistTotal = float(stdPairDistTotal)

regionDict = bedlib.bedIntersectTag(regionDict, windowDict)

filterDict = {}
outHotspotFinal = outFile.replace('.tsv', '.merge.final.tsv')
out = open(outHotspotFinal, 'w')
for (chr,start,end),(startSI,endSI,hotspot,tags) in regionDict.items():
    tags = [ float(tag) for tag in tags.split(',') ]
    segmentPairDist = sum(tags) / len(tags)
    if hotspot == 'hotspot':
        if segmentPairDist > meanPairDistTotal - stdPairDistTotal:
        #if segmentPairDist > percentile25PairDistTotal:
            out.write('\t'.join([chr,start,end,startSI,endSI,hotspot,str(segmentPairDist)]) +'\n')
        else:
            filterDict[((chr,start,end))] = f'{chr}:{start}-{end}'
            logging.info(f'filtered hotspot: {chr}:{start}-{end}')
    else:
        out.write('\t'.join([chr,start,end,startSI,endSI,hotspot,str(segmentPairDist)]) +'\n')
out.close()

# finally filter the window data to remove corresponding hotspot regions
winDict = {}
outWinFinal = outFile.replace('.tsv', '.window.filterHotspot.tsv')
out = open(outWinFinal, 'w')
with open(hotspotWindows) as f:
    for line in f:
        chr,startPointCoor,endPointCoor,startPointIndex,endPointIndex,windowTag,\
            meanRateWindow,meanRateTotal,meanPairDistWindow,meanPairDistTotal,stdPairDistTotal = line.strip().split('\t')
        if chr != CHR:
            continue
        winDict[(chr,startPointCoor,endPointCoor)] = [line]
logging.info(f'number of windows before filtering: {len(winDict)}')
if filterDict:
    winDict = bedlib.bedIntersectRemove(winDict, filterDict)
logging.info(f'number of windows before filtering: {len(winDict)}')
for k,v in winDict.items():
    out.write(v[0])
out.close()

# deprecated in permutation test, since the window data has been shuffled
#####
# finally filter the mutation rate and tree data to remove corresponding hotspot regions
#####
'''
trDict = {}
outRateFinal = outFile.replace('.tsv', '.mutRate.filterHotspot.tsv')
out = open(outRateFinal, 'w')
with open('../s04_rate/all.mutRate.tsv') as f:
    for i,line in enumerate(f):
        if i == 0:
            out.write(line)
            continue
        chr,start,end,startTree,endTree,TR2Tree,nSNP,treeLen = line.strip().split('\t')[:8]
        if chr != CHR:
            continue
        trDict[(chr,start,end)] = [line]
treeDict = {}
with open('../s04_rate/all.tree.tsv') as f:
    for line in f:
        chr,startT,endT,nSNP,treeLen,aveEd,aveDist = line.strip().split('\t')
        if chr != CHR: continue
        treeDict[(chr,startT,endT)] = [line]

# filter TRs in the filtered hotspot regions
logging.info(f'number of TRs before filtering: {len(trDict)}')
if filterDict:
    trDict = bedlib.bedIntersectRemove(trDict, filterDict)
    treeDict = bedlib.bedIntersectRemove(treeDict, filterDict)
logging.info(f'number of TRs after filtering: {len(trDict)}')
for k,v in trDict.items():
    out.write(v[0])
out.close()

outTreeFinal = outFile.replace('.tsv', '.tree.filterHotspot.tsv')
out = open(outTreeFinal, 'w')
for k,v in treeDict.items():
    out.write(v[0])
out.close()
'''

