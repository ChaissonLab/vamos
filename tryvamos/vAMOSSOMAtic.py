#!/usr/bin/env python

import numpy as np

def ParseAnno(anno, refMotifs):
    vals=anno.split(":")
    if len(vals) < 3:
        return None
    hap=int(vals[0])
    trLen=int(vals[1])
    if len(vals[2]) > 0:
        annoMotifs=[int(i) for i in vals[2].split(",")]
    else:
        annoMotifs = []
    trSeqLen = sum(len(refMotifs[i]) for i in annoMotifs)
    return (hap,trLen, trSeqLen)

def ParseLine(line):
    vals=line.split()
    if len(vals) < 4:
        return None
    motifs=vals[3].split(",")
    annos =vals[4].split(";")
    readAnnos = [ParseAnno(a, motifs) for a in annos]
    return readAnnos

def GetVal(annoPairs, p):
    if len(annoPairs) > 0:
        v = [a[p] for a in annoPairs]
        return v
    else:
        return None

def GetLenByMotif(annoPairs):
    return GetVal(annoPairs, 0)

def GetLenBySeq(annoPairs):
    return GetVal(annoPairs, 1)


def GetPairStats(annoPairs, p):
    if len(annoPairs) > 0:
        v = [a[p] for a in annoPairs]
        return (np.mean(v), np.std(v))
    else:
        return 0

def GetMotifStats(annoPairs):
    return GetPairStats(annoPairs,0)

def GetLengthStats(annoPairs):
    return GetPairStats(annoPairs,1)


import sys
import gzip
import argparse 
parser = argparse.ArgumentParser(description="Program to look for somatic variants in TR sequences using population data.")

  # Add positional arguments (required)
  # Add optional arguments with flags and types
parser.add_argument("-s", "--samples", nargs="+", default=[],
                      help="Sample annotation files", required=True)
parser.add_argument("-r", "--ref", required=True,
                    help="Reference annotations, output will be in order of the ref.")

  # Parse arguments 
args = parser.parse_args()

# Access arguments using the parsed object 'args'n

# Your program logic here using the parsed arguments


def AdvanceToPos(chrom, start, end, curAnno, annoFile):
    lineNumber=0
    #
    # Check to see if the state of curAnnoPos is beyond the current start-end.
    # If this is the case, no need to parse a new region.
    #
    if curAnno[0][0] == chrom and curAnno[0][1] >= start:
        return curAnno
    
    while curAnno[0][0] == None or curAnno[0][0] != chrom or curAnno[0][1] < start:


        line = annoFile.readline().decode("utf-8")
        if len(line) > 0 and line[0] == "#":
            continue

        vals=line.split()
        annoChrom = vals[0]
        annoStart = int(vals[1])
        annoEnd   = int(vals[2])
        curAnnoPos = [annoChrom, annoStart, annoEnd]
#
        if annoStart != start:
            print("  Potential out of sync with " + str(vals[0:3]) + "\t" + str(start) + "\t" + str(end))
        if annoChrom is not None and annoChrom < chrom:
            continue
        if annoChrom == chrom and annoStart < start:
            print("  Skipping annotation " + str(vals[0:3]))
            continue
            
        lineNumber+=1
    
        hapAnnos = [[] for i in range(0,3)]
        summaries=ParseLine(line)
        motifStats = [[] for i in range(0,3)]
        lengthStats = [[] for i in range(0,3)]
        motifLengths = [[] for i in range(0,3)]
        seqLengths = [[] for i in range(0,3)]
        if summaries is not None:
            for summary in summaries:
                if summary is None:
                    continue
                if summary[0] < 0 or summary[0] > 2:
                    print("ERROR parsing haplotype at line: " + str(lineNumber) + "\t" + str(summary[0]) + "\n" + line)
                    sys.exit(1)
                hapAnnos[summary[0]].append(summary[1:3])
            for i in range(0,3):
                if len(hapAnnos[i]) > 0:
                    motifStats[i] = GetMotifStats(hapAnnos[i])
                    lengthStats[i] = GetLengthStats(hapAnnos[i])
                    motifLengths[i] = GetLenByMotif(hapAnnos[i])
                    seqLengths[i] = GetLenBySeq(hapAnnos[i])                    
        return (curAnnoPos, motifStats, motifLengths, lengthStats, seqLengths)
    return None
        


nSamples = len(args.samples)  
sampleFiles = [gzip.open(fileName) for fileName in args.samples]

regionFile = open(args.ref)

#
# Each sample has a position in the current file, and motif/length stats at that position.
samplePosAnno = [ [ [None, 0, 0], [], [], [], [] ] for a in sampleFiles ]

#
# Order for sample pos anno: 
# spa[0]  - coordinates of variant
# spa[1]  - motif stats
# spa[1][0,1,2] - unphased, hap1, hap2 motif mean/sd
# spa[2][0,1,2] - unphased, hap1, hap2 motif lengths
# spa[3][0,1,2] - unphased, hap1, hap2 length mean/sd
# spa[4][0,1,2] - unphasdd, hap1, hap2 seq lengths

def IsPhased(s):
    if len(s[1][1]) > 0 and len(s[1][2]) > 0:
        return True
    else:
        return False

def ExtractValue(spa, kept, phased, usePhase, val, meanOrStd, minReads):
    if kept is False:
        return []

    idx = [0]
    if usePhase and phased:
        idx = [1,2]

    readThres = [len(spa[val+1][p]) for p in idx]
    res = []
    for i in range(0,len(readThres)):
        if readThres[i] > minReads:
            res += [spa[val][idx[i]][meanOrStd]]
    return res

def ExtractValues(spaRow, kept, phased, usePhase, val, meanOrStd, minReads):
    allRes = []
    for i in range(0,len(spaRow)):
        allRes += ExtractValue(spaRow[i], kept[i], phased[i], usePhase, val, meanOrStd, minReads)
    return allRes



getMotif = 1
getLen   = 3

getMean = 0
getStd  = 1

for regionLine in regionFile:
    vals=regionLine.split()
    rgnChrom = vals[0]
    rgnStart = int(vals[1])
    rgnEnd   = int(vals[2])
    for i in range(0,nSamples):
        samplePosAnno[i] = AdvanceToPos(rgnChrom, rgnStart, rgnEnd, samplePosAnno[i], sampleFiles[i])
    #
    # Now samples have been loaded for the next
    #
    keep = [False ] * nSamples
    phased = [False ] * nSamples
    nKept=0
    nPhased =0
    for i in range(0,nSamples):
        if samplePosAnno[i][0][0] == rgnChrom and samplePosAnno[i][0][1] == rgnStart:
            keep[i]= True
            nKept+=1
            isPhased = IsPhased(samplePosAnno[i])
            phased[i] = isPhased
            if isPhased:
                nPhased+=1
    if nKept < 0.5 * nSamples:
        continue
    isPhased = False
    if nPhased / nKept > 0.5:
        isPhased = True
    meanVals = ExtractValues(samplePosAnno, keep, phased, isPhased, getMotif, getMean, 5)
    stdVals = ExtractValues(samplePosAnno, keep, phased, isPhased, getMotif, getStd, 5)

    print(str(vals[0:3]))
    print(meanVals)
    print(stdVals)
    print("\n")
    
    
    

                    
            
                    



    
    
