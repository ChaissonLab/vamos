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
    print("Looking for region " + chrom + " " + str(start) + " " + str(end));
    if curAnno[0][0] == chrom and curAnno[0][1] >= start:
        print("Cur anno is beyond this position " + str(curAnno[0]))
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
                    
        print("  Returning annotation at " + str(curAnnoPos))
        return (curAnnoPos, motifStats, motifLengths, lengthStats, seqLengths)
    return None
        


nSamples = len(args.samples)  
sampleFiles = [gzip.open(fileName) for fileName in args.samples]

regionFile = open(args.ref)

#
# Each sample has a position in the current file, and motif/length stats at that position.
samplePosAnno = [ [ [None, 0, 0], [], [], [], [] ] for a in sampleFiles ]


for regionLine in regionFile:
    vals=regionLine.split()
    rgnChrom = vals[0]
    rgnStart = int(vals[1])
    rgnEnd   = int(vals[2])
    for i in range(0,nSamples):
        samplePosAnno[i] = AdvanceToPos(rgnChrom, rgnStart, rgnEnd, samplePosAnno[i], sampleFiles[i])
        if samplePosAnno[i][0] == rgnChrom and samplePosAnno[i] != rgnStart:
            print("Skipping region " + str(vals))
    

                    
            
                    



    
    
