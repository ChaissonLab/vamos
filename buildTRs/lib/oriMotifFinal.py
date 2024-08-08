#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from itertools import count
import os
import sys
import argparse
import re
import datetime
import logging
import seqShift
import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inFiles', 'chrs', 'outPre']
optList = ['inTE', 'inCentro', 'lengthCut', 'uniqueCut']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
TR entries of different samples from "oriMotifLiftDecompose.py" needs to be
aggregated. This program aggregates for each TR locus motifs from different
samples, and does final filtering (length, motif counts, overlap with
transposable elements, centromere) to give the final oriMotif set. Note that
filtering by number of unique motifs is done after cyclic shift adjustment.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput redecomposed TRs,  e.g. /in/Files')
parser.add_argument(posList[1], type=str, \
    help='string\tchr to handle (separated by ","),  e.g. chr1,chr2')
parser.add_argument(posList[2], type=str, \
    help='string\toutput file prefix,  e.g. /out/Pre')
# optional arguments
parser.add_argument('-t', '--'+optList[0], type=str, metavar='', default=None,
    help='string\tTE config for masking, /in/te.bed, default None')
parser.add_argument('-c', '--'+optList[1], type=str, metavar='', default=None,
    help='string\tcentromere config for masking, /in/centro.bed, default None')
parser.add_argument('-l', '--'+optList[2], type=int, metavar='', default=10000,
    help='string\tlength filter cut, e.g. 10000, default 10000')
parser.add_argument('-r', '--'+optList[3], type=str, metavar='', default=500,
    help='string\tnumber of unique motifs filter cut, e.g. 500, default 500')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
logging.info('Parsing Input Arguements...')
for key, value in argsDict.items():
    if key in posList: logging.info(f'Required Argument - {key}: {value}')
    if key in optList: logging.info(f'Optional Argument - {key}: {value}')
    vars()[key] = value # assign values of arguments into global variables
logging.info('Parsing Input Arguements Completed\n')

#--------------------------------------------------------
# global variables and user-defined functions
#--------------------------------------------------------

chrs = chrs.split(',')
def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    logging.info('Reading input files...')
    teDict = {}
    if inTE:
        with open(inTE) as f:
            for line in f:
                chr,start,end,tag = line.strip().split('\t')[:4]
                if 'SINE' in tag or 'LINE' in tag or 'Alu' in tag:
                    teDict[(chr,start,end)] = f'{chr}_{start}-{end}'

    centroDict = {}
    if inCentro:
        with open(inCentro) as f:
            for line in f:
                chr,start,end = line.strip().split('\t')[:3]
                centroDict[chr] = [int(start), int(end)]

    with open(inFiles) as f: files = [ l.strip() for l in f ]
    allDict, consensusDict = {}, {}
    for file in files:
        logging.info(f'Reading file: {file}')
        with open(file) as f:
            for line in f:
                chr,start,end,consensus,motifs = line.strip().split('\t')
                if chr not in chrs: continue

                motifs = motifs.split(',')
                # remove motifs containing N
                motifs = [ m for m in motifs if 'N' not in m ]

                # handle first/last motif
                cut = 0.2*len(consensus)
                if len(motifs) == 2:
                    # remove truncated/excessive edge motifs for 2 motif case
                    # if the two motifs are divergent on length
                    if abs(len(motifs[0])-len(motifs[-1])) >= cut:
                        if abs(len(motifs[0])-len(consensus)) >= cut:
                            motifs = motifs[1:]
                        if abs(len(motifs[-1])-len(consensus)) >= cut:
                            motifs = motifs[:-1]
                if len(motifs) > 2:
                    # remove truncated/excessive edge motifs for >2 motif case
                    if abs(len(motifs[0])-len(consensus)) >= cut:
                        motifs = motifs[1:]
                    if abs(len(motifs[-1])-len(consensus)) >= cut:
                        motifs = motifs[:-1]

                # remove single-base motifs in non-homopolymers
                if len(consensus) != 1:
                    motifs = [ m for m in motifs if len(m) != 1 ]
                # skip locus that has no motifs after filtering
                if motifs == []: continue

                # append motifs
                if (chr,start,end) not in allDict:
                    allDict[(chr,start,end)] = {}
                    consensusDict[(chr,start,end)] = consensus
                for m in motifs:
                    if m not in allDict[(chr,start,end)]:
                        allDict[(chr,start,end)][m] = 0
                    allDict[(chr,start,end)][m] += 1
    logging.info('Reading input files finish')

    logging.info(f'Number of total loci: {len(allDict)}')
    tempDict = {}
    for (chr,start,end),locusDict in allDict.items():
        # filter loci that are in centromere
        if chr in centroDict:
            startC, endC = centroDict[chr]
            opCentro = overlap(int(start), int(end), startC, endC, base=1)
            if opCentro > 0:
                logging.info(f'Centromere: {chr}:{start}-{end}')
                continue
        # sort motifs by counts
        tempDict[(chr,start,end)] = { k:v for k,v in sorted(locusDict.items(), \
                key=lambda item: item[1], reverse=True)}
    allDict = tempDict
    # calculate overlaps with transposable elements
    tempDict = { k:[] for k,v in allDict.items() }
    tempDict = bedlib.bedIntersectTag(tempDict, teDict)
    newDict, opDict = {}, {}
    for (chr,start,end),locusDict in allDict.items():
        if tempDict[(chr,start,end)][-1] != 'NULL':
            op1s, op2s = [0], [0]
            for coor in tempDict[(chr,start,end)][-1].split(','):
                s,e = coor.split('_')[1].split('-')
                start,end,s,e = int(start),int(end),int(s),int(e)
                op = overlap(start,end,s,e)
                op1s.append( round(op/(end-start+1),2) )
                op2s.append( round(op/(e-s+1),2) )
            opDict[(chr,str(start),str(end))] = [max(op1s), max(op2s)]
        else:
            opDict[(chr,str(start),str(end))] = [0, 0]
        newDict[(chr,str(start),str(end))] = locusDict
    # sort the loci by starting pos for output
    allDict = {}
    for (chr,start,end),v in tempDict.items():
        # filter loci that exceed the length cut
        if int(end) - int(start) > lengthCut:
            logging.info(f'Long: {chr}:{start}-{end}')
            continue
        allDict[(chr,start,end)] = newDict[(chr,start,end)]

    out = open(outPre+'.all.tsv', 'w')
    out.write('chr\tstart\tend\tmotifs\tmotifCounts\tconsensus\top1\top2\n')
    logging.info('Shift each sequence...')
    for (chr,start,end),locusDict in allDict.items():

        # shift all motifs to consensus and re-calculate counts
        motifs = list(locusDict.keys())
        consensus = consensusDict[(chr,start,end)]
        _,shiftDict = seqShift.shiftFrame(motifs,use=consensus,lexi=False)
        shifted = [ shiftDict[m][0] for m in motifs ]
        counts = [ locusDict[m] for m in motifs ]
        newCountDict = {}
        for i,m in enumerate(shifted):
            if m not in newCountDict: newCountDict[m] = 0
            newCountDict[m] += counts[i]

        # sort shifted motifs by counts
        tempDict = { k:v for k,v in sorted(newCountDict.items(), \
                key=lambda item: item[1], reverse=True)}
        shifted = ','.join( list(tempDict.keys()) )
        counts = ','.join( map(str,list(tempDict.values())) )

        # filter loci that exceed the unique motifs number cut
        if len(tempDict) > uniqueCut:
            logging.info(f'Diverse: {chr}:{start}-{end}')
            continue

        ops = opDict[(chr,start,end)]
        outList = [chr,start,end,shifted,counts,consensus] + ops
        out.write('\t'.join(map(str,outList)) +'\n')
    logging.info('Decompose each sequence finish')
    out.close()

    logging.info(f'Number of remaining loci: {len(allDict)}')

logging.info('End of Program\n')

