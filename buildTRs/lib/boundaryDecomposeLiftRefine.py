#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
import seqBasic
import logging

from Bio import SeqIO
import seqDecompose

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inLift', 'outPre']
optList = ['byChr']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
After liftover by "boundaryDecomposeLift.py", sequences from different contigs
may be mapped to have overlaps on the same genomic region. This program merges
or removes redundant entries and re-decompose the merged sequence by its motif
compositions using the stringDecomposer algorithm. This program also converts
reversely aligned sequences to the reverse complement.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput liftover result,  e.g. /in/lift.tsv')
parser.add_argument(posList[1], type=str, \
    help='string\toutput tsv file prefix,  e.g. /out/Pre')
# optional arguments
parser.add_argument('-c', '--'+optList[0], action='store_true', default=False, \
    help='bool\tgive separate result by chromosomes,  default False')

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

def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    # read input trf liftover file
    scoreDict, allDict = {}, {}
    with open(inLift) as f:
        for line in f:
            fields = line.strip().split()
            chr,start,end,reverse,supp,contig,contigStart,contigEnd = fields[:8]
            intervals,consensus,premotifs,motifs = fields[8:]
            if chr not in allDict:
                scoreDict[chr] = {}
                allDict[chr] = []
            if contig not in scoreDict[chr]: scoreDict[chr][contig] = 0
            scoreDict[chr][contig] += 1

            # handle reverse complement
            if reverse == 'True':
                premotifs = premotifs.split(',')
                motifs = motifs.split(',')
                consensus = seqBasic.reComDNA(consensus)
                premotifs = reversed([seqBasic.reComDNA(p) for p in premotifs])
                motifs = reversed([seqBasic.reComDNA(m) for m in motifs])
                premotifs = ','.join(premotifs)
                motifs = ','.join(motifs)
            temp = [ chr,start,end,reverse,supp,contig,contigStart,contigEnd ]
            temp += [ intervals,consensus,premotifs,motifs ]
            temp = '\t'.join(temp) + '\n'

            start,end = int(start),int(end)
            if [start,end,end-start,supp,contig,temp] not in allDict[chr]:
                allDict[chr].append([start,end,end-start,supp,contig,temp])

    # grouping different entries
    logging.info('Grouping TRs...')
    groupDict = {chr:[] for chr in allDict.keys()}
    for chr,TRs in allDict.items():
        TRsort = sorted(TRs, key=lambda x: int(x[0]))
        group, checks = [TRsort[0]], 0
        for i,(s1,e1,len1,supp1,contig1,line1) in enumerate(TRsort,start=1):

            # skip TRsort[0] since its already initialized
            if i == 1: continue

            for s2,e2,len2,supp2,contig2,line2 in group:
                op = overlap(s1,e1,s2,e2)
                #if op >= 0.6*len1 or op >= 0.6*len2: checks += 1
                if op > 0: checks += 1

            #if checks == len(group):
            if checks > 0:
                group.append([s1,e1,len1,supp1,contig1,line1])
                checks = 0
            else:
                groupDict[chr].append(group)
                group, checks = [[s1,e1,len1,supp1,contig1,line1]], 0

            # append the last group
            if i == len(TRsort): groupDict[chr].append(group)
    logging.info('Grouping TRs finish')

    # handle each group and output
    logging.info('Removing redundant entries...')
    if byChr:
        for chr,groups in groupDict.items():
            out = open(f'{outPre}.{chr}.tsv', 'w')
            for group in groups:
                if len(group) == 1:
                    out.write(group[0][-1])
                else:
                    scores = [ scoreDict[chr][g[4]] for g in group ]
                    pick = scores.index(max(scores))
                    out.write(group[pick][-1])
            out.close()
    else:
        out = open(f'{outPre}.tsv', 'w')
        for chr,groups in groupDict.items():
            for group in groups:
                if len(group) == 1:
                    out.write(group[0][-1])
                else:
                    scores = [ scoreDict[chr][g[4]] for g in group ]
                    pick = scores.index(max(scores))
                    out.write(group[pick][-1])
        out.close()
    logging.info('Removing redundant entries finish')


logging.info('End of Program\n')

