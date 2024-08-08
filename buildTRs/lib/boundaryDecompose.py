#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
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
posList = ['inFA', 'outFile']
optList = ['inTRF', 'inRM']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
Raw output from trf/rm may contain overlapping entries. This program merges
these entries and re-decompose the merged sequence by its motif compositions
using the stringDecomposer algorithm.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput fasta,  e.g. /in/fasta')
parser.add_argument(posList[1], type=str, \
    help='string\toutput file,  e.g. /out/File')
# optional arguments
parser.add_argument('-t', '--'+optList[0], type=str, metavar='', default=None,
    help='string\tinput trf calls, e.g. /in/trf.raw.tsv, default None')
parser.add_argument('-r', '--'+optList[1], type=str, metavar='', default=None,
    help='string\tinput repeatMasker calls, e.g. /in/rm.raw.tsv, default None')

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

    seqDict = {s.id : str(s.seq).upper() for s in SeqIO.parse(inFA, 'fasta')}
    logging.info('Finish reading input fasta.')

    if inTRF or inRM:
        pass
    else:
        sys.exit('Please input at least one of "inTRF" or "inRM"!')

    # read input trf raw.tsv file
    allDict = {}
    if inTRF:
        logging.info('Reading input trf calls...')
        with open(inTRF) as f:
            for index,line in enumerate(f):
                contig,start,end,consensus = line.strip().split()[:4]
                start,end = int(start),int(end)
                consensus = line.strip().split()[14]

                if contig not in allDict: allDict[contig] = []
                if [start,end,end-start,consensus] not in allDict[contig]:
                    allDict[contig].append([start,end,end-start,consensus])
        logging.info('Reading input trf calls finish')
    # read input repeatMasker raw.tsv file
    if inRM:
        logging.info('Reading input repeatMasker calls...')
        with open(inRM) as f:
            for index,line in enumerate(f):
                if index < 2 or line.strip() == '': continue
                fields = line.strip().split()
                if '*' in fields: fields = fields[:-1]
                _,_,_,_,contig,start,end,left,_,rep,family,_,_,_,_ = fields
                start,end = int(start),int(end)
                # handle only Simple_repeat (STRs)
                if family == 'Simple_repeat':
                    consensus = rep.replace(')n','').replace('(','')
                else:
                    continue
                # remove calls that has <2 motifs
                if end-start < 2*len(consensus): continue

                if contig not in allDict: allDict[contig] = []
                if [start,end,end-start,consensus] not in allDict[contig]:
                    allDict[contig].append([start,end,end-start,consensus])
        logging.info('Reading input repeatMasker calls finish')

    # sort all input entries by starting position on each contig
    for contig,TRs in allDict.items():
        allDict[contig] = sorted(TRs, key=lambda x: x[0])

    # grouping different entries
    logging.info('Grouping TRs...')
    groupDict = {contig:[] for contig in allDict.keys()}
    for contig,TRs in allDict.items():
        group, checks = [TRs[0]], 0
        for i,(s1,e1,len1,consensus1) in enumerate(TRs,start=1):

            # skip TRs[0] since its already initialized
            if i == 1: continue

            for s2,e2,len2,consensus2 in group:
                op = overlap(s1,e1,s2,e2)
                #if op >= 0.6*len1 or op >= 0.6*len2:
                if op > 0: checks += 1

            #if checks == len(group):
            if checks > 0:
                group.append([s1,e1,len1,consensus1])
                checks = 0
            else:
                groupDict[contig].append(group)
                group, checks = [[s1,e1,len1,consensus1]], 0

            # append the last group
            if i == len(TRs): groupDict[contig].append(group)
    logging.info('Grouping TRs finish')

    # redecomposition and output
    logging.info('Re-decompose TRs...')
    nLoci, nNested = 0, 0
    out = open(outFile, 'w')
    for contig, groups in groupDict.items():
        for group in groups:
            # sort all entries by length (longest to shortest) so to pick the
            # longest entry first among all entries having the shortest period
            groupSorted = sorted(group, key=lambda x: x[2], reverse=True)

            start = min([g[0] for g in groupSorted])
            end = max([g[1] for g in groupSorted])
            seq = seqDict[contig][start-1:end]
            # filter sequences that has "N" in the assembly
            if 'N' in seq: continue

            cons = [ g[3] for g in groupSorted ]
            periods = [ len(g[3]) for g in groupSorted ]
            pick = periods.index(min(periods))
            con = cons[pick]

            # handle special mono/di-nucleotide inside larger TRs
            # pick the longest sequence directly
            if len(con) <= 2 and groupSorted[pick][2] / (end - start) < 0.5:
                #pick = periods.index(min([ p for p in periods if p != 1 ]))
                con = cons[0]
                logging.info(f'nested STR found at {contig}:{start}-{end}')
                nNested += 1
            nLoci += 1
            # decompostion
            motifs = seqDecompose.stringDecomposer(seq, [con])

            coors = ','.join([ str(g[0])+'-'+str(g[1]) for g in groupSorted ])
            cons = ','.join(cons)
            temp = [contig,start,end,coors,con,cons,','.join(motifs)]
            out.write('\t'.join(map(str, temp)) +'\n')
    out.close()
    logging.info('Re-decompose TRs finish')
    logging.info(f'nested/total loci: {nNested}/{nLoci}')


logging.info('End of Program\n')

