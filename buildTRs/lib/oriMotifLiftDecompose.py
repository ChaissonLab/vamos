#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from itertools import count
import os
import sys
import argparse
import re
import datetime
import logging
import seqDecompose
from Bio import SeqIO

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inBed', 'inFA', 'outFile']
optList = []
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
Lifted TR sequences from "oriMotifLift.py" are not decomposed to the motif
composition. This program applies the stringDecomposer algorithm and decomposes
the TR sequences, each using a single consensus obtained from
"boundaryDecomposeLiftRefineAggregate.py".
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput TR bed (col1-4:chr,start,end,consensus),  e.g. /in/bed')
parser.add_argument(posList[1], type=str, \
    help='string\tinput fasta for decomposition,  e.g. /in/fa')
parser.add_argument(posList[2], type=str, \
    help='string\toutput file,  e.g. /out/file')
# optional arguments


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

chrList = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9', \
           'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
           'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    logging.info('Reading input files...')
    temp = {s.id:str(s.seq).upper() for s in SeqIO.parse(inFA, 'fasta')}
    seqDict = {chr:{} for chr in chrList}
    for id,seq in temp.items():
        chr,coor = id.split(';')[1].split('_')
        start,end = coor.split('-')
        seqDict[chr][(start,end)] = seq

    bedDict = {chr:{} for chr in chrList}
    with open(inBed) as f:
        for line in f:
            chr,start,end,_,consensus = line.strip().split('\t')[:5]
            # skip alternative chromosomes
            if chr not in chrList: continue
            bedDict[chr][(start,end)] = [consensus]
    logging.info('Reading input files finish')

    # sort all seqs by starting coordinate
    logging.info('Sorting all seqs by starting coordinate...')
    temp = {}
    for chr in chrList:
        temp[chr] = {k:v for k,v in sorted(seqDict[chr].items(), \
                key=lambda x: int(x[0][0]))}
    seqDict = temp
    logging.info('Sorting all seqs by starting coordinate finish')

    out = open(outFile, 'w')
    logging.info('Decompose each sequence...')
    for chr,seqs in seqDict.items():
        for (start,end),seq in seqs.items():
            consensus = bedDict[chr][(start,end)]
            motifs = seqDecompose.stringDecomposer(seq, consensus)
            motifs = ','.join(motifs)
            out.write('\t'.join([chr,start,end,consensus[0],motifs]) +'\n')
    logging.info('Decompose each sequence finish')
    out.close()

logging.info('End of Program\n')

