#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from itertools import count
import os
import sys
import argparse
import re
import datetime
import logging
import seqLift
import pysam
from Bio import SeqIO

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inBed', 'inFA', 'inBam', 'outFile']
optList = ['startAdjust', 'endAdjust', 'skipEmpty']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
TR regions obtained from "boundaryDecompose.py" are only on input sequences
(e.g., assembly). This program returns the corresponding coordinates of these
regions based on the alignment bam (liftover from query to reference).
Reference is dependent on the reference used by the input bam files.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput query bed regions,  e.g. /in/bed')
parser.add_argument(posList[1], type=str, \
    help='string\tinput assembly fasta,  e.g. /in/fa')
parser.add_argument(posList[2], type=str, \
    help='string\tinput assembly bam,  e.g. /in/bam')
parser.add_argument(posList[3], type=str, \
    help='string\toutput file,  e.g. /out/file')
# optional arguments
parser.add_argument('-s', '--'+optList[0], type=int, metavar='', \
    default=0, help='int\tboundary adjust for starting coordinate,  default 0')
parser.add_argument('-e', '--'+optList[1], type=int, metavar='', \
    default=0, help='int\tboundary adjust for ending coordinate,  default 0')
parser.add_argument('-d', '--'+optList[2], action='store_true', default=True, \
    help='bool\tskip empty liftover, default True')

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
    seqDict = {s.id:str(s.seq).upper() for s in SeqIO.parse(inFA, 'fasta')}
    samfile = pysam.AlignmentFile(inBam, 'rb')
    bedDict = {}
    with open(inBed) as f:
        for line in f:
            contig,start,end = line.strip().split('\t')[:3]
            if contig not in bedDict: bedDict[contig] = []
            bedDict[contig].append([int(start),int(end),line])
    logging.info('Reading input files finish')

    # store the cigar of each read for downstream analysis
    logging.info('Getting cigar of each read...')
    readsDict = {}
    for read in samfile.fetch():

        # skip alternative chrs, secondary alignments, alignments with <60 mapQ
        if read.reference_name not in chrList: continue
        if read.is_secondary: continue
        if read.mapping_quality < 60: continue

        id = read.query_name
        chr = read.reference_name
        cigar = read.cigarstring
        startAln = read.reference_start
        reverse = read.is_reverse
        supp = read.is_supplementary
        length = int(read.query_alignment_length)

        # parse cigar
        if id not in readsDict: readsDict[id] = []
        readsDict[id].append([chr, startAln, cigar, reverse, supp, length])
    # close bam file to release memory
    samfile.close()
    logging.info('Getting cigar of each read finish')

    logging.info('Lifting each read...')
    # liftover and output
    out = open(outFile, 'w')
    counter = 0
    for contig,aligns in readsDict.items():
        counter += 1
        if counter % 10 == 0: logging.info(f'{counter} reads processed')
        # sort the alignments of each read id so that the primary alignment
        # comes first and all supplementary alignments are sorted by aligned
        # length (longest to shortest)
        temp = [ a for a in aligns if a[4] ] # all supps
        alignsSorted = [ a for a in aligns if not a[4] ]
        alignsSorted += sorted(temp, key=lambda x: x[5], reverse=True)

        maps = []
        # parse all alignments
        for chr,startAln,cigar,reverse,supp,length in alignsSorted:
            map = seqLift.liftOnePass(cigar, startAln, refAsKey=False) 
            maps.append([chr,map,reverse,supp])

        # liftover each query entry
        # check the primary alignment first
        # otherwise check the supplementary alignment from longest to shortest
        if contig not in bedDict: continue
        for start,end,line in bedDict[contig]:
            for chr,map,reverse,supp in maps:
                startRef, endRef = seqLift.liftQuery2Ref(\
                    seqDict[contig], map, reverse, start, end)
                if startRef != 0 and endRef != 0:
                    temp = f"{chr}\t{startRef}\t{endRef}\t{reverse}\t{supp}\t"
                    out.write(temp+line)
                    break
    out.close()
    logging.info('Lifting each read finish')

logging.info('End of Program\n')
