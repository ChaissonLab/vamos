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
import seqBasic
import bedlib
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
optList = ['reCom', 'startAdjust', 'endAdjust', 'skipEmpty']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
The final TR boundaries obtained from "boundaryDecomposeLiftRefineAggregate.py"
are given on reference. This program returns the corresponding sequences of
these regions on other sequences (e.g., assembly) based on the alignment bam.
Input bam must use the reference of the given bed regions.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput bed regions,  e.g. /in/bed')
parser.add_argument(posList[1], type=str, \
    help='string\tinput assembly fasta,  e.g. /in/fa')
parser.add_argument(posList[2], type=str, \
    help='string\tinput assembly bam,  e.g. /in/bam')
parser.add_argument(posList[3], type=str, \
    help='string\toutput fasta file,  e.g. /out/file.fasta')
# optional arguments
parser.add_argument('-r', '--'+optList[0], action='store_true', default=True, \
    help='bool\tadjust for reverse complement, default True')
parser.add_argument('-s', '--'+optList[1], type=int, metavar='', \
    default=0, help='int\tboundary adjust for starting coordinate,  default 0')
parser.add_argument('-e', '--'+optList[2], type=int, metavar='', \
    default=0, help='int\tboundary adjust for ending coordinate,  default 0')
parser.add_argument('-d', '--'+optList[3], action='store_true', default=True, \
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
    bedDict, inChrList = {}, []
    with open(inBed) as f:
        for line in f:
            chr,start,end = line.strip().split('\t')[:3]
            # skip alternative chromosomes
            if chr not in chrList: continue
            if chr not in inChrList: inChrList.append(chr)
            bedDict[(chr,start,end)] = '_'.join([chr,start,end])
    logging.info('Reading input files finish')

    logging.info('Initial parse of each read...')
    # parse each read once to get the aligned ref start/end coordinates
    readsDict, counter = {}, 0
    for read in samfile.fetch():

        counter += 1
        if counter % 100 == 0: logging.info(f'{counter} reads processed')

        # skip alternative chrs, secondary alignments, alignments with <60 mapQ
        if read.reference_name not in chrList: continue
        if read.is_secondary: continue
        if read.mapping_quality < 60: continue
        # skip reads that doesn't have any bed to map
        if read.reference_name not in inChrList: continue

        id = read.query_name
        chr = read.reference_name
        cigar = read.cigarstring
        startAln = read.reference_start
        reverse = read.is_reverse
        supp = read.is_supplementary

        maps = seqLift.liftOnePass(cigar, startAln, refAsKey=True) # parse cigar
        endAln = max(list(maps.keys()))
        readsDict[(chr,startAln,endAln)] = [id, cigar, startAln, reverse, supp]

    samfile.close() # close bam file to release memory
    logging.info('Initial parse of each read finish')

    # group beds into reads
    logging.info('Mapping ref beds to each read...')
    readsDict = bedlib.bedIntersectTag(readsDict, bedDict)
    logging.info('Mapping ref beds to each read finish')

    # parse each read to get the liftover results for each ref bed
    priDict, suppDict = {}, {}
    for coor,(id,cigar,startAln,reverse,supp,beds) in readsDict.items():

        # skip reads that does not cover any bed
        if beds == 'NULL': continue
        beds = [ b.split('_') for b in beds.split(',') ]

        # parse cigar and liftover each bed covered by the read
        maps = seqLift.liftOnePass(cigar, startAln, refAsKey=True)
        for (chr,start,end) in beds:

            # skip bed that already has a primary liftover
            if (chr,start,end) in priDict: continue

            startQuery, endQuery = seqLift.liftRef2Query(seqDict[id], \
                    maps, reverse, int(start), int(end))
            seq = seqDict[id][(startQuery+startAdjust):(endQuery+endAdjust)]

            # skip empty liftovers
            if seq == '' and skipEmpty: continue
            # reverse complement extracted sequence if instructed
            if reverse and reCom: seq = seqBasic.reComDNA(seq)

            if supp:
                suppDict[(chr,start,end)] = [id,startQuery,endQuery,reverse,seq]
            else:
                priDict[(chr,start,end)] = [id,startQuery,endQuery,reverse,seq]

            # logging the number of processed beds
            temp = len(priDict) + len(suppDict)
            if temp % 50000 == 0: logging.info(f'{temp} beds handled')

    out = open(outFile, 'w')
    # get extracted sequences and output results
    for (chr,start,end),_ in bedDict.items():
        if (chr,start,end) in priDict:
            id, startQ, endQ, reverse, seq = priDict[(chr,start,end)]
            reverse = '-1' if reverse else '1'
            name = f'>{id}_{startQ}-{endQ};{chr}_{start}-{end};{reverse}'
            out.write(f'{name}\n{seq}\n')
        if (chr,start,end) in suppDict:
            id, startQ, endQ, reverse, seq = suppDict[(chr,start,end)]
            reverse = '-1' if reverse else '1'
            name = f'>{id}_{startQ}-{endQ};{chr}_{start}-{end};{reverse}'
            out.write(f'{name}\n{seq}\n')
    out.close()

logging.info('End of Program\n')

