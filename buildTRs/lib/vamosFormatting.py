#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
import logging
import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inMani', 'outPre']
optList = ['genomicAnnotations']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
Aggregated original motif sets and conpressed efficient motif sets may be stored
as separate files. The program reads in all such oriMotif/effMotif files and
formats the entries as oriMotif/effMotif sets for vamos analysis.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput manifect csv (format tag,file),  e.g. /in/Mani.csv')
parser.add_argument(posList[1], type=str, \
    help='string\toutput file prefix,  e.g. /out/Pre')
# optional arguments
parser.add_argument('-a', '--'+optList[0], type=str, metavar='', default=None, \
    help='string\tcsv manifest of extra annotations for tagging,  default None')

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

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":


    logging.info('Processing oriMotif files...')
    with open(inMani) as f: files = [ l.strip().split(',') for l in f ]
    tagsE = [ f[0] for f in files if f[0] != 'oriMotifs' ]

    oDict, taggingDict = {}, {}
    for tag,file in files:
        if tag == 'oriMotifs':
            with open(file) as f:
                for i,line in enumerate(f):
                    if i == 0: continue # skip header
                    chr,start,end,motifs,counts,cons,op1,op2 = \
                            line.strip().split('\t')
                    size = len(cons)
                    cat = 'STR' if size <= 6 else 'VNTR'
                    if chr not in oDict: oDict[chr] = {}
                    oDict[chr][(start,end)] = [motifs,size,cat,op1,op2]
                    taggingDict[(chr,start,end)] = []

    # sort oDict by start
    tempDict = {}
    for chr,chrDict in oDict.items():
        temp = sorted(chrDict.items(), key=lambda x: int(x[0][0]))
        tempDict[chr] = {k:v for k,v in temp}
    oDict = tempDict
    logging.info('Processing oriMotif files finish')

    # tagging extra genomic annotations
    if genomicAnnotations:
        logging.info('Tagging extra genomic annotations...')
        annos = {}
        with open(genomicAnnotations) as f:
            for line in f:
                annoName,path = line.strip().split(',')
                annos[annoName] = path

        for annoName,path in annos.items():
            annoDict = {}
            with open(path) as f:
                for i,line in enumerate(f):
                    if i < 1: continue # skip header
                    chr,start,end = line.strip().split('\t')[:3]
                    annoDict[(chr,start,end)] = f'{chr}_{start}_{end}'

            # tagging TR loci with the extra annotations
            logging.info(f'Tagging extra genomic annotations: {annoName}')
            taggingDict = bedlib.bedIntersectTag(taggingDict, annoDict)

        logging.info('Tagging extra genomic annotations finish')

    # output oriMotifs
    outO = open(f'{outPre}.oriMotifs.tsv', 'w')
    for chr,chrDict in oDict.items():
        for (start,end),(oriMotifs,size,cat,op1,op2) in chrDict.items():
            if genomicAnnotations:
                annoNames = [ annoName for annoName,path in annos.items() ]
                annoTags = taggingDict[(chr,start,end)]
                annoTags = [ annoNames[i] if t != 'NULL' else '.' \
                            for i,t in enumerate(annoTags) ]
            else:
                annoTags = []
            outList = [chr,start,end,oriMotifs,'v-2.1',cat,size,op1,op2]
            outO.write('\t'.join(map(str, outList+annoTags)) +'\n')
    outO.close()

    # output effMotifs
    for tag in tagsE:
        logging.info(f'processing {tag}...')
        eDict = {}
        for t,file in files:
            if t != tag: continue
            with open(file) as f:
                for i,line in enumerate(f):
                    if i == 0: continue
                    chr,start,end,effMotifs = line.strip().split('\t')[:4]
                    size = oDict[chr][(start,end)][1]
                    cat = oDict[chr][(start,end)][2]
                    if chr not in eDict: eDict[chr] = {}
                    eDict[chr][(start,end)] = [effMotifs,size,cat]

        outE = open(f'{outPre}.{tag}.tsv', 'w')
        for chr,chrDict in oDict.items():
            for (start,end),(oriMotifs,size,cat,op1,op2) in chrDict.items():
                effMotifs = eDict[chr][(start,end)][0]
                if genomicAnnotations:
                    annoNames = [ annoName for annoName,path in annos.items() ]
                    annoTags = taggingDict[(chr,start,end)]
                    annoTags = taggingDict[(chr,start,end)]
                    annoTags = [ annoNames[i] if t != 'NULL' else '.' \
                                for i,t in enumerate(annoTags) ]
                else:
                    annoTags = []
                outList = [chr,start,end,effMotifs,'v-2.1',cat,size,op1,op2]
                outE.write('\t'.join(map(str, outList+annoTags)) +'\n')
        outE.close()
        logging.info(f'processing {tag} finish')


logging.info('End of Program\n')

