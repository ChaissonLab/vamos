#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
import logging
import numpy as np
import statistics

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)


#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['prefix', 'motifTable', 'outFile']
optList = []
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2023-10-01)'
description = '\nDescription: The program configs input format for argweaver'

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tprefix of vamos feature table,  e.g. /in/feature/prefix')
parser.add_argument(posList[1], type=str, \
    help='string\tvamos motif config file,  e.g. /vamos/effMotifs-0.1.tsv')
parser.add_argument(posList[2], type=str, \
    help='string\toutput file,  e.g. /out/file.tsv')
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
            'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17', \
            'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

def zygosity(alleles):

    countDict, homo = {}, 0
    for a in alleles:
        if a not in countDict: countDict[a] = 0
        countDict[a] += 1
    for allele,count in countDict.items():
        f = count/len(alleles)
        homo += f * f
    hete = 1 - homo

    return homo,hete

def getStats(alleles):

    # filter uncovered sample
    alleles = [ a for a in alleles if a != 'NA' ]
    lengths = [ len(a.split('-')) for a in alleles ]
    cover = len(alleles)

    # compute zygosity
    homoCom, heteCom = zygosity(alleles)
    homoLen, heteLen = zygosity(lengths)

    # compute variability
    alleles = [ a for a in alleles if a != 'NA' ]
    lengths = [ len(a.split('-')) for a in alleles ]
    cover = len(alleles)
    alleleCom = len(set(alleles))
    alleleLen = len(set(lengths))

    # compute std
    if len(lengths) == 0:
        stdLen = -1
    elif len(lengths) == 1:
        stdLen = 0
    else:
        stdLen = statistics.stdev(lengths)

    # outlier test
    numOutlier = 0
    if lengths != []:
        Q1 = np.percentile(lengths, 25)
        Q3 = np.percentile(lengths, 75)
        IQR = Q3 - Q1
        lower = Q1 - 5 * IQR
        upper = Q3 + 5 * IQR
        outliers = []
        for l in lengths:
            if (l < lower) or (l > upper):
                outliers.append(l)
        numOutlier = len(outliers)

    return cover,stdLen,numOutlier,homoLen,heteLen,homoCom,heteCom,alleleLen,alleleCom

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    locusPeriod, locusType = {}, {}
    with open(motifTable) as f:
        for line in f:
            chr,start,end,motifs,source,ty,period = line.strip().split()[:7]
            locusPeriod[(chr,start,end)] = period
            locusType[(chr,start,end)] = ty
    logging.info(f'reading motif done')

    featureTable = f'../features/{prefix}.tsv'

    out = open(outFile, 'w')
    header = ['chr','start','end','period','type','coverAll','stdLenAll','numOutlierAll','homoLenAll','heteLenAll',\
              'homoComAll','heteComAll','lenAll','comAll']
    batchDict = {'APR-phase1':[],'CPC-phase1':[],'HGSVC-phase3':[],'HPRC-phase2':[],'JSA-phase1':[]}
    for batch in batchDict.keys():
        header += [f'cover_{batch}',f'stdLen_{batch}',f'numOutlier_{batch}',f'homoLen_{batch}',f'heteLen_{batch}',\
                   f'homoCom_{batch}',f'heteCom_{batch}',f'len_{batch}',f'com_{batch}']
    out.write('\t'.join(header) +'\n')

    with open(featureTable) as f:
        for i,line in enumerate(f,start=1):
            if i % 100000 == 0: logging.info(f'{i} loci handled')
            if i == 1:
                samples = line.strip().split('\t')[3:]
                batches = [ s.split('_')[0] for s in samples ]
                continue
            if line.startswith('#'): continue
            chr,start,end = line.strip().split('\t')[:3]
            allelesAll = line.strip().split('\t')[3:]

            # get batch-level alleles
            batchAlleleDict = {'APR-phase1':[],'CPC-phase1':[],'HGSVC-phase3':[],'HPRC-phase2':[],'JSA-phase1':[]}
            for j,sample in enumerate(samples):
                if batches[j] not in batchAlleleDict: continue
                batchAlleleDict[batches[j]].append(allelesAll[j])

            # compute overall stats (cover,zygosity,variability...)
            alleles = [ a for a in allelesAll if a != 'NA' ]
            stats = [ s for s in getStats(alleles) ]

            # compute batch level stats (cover,zygosity,variability...)
            for batch,batchAlleles in batchAlleleDict.items():
                batchAlleles = [ a for a in batchAlleles if a != 'NA' ]
                stats += [ s for s in getStats(batchAlleles) ]

            temp = [ chr,start,end,locusPeriod[(chr,start,end)],locusType[(chr,start,end)] ]
            temp += [ round(s,4) for s in stats ]
            out.write('\t'.join(map(str,temp)) +'\n')
    out.close()

logging.info('End of Program\n')
