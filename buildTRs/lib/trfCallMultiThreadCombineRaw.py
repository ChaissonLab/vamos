#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
import logging

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['trfs', 'outFile']
optList = ['cut']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
Each fasta entry was called separately by "trfCallMulti.py" and processed as
one ".raw.tsv" result file. This program combines all input ".raw.tsv" files
into one single ".tsv" file.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput list of trf results,  e.g. /in/trf.raw.list')
parser.add_argument(posList[1], type=str, \
    help='string\toutput file,  e.g. /out/File')
# optional arguments
parser.add_argument('-c', '--'+optList[0], type=int, metavar='', \
    default=50000, help='integer\tonly keep TRs <= length cut,  default 50000')

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

with open(trfs) as f: files = [ l.strip() for l in f ]


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    # read input trf raw.tsv files
    with open(trfs) as f: files = [ l.strip() for l in f ]

    lociDict = {}
    out = open(outFile, 'w')
    # process
    for file in files:
        with open(file) as f:
            for line in f:
                contig,start,end = line.strip().split()[:3]

                head = contig.split(':')[1].split('-')[0]
                contig = contig.split(':')[0]
                start = int(head) + int(start) - 1
                end = int(head) + int(end) - 1

                if end - start > cut: continue
                #lociDict[(contig,start,end)] = line
                temp = line.strip().split()[3:]
                out.write('\t'.join([contig,str(start),str(end)]+temp) +'\n')

    out.close()


logging.info('End of Program\n')

