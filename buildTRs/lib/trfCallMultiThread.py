#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
import logging

from Bio import SeqIO
from multiprocessing import Pool

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inFA', 'trf', 'parser', 'outPre']
optList = []
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
This program uses trf to call repeats from VNTRs.
Sequence entries from "inFA" is processed by multi-threading with a separate
path named after the sequence index (ordered as in the "inFA" file) under the
directory of "outPre".
Each sequence from "inFA" is written into the temp file ".fa" as trf input. trf
generates tandem repeat calls into temp file(s) ".html" and a final file ".dat".
Contents of individual ".html" files are concatenated into ".aln". ".html" and
".dat" files are parsed by "trfParseHTML.py" to give ".raw.tsv" files.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput fa,  e.g. /in/fa')
parser.add_argument(posList[1], type=str, \
    help='string\tpath of the trf,  e.g. /path/trf')
parser.add_argument(posList[2], type=str, \
    help='string\thtml parser,  e.g. /tool/trfParseHTML.py')
parser.add_argument(posList[3], type=str, \
    help='string\toutput prefix,  e.g. /out/Pre')
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

def call(key):

    logging.info(f'Handling fasta entry: {key}')

    # .fa file containing single VNTR entry for trf calling
    #outDir = os.path.dirname(outPre)
    dat =  f'{outPre}.{key}/{key}.dat'
    fa = f'{outPre}.{key}/{key}.fa'
    aln = f'{outPre}.{key}/{key}.aln'
    raw = f'{outPre}.{key}/{key}.raw.tsv'

    os.system(f'mkdir {outPre}.{key}')
    os.system(f'touch {aln} && rm {aln} && touch {aln}')

    out = open(fa, 'w')
    out.write('\n'.join( ['>'+contigs[key], seqDict[contigs[key]]] ))
    out.close()

    # trf calling
    os.system(f'{trf} {fa} 2 7 7 80 10 50 500 -f -m -ngs > {dat}')
    os.system(f'cat {key}.fa.*.txt.html >> {aln}')

    # clean up intermediate files
    os.system(f'rm -rf {fa} {key}.fa.*html {key}.fa.*mask')

    # parse the html and .dat output
    os.system(f'python {parser} {aln} {dat} {raw}')

    logging.info(f'Finished fasta entry: {key}')


def multi():

    keys = range(len(contigs))

    with Pool(16) as p:
        p.map(call, keys)


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    # read input fasta
    seqDict = {s.id:str(s.seq) for s in SeqIO.parse(inFA, 'fasta')}
    contigs = list(seqDict.keys())
    logging.info('Number of contigs: %s' %(len(contigs)))

    # multi-t
    multi()


logging.info('End of Program\n')

