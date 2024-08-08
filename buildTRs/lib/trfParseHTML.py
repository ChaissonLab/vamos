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
posList = ['inHTML', 'inBed', 'outFile']
optList = []
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
This program parses the trf txt.html & dat output to give a single ".tsv" file
(usually named as ".raw.tsv").
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput HTML,  e.g. /in/HTML')
parser.add_argument(posList[1], type=str, \
    help='string\tinput bed,  e.g. /in/bed')
parser.add_argument(posList[2], type=str, \
    help='string\toutput file,  e.g. /out/File')
# optional arguments


# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
logging.info('Parsing Input Arguements...')
for key, value in argsDict.items():
    if key in posList: logging.info('Required Argument - %s: %s' %(key,value))
    if key in optList: logging.info('Optional Argument - %s: %s' %(key,value))
    vars()[key] = value # assign values of arguments into global variables
logging.info('Parsing Input Arguements Completed\n')

#--------------------------------------------------------
# global variables and user-defined functions
#--------------------------------------------------------

def bedConfig(inBed):

    bedDict = {}
    with open(inBed) as f:
        for line in f:
            if line.startswith('@'):
                cluster = line.strip().replace('@', '')
                bedDict[cluster] = {}
                continue

            start = line.strip().split(' ')[0]
            end = line.strip().split(' ')[1]
            period = line.strip().split(' ')[2]
            consensus = line.strip().split(' ')[13]
            bedDict[cluster]['__'.join([cluster,start,end,period,consensus])] \
                = line.strip().split(' ')

    return bedDict


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    bedDict = bedConfig(inBed)

    oline, monomer, monomers, consensusTag = ['', [], [], 0]

    with open(inHTML) as f:
        for line in f:

            # sline: current line
            # oline: old line
            sline = line.strip()

            if sline.startswith('Sequence: '):
                cluster = sline.replace('Sequence: ', '')
            if sline.startswith('Period size: '):
                period = sline.split(' ')[2]
            if sline.startswith('Indices: '):
                start, end = sline.split(' ')[1].split('--')
                monomer, monomers = [[], []]
                oline = sline
                continue

            # sline being a consensus alignment line
            # (oline being a monomer alignment line)
            if len(sline.split(' ')) >= 2 and len(oline.split(' ')) >= 2:
                s, o = [sline.split(' '), oline.split(' ')]
                if re.match('\d+', s[0]) and re.match('\D+', s[1]) and \
                    re.match('\d+', o[0]) and re.match('\D+', o[1]):

                    if s[0] == '1': # a new monomer
                        # append last monomer
                        if monomer != []: monomers += monomer
                        # start new monomer
                        monomer = ' '.join(o[1:]).replace('-','').split(' ')

                    else: # extending existing monomer
                        temp = ' '.join(o[1:]).replace('-','').split(' ')

                        if len(monomer) == 1: # a compressed trf entry
                            #print(monomer); print(temp)
                            monomer = [ monomer[0] + temp[0] ]
                        else: # a normal trf entry
                            monomer += temp

            # sline being a consensus sequence line
            if oline.startswith('Consensus pattern'): # start new consensus
                monomers += monomer # append the last monomer of this trf seq
                consensusTag = 1
                consensus = sline
                oline = sline
                continue
            if sline.startswith('Left flanking sequence'): # end this consensus
                # handle single-base consensus
                if len(consensus) == 1:
                    temp = []
                    for m in monomers: temp += list(m)
                    monomers = temp
                consensusTag = 0
                temp = '__'.join( [cluster,start,end,period,consensus] )
                bedDict[cluster][temp].append(monomers)
            if consensusTag == 1: consensus += sline # extend this consensus

            oline = sline

    # output
    out = open(outFile, 'w')
    for cluster,value in bedDict.items():
        for key,data in value.items():
            temp = [cluster] + data[:-3] + [','.join(data[-1])]
            out.write('\t'.join( temp ) + '\n')
    out.close()


logging.info('End of Program\n')

