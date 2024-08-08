#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
import logging
import random as rm

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inFile', 'splits', 'outPre', 'outSuff']
optList = ['skip', 'random']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
The final oriMotif set from "oriMotifFinal.py" may contain huge number of TR
entries. This program evenly splits these entries into individual batches for
more efficient processing of the effMotif selection step.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput file,  e.g. /in/File')
parser.add_argument(posList[1], type=int, \
    help='integer\tnumber of splits,  e.g. 100')
parser.add_argument(posList[2], type=str, \
    help='string\toutput prefix,  e.g. /out/Pre')
parser.add_argument(posList[3], type=str, \
    help='string\toutput suffix,  e.g. /out/Suff')
# optional arguments
parser.add_argument('-s', '--'+optList[0], type=int, metavar='', \
    default=0, help='integer\tnumber of lines to skip,  default 0')
parser.add_argument('-r', '--'+optList[1], action='store_true', default=False, \
    help='bool\trandomize file entries, default False')


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


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    # read input file
    with open(inFile) as f: all = [l for i,l in enumerate(f,start=1) if i>skip]

    if random:
        rm.shuffle(all)
        logging.info('Entries randomized.')

    # determine splits
    groupDict = { g+1:[] for g in range(splits) }

    g = 1
    for i in range(len(all)):
        groupDict[g].append(i)
        g += 1
        if g > splits: g = 1

    # split
    for g,index in groupDict.items():
        out = open('%s.%s.%s' %(outPre, g, outSuff), 'w')
        for i in index: out.write(all[i])
        out.close()


logging.info('End of Program\n')

