#!/usr/bin/env python3
import os
import sys
import argparse
import re
import datetime
import logging

import edlib
import seqAlign


logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)


#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inFile', 'outFile']
optList = ['lineRange']
sys.stdout.write('\n'); sys.stdout.flush()
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
The program generates bam coverage statistics.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput annotation feature file.  e.g. /in/sample.annoStr.tsv')
parser.add_argument(posList[1], type=str, \
    help='string\toutput block change stats.  e.g. /out/file.tsv')
# optional arguments
parser.add_argument('-r', '--'+optList[0],
    type=str, metavar='string', default='0',
    help='range of lines to process ("," seperated int), default all')

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

def findCandidatePattern(seq,k=3):

    if len(seq) < k: return([],[])
    posDict, patterns, patternRanges = {}, [], []

    # obtain kmer ending position table
    for start in range(len(seq)-k):
        end = start + k
        kmer = tuple(seq[start:end])
        if kmer not in posDict: posDict[kmer] = []
        posDict[kmer].append(end)

    # obtain pattern list from adjacent kmer matches
    for kmer,ends in posDict.items():

        # skip when the kmer has only two occurrence, this means the block duplication probably happens only once
        if len(ends) < 2: continue

        # test if a kmer match gives promising duplication pattern
        for i in range(1,len(ends)):

            # get the candidate pattern size "d"
            endPrevious, endCurrent = ends[i-1], ends[i]
            d = endCurrent - endPrevious

            # skip cases whose pattern size < k, this happens with kmer matches in nested regions
            if d < k: continue
            # skip edge cases whose pattern extends beyond the sequence
            if endPrevious-d < 0: continue

            # get potential previous and current duplication patterns
            patternPrevious = seq[endPrevious-d:endPrevious]
            patternCurrent = seq[endPrevious:endCurrent]

            # skip homopolymer stretch
            if len(set(patternPrevious)) == 1: continue
            if len(set(patternCurrent)) == 1: continue

            # promissing pattern if the two patterns are similar enough
            # use edlib for cases with <= 256 motifs
            if max([int(s) for s in patternPrevious+patternCurrent]) <= 255:
                anno1 = [chr(int(a)) for a in patternPrevious]
                anno2 = [chr(int(a)) for a in patternCurrent]
                ed = edlib.align(anno1, anno2)['editDistance']
            # use self-written function for cases with > 256 motifs
            else:
                ed = seqAlign.alignGlobal(patternPrevious, patternCurrent, distance=True, match=0, mismatch=1, indel=1)

            # decide as valid block duplication if the two candidate patterns differ by < 20%
            if ed < 0.2 * d:
                patterns.append('-'.join(patternPrevious))
                patterns.append('-'.join(patternCurrent))
                patternRanges.append((endPrevious-d,endPrevious-1))
                patternRanges.append((endPrevious,endCurrent-1))

    return patterns, patternRanges

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    with open(inFile) as f: lineTotal = len(f.readlines())

    if lineRange == '0':
        beginning,ending = 1, lineTotal
    else:
        beginning,ending = [int(n) for n in lineRange.split(',')]

    out = open(outFile, 'w')
    numLoci = 0
    with open(inFile) as f:
        for index,line in enumerate(f,start=1):
            if index == 1:
                samples = line.strip().split()[3:]
                out.write('\t'.join(['chr','start','end','annotated','block','patterns','counts','nAffected'] + samples) +'\n')
                continue
            if line.startswith('#'): continue
            # skip lines not in range of processing
            if index < beginning or index > ending:
                continue
            numLoci += 1
            if numLoci % 1000 == 0:
                logging.info(f'number of loci processed: {numLoci}')
                logging.info(f'number of line passed: {index}')

            fields = line.strip().split()
            chrom,startChr,endChr = fields[:3]
            alleles = fields[3:]

            # skip locus with <20 annotated allele
            annotated = len([ a for a in alleles if a != 'NA' ])
            if annotated < 20: continue

            # run checker on each allele
            allelePatternDict = {}
            for i,allele in enumerate(alleles):

                if allele == 'NA': continue
                if allele not in allelePatternDict:

                    allelePatternDict[allele] = {}
                    patterns, patternRanges = findCandidatePattern(allele.split('-'))
                    if patterns:
                        # pattern range table
                        allelePatternDict[allele] = { patternRanges[j]:p for j,p in enumerate(patterns) }

            # calculate count of each pattern in all samples by all occurrence (count1) or unique occurrence (count2)
            # countDict[pattern] = [count1,count2]
            countDict = {}
            for allele,v in allelePatternDict.items():
                alleleN = alleles.count(allele)
                for pos,pattern in v.items():
                    if pattern not in countDict: countDict[pattern] = [0,0]
                    countDict[pattern][0] += alleleN
                for pattern in set(list(v.values())):
                    countDict[pattern][1] += alleleN
            # total number of genomes that are block positve
            totalAffected = 0
            for allele in alleles:
                if allele == 'NA': continue
                if allelePatternDict[allele]: totalAffected += 1

            # sort pattern by counts in alleles
            countDict = dict(sorted(countDict.items(), key=lambda item: item[1][1], reverse=True))

            # overall decision on all alleles
            block, allelePatternFinalDict = 'False', {}
            for pattern,(count1,count2) in countDict.items():
                if count2 >= 10:
                    block = 'True'
                    break

            # get final list of patterns for each sample
            for allele,v in allelePatternDict.items():
                allelePatternFinalDict[allele] = []
                v = dict(sorted(v.items(), key=lambda item: item[0][0]))
                for (start,end),pattern in v.items():
                    if not allelePatternFinalDict[allele]:
                        allelePatternFinalDict[allele].append(pattern)
                    else:
                        if overlap(startLast, endLast, start, end) > 0:
                            continue
                        else:
                            allelePatternFinalDict[allele].append(pattern)
                    startLast, endLast = start, end

            # output overall block info for this locus
            temp1, temp2 = [], []
            for i,(pattern,(count1,count2)) in enumerate(countDict.items()):
                #if i == 3: break
                temp1.append(pattern)
                temp2.append(str(count2))
            temp1 = ','.join(temp1)
            temp2 = ','.join(temp2)
            outList = [chrom,startChr,endChr,annotated,block,temp1,temp2,totalAffected]
            # output the final (non-overlapping) pattern for each sample of this locus
            for i,sample in enumerate(samples):
                allele = alleles[i]
                if allele != 'NA':
                    patterns = ','.join(allelePatternFinalDict[allele])
                else:
                    patterns = 'NA'
                outList += [patterns]
            outList = [ o if o != '' else 'NA' for o in outList]
            out.write('\t'.join(map(str,outList)) +'\n')
    out.close()

logging.info('End of Program\n')
