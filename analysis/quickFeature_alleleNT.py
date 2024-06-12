#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import datetime
import logging

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)


#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inVCF', 'outFile']
optList = []
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2024-04-15)'
description = '''\nDescription:
This program parses the vamos combined diploid vcf and output annotation alleles
by motif nucleotide sequences for each haplotype.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput vcf file,  e.g. /in/samples.vcf')
parser.add_argument(posList[1], type=str, \
    help='string\toutput directory,  e.g. /out/File')
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

class TR:
    """class for one TR locus.

    Attributes:
        chr (str): chromosome of TR.
        start (str): start coordinate of TR.
        end (str): end coordinate of TR.
        motifsFull (dict[str:int]): indexing table of the full list of motifs
        motifsUsed (dict[str:int]): indexing table of the used list of motifs
        annosByUsed (dict[str:list[str]]): annotations of all samples (by index
            of the used motif list, unannotated allele as ['.'], full deletion
            allele as [])
        annosByFull (dict[str:list[str]]): annotations of all samples (by index
            of the full motif list, unannotated allele as ['.'], full deletion
            allele as [])
        constant (bool): if the locus is constant across all samples
        space (bool): if annotation of any sample is extended to the max length
        edMat (dict[tuple[(str,str)]:int]): edit distance (of annotation string)
            matrix for all samples

    """

    def __init__(self):
        """the __init__ method
        """

        # parsed in method "parseCombinedVCFOneLine"
        self.chr = ''
        self.start = ''
        self.end = ''
        self.motifsUsed = {}
        self.annosByUsed = {}
        self.annosByFull = {}
        self.constant = False
        # parsed in method "readMotifsFull"
        self.motifsFull = {}
        # changed in method "appendLength"
        self.space = False
        # assigned in method "getEditDistance"
        self.edMat = {}


    def readMotifsFull(self, motifsFull:list[str]):
        """Encode the full motif list with numbering

        Args:
            motifs (list[str]): input list of full motifs
        """
        for i,m in enumerate(motifsFull): self.motifsFull[m] = i


    def parseCombinedVCFOneLine(self, samples:list[str], line:str):
        """parse one entry line from a combined vamos vcf

        Args:
            samples (list[str]): ordered list of all sample names for gt parsing
            line (str): one line from the combined vamos vcf
        """

        fields = line.strip().split()
        self.chr,self.start,_,_,_,_,_,info,_ = fields[:9]
        end,ru,_,annos = info.split(';')
        self.end = end.split('=')[1]
        ru = ru.split('=')[1].split(',')
        # encode the list of motifs used for annotation with numbering
        for i,m in enumerate(ru): self.motifsUsed[m] = i

        gts = [ gt.split('/') for gt in fields[9:] ]
        alleles = annos.split('=')[1].split(',')

        # check if the locus is constant over all samples
        if alleles.count(alleles[0]) == len(alleles): self.constant = True

        # parse the h1/h2 allele
        for i,sample in enumerate(samples):

            if gts[i][0] == '.':
                temp = ['.']
            elif gts[i][0] == 'DEL':
                temp = []
            else:
                temp = alleles[int(gts[i][0])-1].split('-')

            self.annosByUsed[sample+'_h1'] = temp
            if self.annosByFull:
                if temp != ['.']:
                    temp = [ self.motifsFull[ru[int(t)]] for t in temp ]
                self.annosByFull[sample+'_h1'] = temp

            if gts[i][1] == '.':
                temp = ['.']
            elif gts[i][1] == 'DEL':
                temp = []
            else:
                temp = alleles[int(gts[i][1])-1].split('-')

            self.annosByUsed[sample+'_h2'] = temp
            if self.annosByFull:
                if temp != ['.']:
                    temp = [ self.motifsFull[ru[int(t)]] for t in temp ]
                self.annosByFull[sample+'_h2'] = temp

    def getEditDistance(self):
        """calculate pairwise edit distance of all annotation strings
        apply ascii encoding if <=256 unique motifs and edlib for speeding up
        """

        checkDict = {}
        for i,s in enumerate(self.annosByUsed.keys()):
            for j,t in enumerate(self.annosByUsed.keys()):
                if i >= j: continue
                anno1, anno2 = self.annosByUsed[s], self.annosByUsed[t]
                temp1, temp2 = '-'.join(anno1), '-'.join(anno2)
                if anno1 in [['.'],[]] or anno2 in [['.'],[]]:
                    ed = 'NULL'
                else:
                    if (temp1, temp2) in checkDict:
                        ed = checkDict[(temp1, temp2)]
                    else:
                        if len(self.motifsUsed) <= 256:
                            anno1 = [chr(int(a)) for a in anno1]
                            anno2 = [chr(int(a)) for a in anno2]
                            ed = edlib.align(anno1, anno2)['editDistance']
                        else:
                            ed = seqAlign.alignGlobal(anno1,anno2,True,0,1,1)
                        checkDict[(temp1, temp2)] = ed
                        checkDict[(temp2, temp1)] = ed
                self.edMat[(s,t)] = ed
                self.edMat[(t,s)] = ed


def parseCombinedVCF(vcf, outFile):

    out = open(outFile, 'w')

    with open(vcf) as f:
        for i,line in enumerate(f,start=1):
            if line.startswith('##'): continue
            if line.startswith('#'):
                samples = line.strip().split()[9:]
                temp = []
                for s in samples: temp += s.split('/')
                out.write('\t'.join(['#chr','start','end']+temp)+'\n')
                continue

            tr = TR()
            tr.parseCombinedVCFOneLine(samples, line)

            # skip unwanted loci
            #if not (tr.chr == 'chr1' and tr.start == '11636'): continue

            alleles = []
            for s,a in tr.annosByUsed.items():
                if a == ['.']:
                    alleles.append('NA')
                else:
                    motifsDict = {v:k for k,v in tr.motifsUsed.items()}
                    alleleNT = '-'.join( [motifsDict[int(i)] for i in a] )
                    alleles.append(alleleNT)

            out.write('\t'.join([tr.chr,tr.start,tr.end]+alleles) +'\n')

            if i % 100000 == 0: logging.info(f'{i} of lines processed.')

    out.close()


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    logging.info('Reading input vcf...')
    parseCombinedVCF(inVCF, outFile)

logging.info('End of Program\n')

