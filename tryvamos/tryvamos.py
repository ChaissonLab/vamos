#!/usr/bin/env python3
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
parserDict = {}
# author and version info
usage = ' <mode> <options> \n'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
This program provides a tool set for various downstream analysis using TR
annotation output from the Vamos software.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers(dest='command')

######################### combineVCF mode #########################
parser_combineVCF = subparsers.add_parser('combineVCF', help=
'''quickFeature command generates diploid multi-sample
vcf from given haploid or diploid single sample vcfs.
Each line of the input csv file can be a single diploid
vcf or two haploid vcf as "hap1.vcf,hap2.vcf".
''')
combineVCFPosList = ['inVCFs', 'outVCF']
combineVCFOptList = []
parserDict['combineVCF'] = [combineVCFPosList, combineVCFOptList]
# positional arguments
parser_combineVCF.add_argument(combineVCFPosList[0], type=str, \
    help='string\tinput list of vcfs,  e.g. /in/samples.csv')
parser_combineVCF.add_argument(combineVCFPosList[1], type=str, \
    help='string\toutput combined vcf,  e.g. /out/combined.vcf')
# optional arguments


######################### quickFeature mode #########################
parser_quickFeature = subparsers.add_parser('quickFeature', help=
'''quickFeature command generates feature matrix for each
sample and alleles from the input vamos diploid vcf.
Currently 4 feature types are supported:
    annoLen: allele length in motif unit
    annoStr: allele by motif string
    topCount: allele by count of the most frequent motif
              (i.e., motif "0" in the annotation string)
    nt: allele by annotated nucleotide string
''')
quickFeaturePosList = ['inVCF', 'outFile']
quickFeatureOptList = ['feature', 'demographics', 'skipLoci']
parserDict['quickFeature'] = [quickFeaturePosList, quickFeatureOptList]
# positional arguments
parser_quickFeature.add_argument(quickFeaturePosList[0], type=str, \
    help='string\tinput combined diploid vamos vcf,  e.g. /in/samples.vcf')
parser_quickFeature.add_argument(quickFeaturePosList[1], type=str, \
    help='string\toutput tsv file,  e.g. /out/File.tsv')
# optional arguments
parser_quickFeature.add_argument('-f', '--'+quickFeatureOptList[0], \
    type=str, metavar='', default='annoLen', \
    choices=['annoLen','annoStr','topCount','nt'], \
    help='{annoLen,annoStr,topCount,nt}\tfeature type, default annoLen')
parser_quickFeature.add_argument('-d', '--'+quickFeatureOptList[1], \
    type=str, metavar='', default=None, \
    help='string\tdemographics info,  default None')
parser_quickFeature.add_argument('-s', '--'+quickFeatureOptList[2], \
    type=str, metavar='', default=None, \
    help='string\tlist of loci to skip (in bed format),  default None')


######################### waterfallPlot mode #########################
waterfallPlotPosList = ['inVCF', 'outDir']
waterfallPlotOptList = ['useLoci', 'sort']
parserDict['waterfallPlot'] = [waterfallPlotPosList, waterfallPlotOptList]
parser_waterfallPlot = subparsers.add_parser('waterfallPlot', help=
'''waterfallPlot command generates waterfall plots of all
loci or selected loci in the input vamos diploid vcf.
''')
# positional arguments
parser_waterfallPlot.add_argument(waterfallPlotPosList[0], type=str, \
    help='string\tinput combined diploid vamos vcf,  e.g. /in/samples.vcf')
parser_waterfallPlot.add_argument(waterfallPlotPosList[1], type=str, \
    help='string\toutput directory,  e.g. /out/Dir')
# optional arguments
parser_waterfallPlot.add_argument('-l', '--'+waterfallPlotOptList[0], 
    type=str, metavar='', default=None, \
    help='string\tinput bed of loci for plotting, default None')
parser_waterfallPlot.add_argument('-s', '--'+waterfallPlotOptList[1], \
    action='store_true', default=False, \
    help='bool\tshould the plot alleles be sorted, default False')


######################### testTwoPanels mode #########################
testTwoPanelsPosList = ['inVCF1', 'inVCF2', 'outFile']
testTwoPanelsOptList = ['skipLoci', 'testType', 'varCut', 'outPlotDir']
parserDict['testTwoPanels'] = [testTwoPanelsPosList, testTwoPanelsOptList]
parser_testTwoPanels = subparsers.add_parser('testTwoPanels', help=
'''testTwoPanels command performs statistical tests to
compare alleles for each TR locus on two given panel of
samples. Currently 2 test types are supported:
    tf: t-test and f-test of annotation length
    ks: ks-test of motif count distribution
''')
# positional arguments
parser_testTwoPanels.add_argument(testTwoPanelsPosList[0], type=str, \
    help='string\tinput vcf file 1,  e.g. /in/1.vcf')
parser_testTwoPanels.add_argument(testTwoPanelsPosList[1], type=str, \
    help='string\tinput vcf file 2,  e.g. /in/2.vcf')
parser_testTwoPanels.add_argument(testTwoPanelsPosList[2], type=str, \
    help='string\toutput tsv file,  e.g. /out/File.tsv')
# optional arguments
parser_testTwoPanels.add_argument('-s', '--'+testTwoPanelsOptList[0], \
    type=str, metavar='', default=None, \
    help='string\tlist of loci to skip (in bed format),  default None')
parser_testTwoPanels.add_argument('-t', '--'+testTwoPanelsOptList[1], \
    type=str, metavar='', default='tf', choices=['tf','ks'], \
    help='{tf,ks}\ttype of statistical test,  default tf')
parser_testTwoPanels.add_argument('-c', '--'+testTwoPanelsOptList[2], \
    type=float, metavar='', default=1, \
    help='float\tvariance lower bound for tested (tf only),  default 1')
parser_testTwoPanels.add_argument('-o', '--'+testTwoPanelsOptList[3], \
    type=str, metavar='', default=None, \
    help='string\tdir of histograms of sig loci (tf only),  default None')


######################### pairwiseCompare mode #########################
pairwiseComparePosList = ['inVCF1', 'inVCF2', 'outFile']
pairwiseCompareOptList = ['skipLoci']
parserDict['pairwiseCompare'] = [pairwiseComparePosList, pairwiseCompareOptList]
parser_pairwiseCompare = subparsers.add_parser('pairwiseCompare', help=
'''pairwiseCompare command generates allele comparison for
TR loci of given pair of single haplotype vamos vcfs.
''')
# positional arguments
parser_pairwiseCompare.add_argument(pairwiseComparePosList[0], type=str, \
    help='string\tinput vcf file 1,  e.g. /in/1.vcf')
parser_pairwiseCompare.add_argument(pairwiseComparePosList[1], type=str, \
    help='string\tinput vcf file 2,  e.g. /in/2.vcf')
parser_pairwiseCompare.add_argument(pairwiseComparePosList[2], type=str, \
    help='string\toutput txt file,  e.g. /out/File.txt')
# optional arguments
parser_pairwiseCompare.add_argument('-s', '--'+pairwiseCompareOptList[0], \
    type=str, metavar='', default=None, \
    help='string\tlist of loci to skip (in bed format),  default None')


########## strict positional/optional arguments checking ##########
argsDict = vars(parser.parse_args())
posList, optList = parserDict[argsDict['command']]
logging.info('Parsing Input Arguements...')
logging.info(f'Required Argument - mode: {argsDict["command"]}')
for key, value in argsDict.items():
    if key in posList: logging.info(f'Required Argument - {key}: {value}')
    if key in optList: logging.info(f'Optional Argument - {key}: {value}')
    vars()[key] = value # assign values of arguments into global variables
logging.info('Parsing Input Arguements Completed\n')

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    ######################### combineVCF mode #########################
    if command == 'combineVCF':
        import lib.combineVCF as combineVCF

        logging.info(f'Combing vcfs...')
        combineVCF.combineVCF(inVCFs, outVCF)

    ######################### quickFeature mode #########################
    if command == 'quickFeature':

        import lib.quickFeature as quickFeature

        logging.info(f'Reading input vcf to generate feature matrix...')
        quickFeature.writeFeature(inVCF, feature, outFile, \
                                  skipLoci, demographics)

    ######################### waterfallPlot mode #########################
    elif command == 'waterfallPlot':
        import lib.waterfallPlot as waterfallPlot

        logging.info(f'Gnerating waterfall plot...')
        waterfallPlot.plot(inVCF, outDir, useLoci, sort)


    ######################### testTwoPanels mode #########################
    elif command == 'testTwoPanels':

        import lib.testTwoPanels as testTwoPanels

        logging.info(f'Reading input vcfs for statistical tests...')
        testTwoPanels.test(inVCF1, inVCF2, testType, skipLoci, \
                           varCut, outFile, outPlotDir)


    ######################### pairwiseCompare mode #########################
    elif command == 'pairwiseCompare':

        import lib.pairwiseCompare as pairwiseCompare

        logging.info(f'Reading input vcfs for pairwise comparison...')
        pairwiseCompare.pairwiseCompare(inVCF1, inVCF2, skipLoci, outFile)

    logging.info('End of Program\n')

