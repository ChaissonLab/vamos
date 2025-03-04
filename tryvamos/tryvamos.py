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
version = 'Version: 1.2.0.0'
description = '''\nDescription:
The tryvamos program provides a tool set for various downstream analysis using TR annotation output from the Vamos 
software.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers(dest='command')

######################### combineVCF mode #########################


parser_combineVCF = subparsers.add_parser('combineVCF', description=
'''combineVCF:
This command generates diploid multi-sample vcf from given haploid or diploid single sample vcfs.
Each line of the input csv file can be a single diploid vcf as "dip.vcf" or two comma delimited haploid vcfs as 
"hap1.vcf,hap2.vcf".
The best chromosome orders in the combined vcf are determined from the chromosome orders of all input vcf (as ordered
in the vcf headers).
For each locus, allele index is ordered by the order it appears in all samples.
Samples are ordered as in the input csv file.
''', formatter_class=argparse.RawTextHelpFormatter)
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
parser_quickFeature = subparsers.add_parser('quickFeature', description=
'''quickFeature:
This command generates feature matrix for each haplotype/sample and alleles from the input vamos diploid vcf.
Currently 4 feature types are supported:
    annoLen:    allele length in motif unit
    annoStr:    allele by motif string
    topCount:   allele by count of the most frequent motif (i.e., motif "0" in
                the annotation string)
    nt:         allele by annotated nucleotide string
''', formatter_class=argparse.RawTextHelpFormatter)
quickFeaturePosList = ['inVCF', 'outFile']
quickFeatureOptList = ['feature', 'byDip', 'demographics', 'skipLoci']
parserDict['quickFeature'] = [quickFeaturePosList, quickFeatureOptList]
# positional arguments
parser_quickFeature.add_argument(quickFeaturePosList[0], type=str,
    help='string\tinput combined diploid vamos vcf,  e.g. /in/samples.vcf')
parser_quickFeature.add_argument(quickFeaturePosList[1], type=str,
    help='string\toutput tsv file,  e.g. /out/File.tsv')
# optional arguments
parser_quickFeature.add_argument('-f', '--'+quickFeatureOptList[0],
    type=str, metavar='{annoLen,annoStr,topCount,nt}', default='annoLen',
    choices=['annoLen','annoStr','topCount','nt'],
    help='feature type, default annoLen')
parser_quickFeature.add_argument('-D', '--'+quickFeatureOptList[1],
    action='store_true',
    help='output features by diploid sample, homozygote as "feature/-"')
parser_quickFeature.add_argument('-d', '--'+quickFeatureOptList[2],
    type=str, metavar='string', default=None,
    help='demographics info,  default None')
parser_quickFeature.add_argument('-s', '--'+quickFeatureOptList[3],
    type=str, metavar='string', default=None,
    help='list of loci to skip (in bed format),  default None')


######################### waterfallPlot mode #########################
waterfallPlotPosList = ['inVCF', 'outDir']
waterfallPlotOptList = ['useLoci', 'width', 'height', 'ylabel', 'sort']
parserDict['waterfallPlot'] = [waterfallPlotPosList, waterfallPlotOptList]
parser_waterfallPlot = subparsers.add_parser('waterfallPlot', description=
'''waterfallPlot:
This command generates waterfall plots of all loci or selected loci in the input vamos diploid vcf.
''', formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser_waterfallPlot.add_argument(waterfallPlotPosList[0], type=str,
    help='string\tinput combined diploid vamos vcf,  e.g. /in/samples.vcf')
parser_waterfallPlot.add_argument(waterfallPlotPosList[1], type=str,
    help='string\toutput directory,  e.g. /out/Dir')
# optional arguments
parser_waterfallPlot.add_argument('-l', '--'+waterfallPlotOptList[0], 
    type=str, metavar='string', default=None,
    help='input bed of loci for plotting, default None')
parser_waterfallPlot.add_argument('-W', '--'+waterfallPlotOptList[1],
    type=float, metavar='float', default=None,
    help='use-specified plotting width to replace the program-determined one, default None')
parser_waterfallPlot.add_argument('-H', '--'+waterfallPlotOptList[2],
    type=float, metavar='float', default=None,
    help='use-specified plotting height to replace the program-determined one, default None')
parser_waterfallPlot.add_argument('-y', '--'+waterfallPlotOptList[3],
    type=str, metavar='{empty,id,index}', default='empty',
    choices=['empty','id','index'],
    help='how should y-axis be labeled, default empty')
parser_waterfallPlot.add_argument('-s', '--'+waterfallPlotOptList[4],
    action='store_true',
    help='should the plot alleles be sorted by length')


######################### testTwoPanels mode #########################
testTwoPanelsPosList = ['inVCF1', 'inVCF2', 'outFile']
testTwoPanelsOptList = ['skipLoci', 'testType', 'varCut', 'outPlotDir']
parserDict['testTwoPanels'] = [testTwoPanelsPosList, testTwoPanelsOptList]
parser_testTwoPanels = subparsers.add_parser('testTwoPanels', description=
'''testTwoPanels:
This command performs statistical tests to compare alleles for each TR locus on two given panel of samples.
Currently 2 test types are supported:
    tf: t-test and f-test of annotation length
    ks: ks-test of motif count distribution
''', formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser_testTwoPanels.add_argument(testTwoPanelsPosList[0], type=str,
    help='string\tinput vcf file 1,  e.g. /in/1.vcf')
parser_testTwoPanels.add_argument(testTwoPanelsPosList[1], type=str,
    help='string\tinput vcf file 2,  e.g. /in/2.vcf')
parser_testTwoPanels.add_argument(testTwoPanelsPosList[2], type=str,
    help='string\toutput tsv file,  e.g. /out/File.tsv')
# optional arguments
parser_testTwoPanels.add_argument('-s', '--'+testTwoPanelsOptList[0],
    type=str, metavar='string', default=None,
    help='list of loci to skip (in bed format),  default None')
parser_testTwoPanels.add_argument('-t', '--'+testTwoPanelsOptList[1],
    type=str, metavar='{tf,ks}', default='tf', choices=['tf','ks'],
    help='type of statistical test,  default tf')
parser_testTwoPanels.add_argument('-c', '--'+testTwoPanelsOptList[2],
    type=float, metavar='float', default=1,
    help='variance lower bound for tested (tf only),  default 1')
parser_testTwoPanels.add_argument('-o', '--'+testTwoPanelsOptList[3],
    type=str, metavar='string', default=None,
    help='dir of histograms of sig loci (tf only),  default None')


######################### pairwiseCompare mode #########################
pairwiseComparePosList = ['inVCF1', 'inVCF2', 'outFile']
pairwiseCompareOptList = ['skipLoci']
parserDict['pairwiseCompare'] = [pairwiseComparePosList, pairwiseCompareOptList]
parser_pairwiseCompare = subparsers.add_parser('pairwiseCompare', description=
'''pairwiseCompare:
This command generates allele comparison for TR loci of given pair of single haplotype vamos vcfs.
''', formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser_pairwiseCompare.add_argument(pairwiseComparePosList[0], type=str,
    help='string\tinput vcf file 1,  e.g. /in/1.vcf')
parser_pairwiseCompare.add_argument(pairwiseComparePosList[1], type=str,
    help='string\tinput vcf file 2,  e.g. /in/2.vcf')
parser_pairwiseCompare.add_argument(pairwiseComparePosList[2], type=str,
    help='string\toutput txt file,  e.g. /out/File.txt')
# optional arguments
parser_pairwiseCompare.add_argument('-s', '--'+pairwiseCompareOptList[0],
    type=str, metavar='string', default=None,
    help='list of loci to skip (in bed format),  default None')


########## strict positional/optional arguments checking ##########
argsDict = vars(parser.parse_args())
if argsDict['command'] is None:
    parser.print_help()
    sys.exit(0)
    
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
        combineVCF.combineVcfs(inVCFs, outVCF)

    ######################### quickFeature mode #########################
    if command == 'quickFeature':

        import lib.quickFeature as quickFeature

        logging.info(f'Reading input vcf to generate feature matrix...')
        quickFeature.writeFeature(inVCF, feature, byDip, outFile, \
                                  skipLoci, demographics)

    ######################### waterfallPlot mode #########################
    elif command == 'waterfallPlot':
        import lib.waterfallPlot as waterfallPlot

        logging.info(f'Generating waterfall plot...')
        waterfallPlot.plot(inVCF, outDir, useLoci, width, height, ylabel, sort)


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

