# -*- coding: UTF-8 -*-
import os
import scipy.stats as stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import datetime
import logging

from lib.tr import TR
import lib.general as general


logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)


def twoSampleTF(data1:list[int], data2:list[int], equalVar:bool=True) \
    -> tuple[int,int]:
    """function to perform t-test & f-test for one locus on two sample panels

    Args:
        data1 (list[int]): input tr data of panel 1 samples
        data2 (list[int]): input tr data of panel 2 samples
        equalVar (bool, optional): equal or unequal variance for two sample
            t-test. Defaults to True.

    Returns:
        tuple[int,int]: returned t test and f test p-value
    """

    t_stat,tp = stats.ttest_ind(a=data1, b=data2, equal_var=equalVar)

    F = np.var(data1) / np.var(data2)
    df1, df2 = len(data1)-1, len(data2)-1
    fp = stats.f.cdf(F, df1, df2)

    return tp, fp


def twoSampleKS(data1:list[int], data2:list[int]) -> int:
    """function to perform KS test for one locus on two sample panels

    Args:
        data1 (list[int]): input tr data of panel 1 samples
        data2 (list[int]): input tr data of panel 2 samples

    Returns:
        int: returned KS test p-value
    """

    ks_stat,ksp = stats.kstest(data1, data2)

    return ksp


def test(inVCF1, inVCF2, testType, skipLoci, varCut, outFile, outPlotDir):

    if skipLoci:
        logging.info('Reading input loci for skipping...')
        skipDict = general.parseBed(skipLoci)
    else:
        skipDict = {}

    if outPlotDir:
        if not os.path.exists(outPlotDir): os.makedirs(outPlotDir)

    chrList = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9', \
            'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17', \
            'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

    file1 = open(inVCF1, 'r')
    file2 = open(inVCF2, 'r')
    line1 = file1.readline()
    line2 = file2.readline()

    out = open(outFile, 'w')

    counter = 1
    while line1 and line2:
        # handle headers
        if line1.startswith('##'):
            line1 = file1.readline()
            continue
        if line2.startswith('##'):
            line2 = file2.readline()
            continue
        if line1.startswith('#'):
            samples1 = line1.strip().split()[9:]
            line1 = file1.readline()
            continue
        if line2.startswith('#'):
            samples2 = line2.strip().split()[9:]
            line2 = file2.readline()
            continue

        if counter % 100000 == 0: logging.info(f'{counter} loci tested.')
        tr1, tr2 = TR(), TR()
        tr1.parseDiploidVCFOneLine(samples1, line1)
        tr2.parseDiploidVCFOneLine(samples2, line2)

        # different chr
        if tr1.chr != tr2.chr:
            # line1 go past line2, increase line2
            if chrList.index(tr1.chr) > chrList.index(tr2.chr):
                line2 = file2.readline()
            # line2 go past line1, increase line2
            else:
                line1 = file1.readline()

        # same chr but unmatched coordinate
        elif tr1.start != tr2.start:
            # line1 goes beyond, skip this line2
            if int(tr1.start) > int(tr2.start):
                line2 = file2.readline()
            # line2 goes beyond, skip this line1
            else:
                line1 = file1.readline()

        # same chr and same coordinate (proceed to test)
        else:

            if tr1.chr in skipDict:
                if (tr1.start, tr1.end) in skipDict[tr1.chr]:
                    line1 = file1.readline()
                    line2 = file2.readline()
                    continue

            if testType == 'tf':
            # test by annotation length
                anno1 = [ a for s,a in tr1.annosByUsed.items() if a != ['.'] ]
                anno2 = [ a for s,a in tr2.annosByUsed.items() if a != ['.'] ]

                data1 = [ len(a) for a in anno1 ]
                data2 = [ len(a) for a in anno2 ]

                if np.var(data1) > varCut and np.var(data2) > varCut:
                    tp,fp = twoSampleTF(data1, data2)
                    # plot histogram contrasting two groups
                    if outPlotDir:
                        lower,upper = min(data1+data2), max(data1+data2)
                        n = round(len(data1+data2)/10)
                        bins = np.linspace(lower,upper,n)
                        outPlot = f'{outPlotDir}/{tr1.chr}_{tr1.start}-{tr1.end}.png'
                        plt.figure(figsize=(15, 12))
                        plt.hist(data1, bins, alpha=0.5, label='panel1')
                        plt.hist(data2, bins, alpha=0.5, label='panel2')
                        plt.legend(loc='upper right')
                        plt.savefig(outPlot, dpi=300, format='png')
                    data1 = ','.join(map(str,data1))
                    data2 = ','.join(map(str,data2))
                    sig = tp < 0.00000005
                    temp = [tr1.chr,tr1.start,tr1.end,sig,tp,fp,data1,data2]
                    out.write('\t'.join(map(str,temp)) +'\n')

            # test by motif count distribution (ks-test)
            elif testType == 'ks':
                anno1 = [ a for s,a in tr1.annosByUsed.items() if a != ['.'] ]
                anno2 = [ a for s,a in tr2.annosByUsed.items() if a != ['.'] ]

                numMotifs = len(tr1.motifsUsed)
                data1, data2 = [], []
                for i in range(numMotifs):
                    data1.append( sum([ a.count(str(i)) for a in anno1 ]) )
                for i in range(numMotifs):
                    data2.append( sum([ a.count(str(i)) for a in anno2 ]) )

                data1 = [ round(d/len(anno1),2) for d in data1 ]
                data2 = [ round(d/len(anno2),2) for d in data2 ]

                if len(anno1) > 5 and len(anno2) > 5:
                    #if np.var(data1) != 0 and np.var(data2) != 0:
                    ksp = twoSampleKS(data1, data2)
                    data1 = ','.join(map(str,data1))
                    data2 = ','.join(map(str,data2))
                    sig = ksp < 0.00000005
                    temp = [tr1.chr,tr1.start,tr1.end,sig,ksp,data1,data2]
                    out.write('\t'.join(map(str,temp)) +'\n')

            line1 = file1.readline()
            line2 = file2.readline()

        counter += 1
    out.close()

