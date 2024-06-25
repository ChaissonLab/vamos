import datetime
import logging

import lib.general as general


logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)


def parseHapVCF(vcf, skipDict):

    outDict = {}

    with open(vcf) as f:
        for line in f:
            if line.startswith('#'): continue

            chr,start,_,_,_,_,_,info,_,gt = line.strip().split()
            end,ru,_,anno,length,_ = info.split(';')
            end = end.split('=')[1]
            anno = anno.split('=')[1].split(',')

            # skip unwanted loci
            if chr in skipDict:
                if (start,end) in skipDict[chr]: continue

            outDict[(chr,start,end)] = anno

    return outDict


def pairwiseCompare(inVCF1, inVCF2, skipLoci, outFile):

    if skipLoci:
        logging.info('Reading input loci for skipping...')
        skipDict = general.parseBed(skipLoci)
    else:
        skipDict = {}

    logging.info('Reading input vcf1...')
    vcfDict1 = parseHapVCF(inVCF1, skipDict)
    logging.info('Reading input vcf2...')
    vcfDict2 = parseHapVCF(inVCF2, skipDict)

    common = [ k for k in vcfDict1.keys() if k in vcfDict2 ]
    logging.info(f'Number of loci in both vcf1 and vcf2: {len(common)}')


    logging.info(f'Calculating edit distance of divergent loci...')
    dists = {}
    for i,key in enumerate(common,start=1):
        if i % 100000 == 0: logging.info(f'{i} loci processed.')
        anno1, anno2 = vcfDict1[key], vcfDict2[key]
        # only output loci with divergent annotation
        if anno1 != anno2:
            dists[key] = general.alignGlobal(anno1, anno2, True, 0, 1, 1)


    logging.info(f'Writting output...')
    out = open(outFile, 'w')

    out.write(f'vcf1: {inVCF1}\n')
    out.write(f'vcf2: {inVCF2}\n')

    out.write(f'Number of loci in both vcf1 and vcf2: {len(common)}\n')
    temp = len(vcfDict1.keys()) - len(common)
    out.write(f'Number of unique loci in vcf1: {temp}\n')
    temp = len(vcfDict2.keys()) - len(common)
    out.write(f'Number of unique loci in vcf2: {temp}\n')

    out.write(f'Number of divergent loci: {len(dists)}\n')
    temp = sum(dists.values())/len(dists.values())
    out.write(f'Mean edit distance of divergent loci: {temp}\n\n')

    # output all distance data
    for (chr,start,end),dist in dists.items():
        out.write('\t'.join([chr,start,end,str(dist)]) +'\n')
    out.close()

