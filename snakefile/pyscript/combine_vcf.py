# -*- coding: UTF-8 -*-

import os
import sys
import getopt
import datetime
import logging

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

def parse(argv):

    try:
        opts, args = getopt.getopt(argv, 'hi:o:', ['inVCFs=', 'outVCF='])
    except getopt.GetoptError as err:
        print(err)
        print('\ntest.py -i <inVCFs.csv> -o <outVCF.vcf>\n')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('\ntest.py -i <inVCFs.csv> -o <outVCF.vcf>\n')
            sys.exit()

    if len(opts) != 2:
        print('\ntest.py -i <inVCFs.csv> -o <outVCF.vcf>\n')
        sys.exit()

    logging.info('Parsing Input Arguements...')
    for opt, arg in opts:
        if opt in ('-i', '--inVCFs'):
            inVCFs = arg
            logging.info('Required Argument - inVCFs: ' + str(inVCFs))
        elif opt in ('-o', '--outVCF'):
            outVCF = arg
            logging.info('Required Argument - outVCF: ' + str(outVCF))
    logging.info('Parsing Input Arguements Completed\n')

    return(inVCFs, outVCF)


def getVCFSampleName(vcf):
    with open(vcf) as f:
        for line in f:
            if line.startswith('#CHROM'):
                sample = line.strip().split('\t')[-1]
                break

    return sample


def readOneVCFLine(line,vcfDict):

    chr,pos,id,ref,alt,qual,filter,info,_,_ = line.strip().split('\t')
    if 'ALTANNO_H2' in info:
        end,ru,_,altannoH1,_,altannoH2,_ = info.split(';')[:-1]
    else:
        end,ru,_,altannoH1,_ = info.split(';')[:-1]
        altannoH2 = 'ALTANNO_H2=.'
    altannoH1 = altannoH1.split('=')[1]
    altannoH2 = altannoH2.split('=')[1]

    coor = '__'.join([chr,pos,end])
    constant = [chr, pos, id, ref, alt, qual, filter, end, ru]

    if not coor in vcfDict:
        vcfDict[coor] = { 'constant':constant, 'sample':[], \
        'altannoH1':[], 'altannoH2':[], 'genotype':[], 'alleles':[] }

    return(coor, constant, altannoH1, altannoH2)


# read one haploid/diploid vcf
def readOneVCF(vcf,vcfDict):

    with open(vcf) as f:
        for line in f:

            if line.startswith('#CHROM'):
                sample = line.strip().split('\t')[-1]
            if line.startswith('#'): continue

            coor, constant, altannoH1, altannoH2 = readOneVCFLine(line,vcfDict)

            vcfDict[coor]['sample'].append(sample)
            vcfDict[coor]['altannoH1'].append(altannoH1)
            vcfDict[coor]['altannoH2'].append(altannoH2)
            vcfDict[coor]['genotype'].append('NULL')


# read two haploid vcfs (note that "altannoH1" is always taken)
def readPairedVCF(vcf1,vcf2,vcfDict):

    tempDict = {}
    with open(vcf1) as f:
        for line in f:

            if line.startswith('#CHROM'):
                sample1 = line.strip().split('\t')[-1]
            if line.startswith('#'): continue
            coor, constant, altannoH1, altannoH2 = readOneVCFLine(line,vcfDict)
            tempDict[coor] = [altannoH1, '.']

    with open(vcf2) as f:
        for line in f:

            if line.startswith('#CHROM'):
                sample2 = line.strip().split('\t')[-1]
            if line.startswith('#'): continue
            coor, constant, altannoH1, altannoH2 = readOneVCFLine(line,vcfDict)
            if coor in tempDict:
                tempDict[coor][1] = altannoH1
            else:
                tempDict[coor] = ['.', altannoH1]

    sample = sample1+'/'+sample2
    for coor,(altannoH1,altannoH2) in tempDict.items():

        vcfDict[coor]['sample'].append(sample)
        vcfDict[coor]['altannoH1'].append(altannoH1)
        vcfDict[coor]['altannoH2'].append(altannoH2)
        vcfDict[coor]['genotype'].append('NULL')


def readAllVCF(inVCFs):

    with open(inVCFs) as f:
        vcfList = [ line.strip() for line in f ]
    vcfList = [ vcf.split(',') for vcf in vcfList ]

    # get all sample names
    samplesAll = []
    for vcfs in vcfList:
        if len(vcfs) == 1:
            samplesAll.append(getVCFSampleName(vcfs[0]))
        else:
            sample = getVCFSampleName(vcfs[0]) +'/'+ getVCFSampleName(vcfs[1])
            samplesAll.append(sample)

    # make the header for the combined vcf
    with open(vcfList[0][0]) as f:
        header = [line for line in f if line.startswith('#')]
    header = [ h for h in header[:len(header)-1] if 'INFO=<ID=' not in h ]
    header.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
    header.append('##INFO=<ID=RU,Number=1,Type=String,Description="Comma separated motif sequences list in the reference orientation">\n')
    header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    header.append('##INFO=<ID=ALTANNO,Number=A,Type=String,Description="Motif representation for all alleles">\n')
    header.append('#' + '\t'.join(['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samplesAll) +'\n')

    vcfDict = {}
    for vcfs in vcfList:
        if len(vcfs) == 1:
            readOneVCF(vcfs[0],vcfDict)
        else:
            readPairedVCF(vcfs[0],vcfs[1],vcfDict)

    return(samplesAll, header, vcfDict)


def allele(vcfDict):

    for coor,coorDict in vcfDict.items():

        for a in coorDict['altannoH1']:
            if a not in coorDict['alleles'] and a != '.':
                coorDict['alleles'].append(a)
        for a in coorDict['altannoH2']:
            if a not in coorDict['alleles'] and a != '.':
                coorDict['alleles'].append(a)

        for i,sample in enumerate(coorDict['sample']):

            if not coorDict['altannoH1'][i] == '.':
                gt1 = str(coorDict['alleles'].index(coorDict['altannoH1'][i])+1)
            else:
                gt1 = '.'

            if not coorDict['altannoH2'][i] == '.':
                gt2 = str(coorDict['alleles'].index(coorDict['altannoH2'][i])+1)
            else:
                gt2 = '.'

            coorDict['genotype'][i] = gt1 +'/'+ gt2


if __name__ == "__main__":

    inVCFs, outVCF = parse(sys.argv[1:])

    samplesAll, header, vcfDict = readAllVCF(inVCFs)

    allele(vcfDict)

    out = open(outVCF, 'w')
    for h in header: out.write(h)
    for coor,coorDict in vcfDict.items():

        alleles = ','.join([ a.replace(',','-') \
                            for i,a in enumerate(coorDict['alleles']) ])

        info = '%s;%s;SVTYPE=VNTR;ALTANNO=%s' \
            %(coorDict['constant'][7], coorDict['constant'][8], alleles)

        genotypes = [coorDict['genotype'][coorDict['sample'].index(sample)] \
            if sample in coorDict['sample'] else './.' for sample in samplesAll]

        out.write('\t'.join(coorDict['constant'][:7] + [info,'GT'] + genotypes))
        out.write('\n')

    out.close()

logging.info('End of Program\n')

