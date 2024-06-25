# -*- coding: UTF-8 -*-
import sys
import datetime
import logging

import lib.general as general

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)


def getVCFSampleName(vcf):
    with open(vcf) as f:
        for line in f:
            if line.startswith('#CHROM'):
                sample = line.strip().split('\t')[-1]
                break

    return sample


def parseSingleVCFOneLine(line):

    chr,pos,id,ref,alt,qual,filter,info,_,gt = line.strip().split('\t')
    if 'ALTANNO_H2' in info:
        end,ru,_,altH1,_,altH2,_ = info.split(';')[:-1]
    else:
        end,ru,_,altH1,_ = info.split(';')[:-1]
        altH2 = 'ALTANNO_H2=.'
    altH1 = altH1.split('=')[1]
    altH2 = altH2.split('=')[1]
    gt1,gt2 = gt.split('/')
    gt1 = [altH1, altH2][int(gt1)-1]
    gt2 = [altH1, altH2][int(gt2)-1]

    constant = [chr, pos, id, ref, alt, qual, filter, end, ru]

    return(chr, pos, constant, gt1, gt2)


# read one haploid/diploid vcf
def readSingleVCF(vcf,vcfDict):

    with open(vcf) as f:
        for line in f:

            if line.startswith('#CHROM'):
                sample = line.strip().split('\t')[-1]
            if line.startswith('#'): continue

            chr, start, constant, altH1, altH2 = parseSingleVCFOneLine(line)

            if chr not in vcfDict: vcfDict[chr] = {}
            if start not in vcfDict[chr]:
                vcfDict[chr][start] = { 'constant':constant, 'sample':[], \
                    'altH1':[], 'altH2':[], 'genotype':[], 'alleles':[] }

            vcfDict[chr][start]['sample'].append(sample)
            vcfDict[chr][start]['altH1'].append(altH1)
            vcfDict[chr][start]['altH2'].append(altH2)
            vcfDict[chr][start]['genotype'].append('NULL')


# read two haploid vcfs (note that "altH1" is always taken for each vcf)
def readPairedVCF(vcf1,vcf2,vcfDict):

    tempDict = {}
    with open(vcf1) as f:
        for line in f:

            if line.startswith('#CHROM'):
                sample1 = line.strip().split('\t')[-1]
            if line.startswith('#'): continue
            chr, start, constant, altH1, altH2 = parseSingleVCFOneLine(line)
            tempDict[(chr,start)] = [altH1, '.']

            if chr not in vcfDict: vcfDict[chr] = {}
            if start not in vcfDict[chr]:
                vcfDict[chr][start] = { 'constant':constant, 'sample':[], \
                    'altH1':[], 'altH2':[], 'genotype':[], 'alleles':[] }

    with open(vcf2) as f:
        for line in f:

            if line.startswith('#CHROM'):
                sample2 = line.strip().split('\t')[-1]
            if line.startswith('#'): continue
            chr, start, constant, altH1, altH2 = parseSingleVCFOneLine(line)
            if (chr,start) in tempDict:
                tempDict[(chr,start)][1] = altH1
            else:
                tempDict[(chr,start)] = ['.', altH1]

            if chr not in vcfDict: vcfDict[chr] = {}
            if start not in vcfDict[chr]:
                vcfDict[chr][start] = { 'constant':constant, 'sample':[], \
                    'altH1':[], 'altH2':[], 'genotype':[], 'alleles':[] }

    sample = sample1+'/'+sample2
    for (chr,start),(altH1,altH2) in tempDict.items():

        vcfDict[chr][start]['sample'].append(sample)
        vcfDict[chr][start]['altH1'].append(altH1)
        vcfDict[chr][start]['altH2'].append(altH2)
        vcfDict[chr][start]['genotype'].append('NULL')


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

    # make header for the combined vcf
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
            logging.info(f'Reading single vcf: {vcfs[0]}')
            readSingleVCF(vcfs[0], vcfDict)
        elif len(vcfs) == 2:
            logging.info(f'Reading pair of vcfs: {vcfs[0]} & {vcfs[1]}')
            readPairedVCF(vcfs[0], vcfs[1], vcfDict)
        else:
            sys.exit('Error on input! Please check the input instructions.')

    # sort vcfDict by start for each chr
    for chr,chrDict in vcfDict.items():
        vcfDict[chr] = {k:v for k,v in sorted(chrDict.items(), \
                                              key=lambda x: int(x[0]))}

    return(samplesAll, header, vcfDict)


def allele(vcfDict):

    for chr,chrDict in vcfDict.items():
        for start,coorDict in chrDict.items():

            # save all alleles to coorDict['alleles'] for allele indexing
            for a in coorDict['altH1']:
                if a not in coorDict['alleles'] and a != '.':
                    coorDict['alleles'].append(a)
            for a in coorDict['altH2']:
                if a not in coorDict['alleles'] and a != '.':
                    coorDict['alleles'].append(a)

            # allele indexing
            for i,sample in enumerate(coorDict['sample']):

                if not coorDict['altH1'][i] == '.':
                    gt1 = str(coorDict['alleles'].index(coorDict['altH1'][i])+1)
                else:
                    gt1 = '.'

                if not coorDict['altH2'][i] == '.':
                    gt2 = str(coorDict['alleles'].index(coorDict['altH2'][i])+1)
                else:
                    gt2 = '.'

                coorDict['genotype'][i] = gt1 +'/'+ gt2


def combineVCF(inVCFs, outVCF):

    samplesAll, header, vcfDict = readAllVCF(inVCFs)

    logging.info(f'Summarizing alleles for each locus...')
    allele(vcfDict)

    logging.info(f'Writting combined vcf...')
    out = open(outVCF, 'w')
    for h in header: out.write(h)
    for chr,chrDict in vcfDict.items():
        for start,coorDict in chrDict.items():

            alleles = ','.join([ a.replace(',','-') \
                                for i,a in enumerate(coorDict['alleles']) ])

            info = '%s;%s;SVTYPE=VNTR;ALTANNO=%s' \
                %(coorDict['constant'][7], coorDict['constant'][8], alleles)

            # if the sample is present in 
            genotypes = [coorDict['genotype'][coorDict['sample'].index(s)] \
                if s in coorDict['sample'] else './.' for s in samplesAll]

            temp = coorDict['constant'][:7] + [info,'GT'] + genotypes
            out.write('\t'.join(temp) + '\n')

    out.close()

