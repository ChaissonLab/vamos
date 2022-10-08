# -*- coding: UTF-8 -*-

import os, sys, getopt, datetime

def parse(argv):

    try:
        opts, args = getopt.getopt(argv, 'hi:o:', ['inVCF=', 'outVCF='])
    except getopt.GetoptError as err:
        print(err)
        sys.stdout.write('\ntest.py -i <inVCFs> -o <outVCF>\n'); sys.stdout.flush()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            sys.stdout.write('\ntest.py -i <inVCFs> -o <outVCF>\n'); sys.stdout.flush()
            sys.exit()

    if len(opts) != 2:
        sys.stdout.write('\ntest.py -i <inVCFs> -o <outVCF>\n'); sys.stdout.flush()
        sys.exit()

    sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - Parsing Input Arguements...\n\n'); sys.stdout.flush()
    for opt, arg in opts:
        if opt in ('-i', '--inVCFs'):
            inVCFs = arg
            sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Required Argument - ' + 'inVCFs' +': '+ str(inVCFs) + '\n'); sys.stdout.flush()
        elif opt in ('-o', '--outVCF'):
            outVCF = arg
            sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Required Argument - ' + 'outVCF' +': '+ str(outVCF) + '\n'); sys.stdout.flush()

    return(inVCFs, outVCF)


def readVCF(inVCFs):

    with open(inVCFs) as f: vcfList = [line.strip() for line in f]
    with open(vcfList[0]) as f: header = [line for line in f if line.startswith('#')]
    header = [ h for h in header[:len(header)-1] if 'INFO=<ID=' not in h ]
    header.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
    header.append('##INFO=<ID=RU,Number=1,Type=String,Description="Comma separated motif sequences list in the reference orientation">\n')
    header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    header.append('##INFO=<ID=ALTANNO,Number=A,Type=String,Description=""Motif representation for all alleles>\n')

    vcfDict = {}
    samplesAll = []
    for vcf in vcfList:
        with open(vcf) as f:
            for line in f:

                if line.startswith('#CHROM'):
                    sample = line.strip().split('\t')[-1]
                    samplesAll.append(sample)
                if line.startswith('#'): continue

                chr, pos, id, ref, alt, qual, filter, info, format, genotype = line.strip().split('\t')
                if 'ALTANNO_H2' in info:
                    end, ru, svtype, altanno_h1, len_h1, altanno_h2, len_h2 = info.split(';')[:-1]
                else:
                    end, ru, svtype, altanno_h1, len_h1, altanno_h2, len_h2 = info.split(';')[:-1] + ['ALTANNO_H2=', 'LEN=0']
                coor = '__'.join([chr,pos,end])

                if not coor in vcfDict: vcfDict[coor] = {'constant':[chr, pos, id, ref, alt, qual, filter, end, ru], 'sample':[], 'info':[], 'altanno_h1':[], 'altanno_h2':[], 'genotype':[], 'alleles':[]}
                vcfDict[coor]['sample'].append(samplesAll[-1])
                vcfDict[coor]['info'].append(info)
                vcfDict[coor]['altanno_h1'].append(altanno_h1.split('=')[1])
                vcfDict[coor]['altanno_h2'].append(altanno_h2.split('=')[1])
                vcfDict[coor]['genotype'].append('NULL')

    header.append( '#' + '\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samplesAll) + '\n' )

    return(samplesAll, header, vcfDict)


def allele(vcfDict):

    for coor,coorDict in vcfDict.items():

        alleles = coorDict['altanno_h1'] + [a for a in coorDict['altanno_h2'] if a != '']
        coorDict['alleles'] = list(set(alleles))
        for i,sample in enumerate(coorDict['sample']):
            if coorDict['altanno_h2'][i] == '':
                coorDict['genotype'][i] = str(coorDict['alleles'].index(coorDict['altanno_h1'][i])+1) +'/'+ str(coorDict['alleles'].index(coorDict['altanno_h1'][i])+1)
                #coorDict['genotype'][i] = str(coorDict['alleles'].index(coorDict['altanno_h1'][i])) +'/'+ str(coorDict['alleles'].index(coorDict['altanno_h1'][i]))
            else:
                coorDict['genotype'][i] = str(coorDict['alleles'].index(coorDict['altanno_h1'][i])+1) +'/'+ str(coorDict['alleles'].index(coorDict['altanno_h2'][i])+1)
                #coorDict['genotype'][i] = str(coorDict['alleles'].index(coorDict['altanno_h1'][i])) +'/'+ str(coorDict['alleles'].index(coorDict['altanno_h2'][i]))


if __name__ == "__main__":

    inVCFs, outVCF = parse(sys.argv[1:])

    samplesAll, header, vcfDict = readVCF(inVCFs)

    allele(vcfDict)

    out = open(outVCF, 'w')
    for h in header: out.write(h)
    for coor,coorDict in vcfDict.items():

        alleles = ','.join([ a.replace(',','-') for i,a in enumerate(coorDict['alleles']) ])
        #alleles = ';'.join([ 'ALTANNO_H%s=%s' %(i,a) for i,a in enumerate(coorDict['alleles']) ])
        info = '%s;%s;SVTYPE=VNTR;ALTANNO=%s' %(coorDict['constant'][7], coorDict['constant'][8], alleles)
        genotypes = [ coorDict['genotype'][coorDict['sample'].index(sample)] if sample in coorDict['sample'] else '.' for sample in samplesAll ]

        out.write('\t'.join(coorDict['constant'][:7] + [info, 'GT'] + genotypes) +'\n')

    out.close()

