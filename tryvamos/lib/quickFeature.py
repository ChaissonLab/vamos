import datetime
import logging

from lib.tr import TR
import lib.general as general

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)


def parseDemo(demographics):

    demoDict = {}
    with open(demographics) as f:
        for line in f:

            id = line.strip('\t').split()[0]
            demoDict[id] = line.strip('\t').split()[1:]

    return demoDict


def writeFeature(vcf, feature, byDip, outFile, skipLoci, demographics):

    if demographics:
        logging.info('Reading input demographics...')
        demoDict = parseDemo(demographics)
    else:
        demoDict = {}

    if skipLoci:
        logging.info('Reading input loci for skipping...')
        skipDict = general.parseBed(skipLoci)
    else:
        skipDict = {}

    logging.info('Generating features...')
    out = open(outFile, 'w')
    with open(vcf) as f:
        for i,line in enumerate(f,start=1):
            if line.startswith('##'): continue
            if line.startswith('#'):
                samples = line.strip().split()[9:]
                haps = []
                for s in samples: haps += s.split('/')
                if byDip:
                    out.write('\t'.join(['#chr','start','end']+samples)+'\n')
                else:
                    out.write('\t'.join(['#chr','start','end']+haps)+'\n')

                # writing demographic tags for each sample (each tag as one row)
                if demoDict:
                    for j in range(len(demoDict[haps[0]])):
                        if byDip:
                            tags = [ demoDict[s][j] for s in samples ]
                        else:
                            tags = [ demoDict[s][j] for s in haps ]
                        out.write('\t'.join(['#chr','start','end']+tags)+'\n')
                continue

            if i % 100000 == 0: logging.info(f'Processing line: {i}')

            tr = TR()
            tr.parseDiploidVCFOneLine(samples, line)

            # skip unwanted loci
            if tr.chr in skipDict:
                if (tr.start,tr.end) in skipDict[tr.chr]: continue

            alleles = {}
            for s,a in tr.annosByUsed.items():
                if a == ['.']:
                    alleles[s] = 'NA'
                else:
                    if feature == 'annoLen':
                        alleles[s] = str(len(a))
                    elif feature == 'annoLenNT':
                        motifsDict = {v:k for k,v in tr.motifsUsed.items()}
                        alleles[s] = str(sum([ len(motifsDict[int(i)]) for i in a ]))
                    elif feature == 'annoStr':
                        alleles[s] = '-'.join(a)
                    elif feature == 'topCount':
                        alleles[s] = str(a.count('0'))
                    else:
                        motifsDict = {v:k for k,v in tr.motifsUsed.items()}
                        alleles[s] = '-'.join( [motifsDict[int(i)] for i in a] )

            if byDip:
                allelesOut = []
                for s in samples:
                    allele1, allele2 = alleles[s+'_h1'], alleles[s+'_h2']
                    if allele1 == allele2:
                        allelesOut.append(f"{allele1}/-")
                    else:
                        allelesOut.append(f"{allele1}/{allele2}")
            else:
                allelesOut = [ a for s,a in alleles.items() ]

            out.write('\t'.join([tr.chr,tr.start,tr.end]+allelesOut) +'\n')

    out.close()

