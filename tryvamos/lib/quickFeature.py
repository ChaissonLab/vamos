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


def writeFeature(vcf, feature, outFile, skipLoci, demographics):

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
                temp = []
                for s in samples: temp += s.split('/')
                out.write('\t'.join(['#chr','start','end']+temp)+'\n')

                # writing demographic tags for each sample (each tag as one row)
                if demoDict:
                    for j in range(len(demoDict[temp[0]])):
                        tags = [ demoDict[s][j] for s in temp ]
                        out.write('\t'.join(['#chr','start','end']+tags)+'\n')
                continue

            if i % 100000 == 0: logging.info(f'Processing line: {i}')

            tr = TR()
            tr.parseDiploidVCFOneLine(samples, line)

            # skip unwanted loci
            if tr.chr in skipDict:
                if (tr.start,tr.end) in skipDict[tr.chr]: continue

            alleles = []
            for s,a in tr.annosByUsed.items():
                if a == ['.']:
                    alleles.append('NA')
                else:
                    if feature == 'annoLen':
                        alleles.append( str(len(a)) )
                    elif feature == 'annoStr':
                        alleles.append( '-'.join(a) )
                    elif feature == 'topCount':
                        alleles.append( str(a.count('0')) )
                    else:
                        motifsDict = {v:k for k,v in tr.motifsUsed.items()}
                        alleleNT = '-'.join( [motifsDict[int(i)] for i in a] )
                        alleles.append( alleleNT )

            out.write('\t'.join([tr.chr,tr.start,tr.end]+alleles) +'\n')

    out.close()

