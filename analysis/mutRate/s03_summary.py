import sys
import os
import pysam
import numpy as np
import logging

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

windows = sys.argv[1]

with open(windows) as f:
    for line in f:
        CHR,startW,endW,_,samplesFile = line.strip().split(',')
        coorW = f'{CHR}_{startW}_{endW}'
        logging.info(f'handling {coorW}...')
        with open(samplesFile) as g:
            samples = [ l.strip().split('\t')[0] for l in g ]

        # record all sample, all bed data
        bedSampleDict = {}
        for sample in samples:
            sampleFile = f'alignSegdup/{sample}/{sample}.diff.tsv'
            with open(sampleFile) as g:
                for l in g:
                    coor,coorNew,_,tag,_,_,_,lengthRef,lengthContig,diff,prop = l.strip().split('\t')
                    chr,start,end = coor.split('__')

                    # skip segdup entries not in this window
                    if chr != CHR: continue
                    if int(end) < int(startW) or int(start) > int(endW): continue

                    if (chr,start,end) not in bedSampleDict:
                        bedSampleDict[(chr,start,end)] = {}
                    bedSampleDict[(chr,start,end)][sample] = prop

        os.makedirs(f'summary/{CHR}/{coorW}', exist_ok=True)
        out = open(f'summary/{CHR}/{coorW}/{coorW}.summary.tsv', 'w')
        for (chr,start,end),dataDict in bedSampleDict.items():
            outData = []
            numSuspecious = 0
            for sample in samples:
                if sample in dataDict:
                    outData.append(dataDict[sample])
                    if abs(float(dataDict[sample])) > 5:
                        numSuspecious += 1
                else:
                    outData.append('NULL')
            validData = [ abs(float(d)) for d in outData if d != 'NULL' ]
            numMapped = len(validData)
            meanAbs = sum(validData) / len(validData)
            medianAbs = np.median(np.array(validData))

            # output suspecious samples
            if numSuspecious > 0 and numSuspecious < 5000:
                logging.info(f'locus mean diff: {meanAbs}')
                logging.info(f'locus median diff: {medianAbs}')
                for sample in samples:
                    if sample in dataDict:
                        if abs(float(dataDict[sample])) > 5:
                            logging.info(f'suspecious hit: {chr}:{start}-{end}, {sample}, {dataDict[sample]}')

            out.write('\t'.join([CHR,start,end,str(numMapped),str(len(samples)),str(meanAbs),str(medianAbs)]+outData) +'\n')
        out.close()

os.system('cat summary/chr*/*/*summary.tsv > all.segdup.tsv')
os.system('cut -d$"\t" -f1-7 all.segdup.tsv > temp')
os.system('bedtools sort -g /project/mchaisso_100/cmb-16/bidagu/databases/references/name.txt -header -i temp > all.segdup.sort.tsv')
os.system('rm temp')

# read segdups
segdups = 'grch38_superdups.max_identity.bed'
#segdups = '/project/mchaisso_100/shared/references/hg38/regions/segdups/grch38_superdups.bed'
segdupDict = {}
with open(segdups) as f:
    for line in f:
        chr,start,end,identity = line.strip().split()
        if chr not in segdupDict: segdupDict[chr] = []
        segdupDict[chr].append((start,end,identity))

# read recomb hotspot
recombFemale = '/project/mchaisso_100/mchaisso/projects/trcompdb/decode_recomb_hotspot_female.tsv.hg38.bed'
recombMale = '/project/mchaisso_100/mchaisso/projects/trcompdb/decode_recomb_hotspot_male.tsv.hg38.bed'
recombBoth = 'decode_recomb_hotspot_female_male.tsv.hg38.sorted.merge.bed'

recombFemaleDict = {}
with open(recombFemale) as f:
    for line in f:
        chr,start,end,recombRate = line.strip().split()
        if chr not in recombFemaleDict: recombFemaleDict[chr] = []
        recombFemaleDict[chr].append((start,end,recombRate))

recombMaleDict = {}
with open(recombMale) as f:
    for line in f:
        chr,start,end,recombRate = line.strip().split()
        if chr not in recombMaleDict: recombMaleDict[chr] = []
        recombMaleDict[chr].append((start,end,recombRate))

recombBothDict = {}
with open(recombBoth) as f:
    for line in f:
        chr,start,end = line.strip().split()
        if chr not in recombBothDict: recombBothDict[chr] = []
        recombBothDict[chr].append((start,end,'NULL'))

def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)

out = open('all.segdup.sort.insersected.tsv', 'w')
with open('all.segdup.sort.tsv') as f:
    for line in f:
        chr,start,end,_,_,mean,median = line.strip().split()
        hotspotLen = int(end) - int(start)
        # map to segdups
        identities, identityMax = [], 0
        for (s,e,identity) in segdupDict[chr]:
            segdupLen = int(e) - int(s)
            op = overlap(int(start),int(end),int(s),int(e))
            if op > 0:
            #if op > 0.5*hotspotLen or op > 0.5*segdupLen:
                identities.append(identity)
        if identities:
            identityMax = max([float(d) for d in identities])
        else:
            identities = ['NULL']
        # map to recomb hotspots
        opFemale, opMale, opBoth = 0, 0, 0
        for (s,e,recomb) in recombFemaleDict[chr]:
            opFemale += overlap(int(start),int(end),int(s),int(e))
        for (s,e,recomb) in recombMaleDict[chr]:
            opMale += overlap(int(start),int(end),int(s),int(e))
        for (s,e,recomb) in recombBothDict[chr]:
            opBoth += overlap(int(start),int(end),int(s),int(e))
        opFemaleProp = str(opFemale / hotspotLen)
        opMaleProp = str(opMale / hotspotLen)
        opBothProp = str(opBoth / hotspotLen)
        out.write('\t'.join([chr,start,end,mean,median,str(identityMax),','.join(identities),opFemaleProp,opMaleProp,opBothProp]) +'\n')
out.close()

