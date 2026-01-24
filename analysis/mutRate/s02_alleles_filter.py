import os
import sys
import logging
import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

inCHRSamples = sys.argv[1] # kept sample for this chromosome
inCHRBed = sys.argv[2] # mask bed for this chromosome
inSNPVCF = sys.argv[3] # input master SNP VCF
inTRVCF = sys.argv[4] # input master TR VCF
coor = sys.argv[5] # coordinate chr_start_end
outCHRSNP = sys.argv[6]
outCHRTR = sys.argv[7]

def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)

CHR,START,END = coor.split('_')

# function to parse input dipcall vcfs
def parseDipcall(inSNPVCF, CHR):

    snpDict = {}
    with open(inSNPVCF) as f:
        for idx,line in enumerate(f,start=1):
            if idx % 5000000 == 0: logging.info(f'{idx} lines processed')
            if line.startswith('#'): continue
            chr,start,_,ref,alt,_,_,_,_ = line.strip().split()[:9]
            if CHR == 'chr23' and chr == 'chrX': chr = 'chr23'
            if CHR == 'chr24' and chr == 'chrY': chr = 'chr24'
            ref, alt = ref.upper(), alt.upper()

            gts = line.strip().split()[9:]
            #gts = gts[:useSamples+1]

            # check for chromosome match
            if chr != CHR: continue
            # only keep in-range SNPs
            if int(start) < int(START) or int(start) > int(END): continue

            # keep only single base bialleleic SNPs
            if len(ref) != 1 or len(alt) != 1: continue
            if ref == 'N' or alt == 'N': continue

            gtsByAncient = []
            for i,gt in enumerate(gts):
                # get the gts by base for this sample
                if '/' in gt: gt1,gt2 = gt.split(':')[0].split('/')
                if '|' in gt: gt1,gt2 = gt.split(':')[0].split('|')
                gt1Base, gt2Base = '.', '.'
                if gt1 != '.': gt1Base = [ref,alt][int(gt1)]
                if gt2 != '.': gt2Base = [ref,alt][int(gt2)]
                # find the ancient allele
                if i == 0:
                    ancient = '.'
                    if gt1Base == gt2Base: ancient = gt1Base
                    mutation = ref if ancient == alt else alt
                    # skip SNP that has ambiguous ancient allele
                    if ancient == '.': break
                else:
                    # decide gts as ancient allele or not
                    gt1ByAncient, gt2ByAncient = '1', '1'
                    if gt1Base == ancient: gt1ByAncient = '0'
                    if gt1Base == '.': gt1ByAncient = '.'
                    if gt2Base == ancient: gt2ByAncient = '0'
                    if gt2Base == '.': gt2ByAncient = '.'
                    gtsByAncient += [gt1ByAncient, gt2ByAncient]
            # skip SNP that has ambiguous ancient allele
            if ancient == '.': continue
            snpDict[(chr,start,ancient,mutation)] = gtsByAncient
            if len(snpDict) % 500000 == 0: logging.info(f'{len(snpDict)} SNPs read')

    return snpDict

# get kept sample IDs for this window
with open(inCHRSamples) as f:
    samplesCHR = [ l.strip().split('\t')[0] for l in f ]

# get all masking bed for this window
maskBed = {}
with open(inCHRBed) as f:
    for line in f:
        chr,start,end = line.strip().split('\t')[:3]
        maskBed[(chr,start,end)] = f'{chr}_{start}_{end}'

# get sample IDs for all sample from the master TR VCF
samplesDip, samplesHap = [], []
with open(inTRVCF) as f:
    for line in f:
        if line.startswith('##'): continue
        if line.startswith('#'):
            samplesDip = line.strip().split()[9:]
            for sample in samplesDip: samplesHap += sample.split('/')
            break

# get the kept index (in hap for SNP, in dip for TR)
logging.info(f'configuring samples...')
keptHaps, keptHapIndex, keptDips, keptDipIndex = [], [], [], []
for i,sample in enumerate(samplesHap):
    if sample in samplesCHR:
        keptHaps.append(sample)
        keptHapIndex.append(i)
for i,sample in enumerate(samplesDip):
    hap1,hap2 = sample.split('/')
    if hap1 in samplesCHR and hap2 in samplesCHR:
        keptDips.append(sample)
        keptDipIndex.append(i)
logging.info(f'{len(keptHaps)} samples remained for {coor}')

##### write output TR #####
logging.info(f'filtering TRs...')
out = open(outCHRTR, 'w')
with open(inTRVCF) as f:
    for line in f:
        if line.startswith('##'):
            out.write(line)
            continue
        if line.startswith('#'):
            fields = line.strip().split('\t')
            out.write('\t'.join(fields[:9]+keptDips) +'\n')
            continue
        fields = line.strip().split('\t')
        chr = fields[0]
        start = fields[1]
        end = fields[7].split(';')[0].split('=')[1]
        if chr != CHR: continue
        # keep only overlap TRs for this window
        op = overlap(int(START),int(END),int(start),int(end))
        if op == 0: continue
        # keep only kept samples
        data = [ d for i,d in enumerate(fields[9:]) if i in keptDipIndex ]
        out.write('\t'.join(map(str, fields[:9]+data)) +'\n')
logging.info(f'filtering TRs done')

##### filter SNP #####
# read master SNP
logging.info(f'reading SNPs...')
snpDict = parseDipcall(inSNPVCF, CHR)
logging.info(f'{len(snpDict)} SNPs read in')

# filter mask bed
logging.info(f'masking SNPs...')
keyDict = { (chr,start,start):[] for chr,start,ancient,mutation in snpDict.keys() }
keyDict = bedlib.bedIntersectRemove(keyDict, maskBed)
logging.info(f'{len(keyDict)} SNPs remained after masking')

# write output SNP
logging.info(f'filtering SNPs...')
out = open(outCHRSNP, 'w')
out.write('\t'.join(['chr','pos','ancient','mutation']+keptHaps) +'\n')
for (chr,start,ancient,mutation),data in snpDict.items():
    if (chr,start,start) not in keyDict: continue
    # keep only kept samples
    temp = [ d for i,d in enumerate(data) if i in keptHapIndex ]
    out.write('\t'.join(map(str, [chr,start,ancient,mutation]+temp)) +'\n')
out.close()
logging.info(f'filtering SNPs done')
