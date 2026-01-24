import os
import logging
import numpy as np
import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

# mask1: tr & centromere
mask1 = '/project2/mchaisso_100/bidagu/databases/references/GRCh38/masks/vamosExpanded_v3.0_and_centro.concat.sorted.merged.bed'
centromere = '/project2/mchaisso_100/bidagu/databases/references/GRCh38/masks/para_centro_grch38_from_quentin.bed'
genomeBed = '/project2/mchaisso_100/bidagu/databases/references/GRCh38/chrs.1-Y.segCentro.bed'

bamStatsDir = '/project2/mchaisso_100/bidagu/databases/main_collections'
sampleList = 'pick_rmDup_rmKids.list'
coverStats = '/project2/mchaisso_100/bidagu/databases/main_collections/z_bamCoveredRegions_GRCh38_summary_all.tsv'

def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)

# get sample list and haploid pairing info
# the pairing info is used to ensure both haploids of a diploid is filtered if at least one haploid is filtered
# This is because Relate needs both haploid pairs for running
ids, dipDict1, dipDict2, idLast = [], {}, {}, ''
with open(sampleList) as f:
    for i,line in enumerate(f,start=1):
        id = line.strip().split()[0]
        ids.append(id)
        if i % 2 == 0:
            dipDict1[idLast] = id
            dipDict2[id] = idLast
        idLast = id

# seperate the whole genome into windows (skip centromere)
logging.info(f'configuring genomic windows...')
genomeWinDict, chrs = {}, []
with open(genomeBed) as f:
    for line in f:
        chr,start,centroStart,centroEnd,end = line.strip().split()
        if chr not in chrs: chrs.append(chr)
        if start == centroStart:
            chrMainSegs = [(chr,int(centroEnd),int(end))]
        else:
            chrMainSegs = [(chr,int(start),int(centroStart)),(chr,int(centroEnd),int(end))]
        for (chrSeg,startSeg,endSeg) in chrMainSegs:
            windows, startWin, step = [], int(startSeg), 2000000
            for i in range(startSeg,endSeg+1):
                if i % step == 0 and i != startSeg:
                    genomeWinDict[(chr,startWin,i)] = {}
                    startWin = i+1
                elif i == endSeg:
                    genomeWinDict[(chr,startWin,endSeg)] = {}

out = open('genomeWindows.bed', 'w')
for (chrW,startW,endW) in genomeWinDict.keys():
    out.write('\t'.join([chrW,str(startW),str(endW)]) +'\n')
out.close()

genomeWinDict, chrs, flank = {}, [], 50000
with open('genomeWindows.bed') as f:
    for line in f:
        chr,start,end = line.strip().split()
        if chr not in chrs: chrs.append(chr)
        if int(start) - flank > 0:
            start = int(start) - flank
        end = int(end) + flank
        genomeWinDict[(chr,int(start),int(end))] = {}

for chr in chrs:
    count = 0
    for (chrW,startW,endW) in genomeWinDict.keys():
        if chrW == chr: count += 1
    logging.info(f'total number of windows on {chr}: {count}')

# calculate the coverage stats for each window each sample (uncovered bases)
logging.info(f'calculating coverage stats for each genomic window each sample...')
for id in ids:
    batch,sample = id.split('__')
    path = f'{bamStatsDir}/{batch}/bam_mm2_2.23_GRCh38/{sample}/{sample}_bamCoveredRegions.bed'

    sampleBedDict = {}
    with open(path) as f:
        for line in f:
            chr,start,end,cover = line.strip().split()[:4]
            if cover == 'uncovered':
                sampleBedDict[(chr,start,end)] = f'{chr}_{start}_{end}'
    genomeMapDict = {k:[] for k,v in genomeWinDict.items()}
    genomeMapDict = bedlib.bedIntersectTag(genomeMapDict, sampleBedDict)

    for (chrW,startW,endW),tags in genomeMapDict.items():
        genomeWinDict[(chrW,startW,endW)][id] = 0
        tags = tags[-1].split(',')
        if tags == ['NULL']:
            continue
        for tag in tags:
            chrB,startB,endB = tag.split('_')
            op = overlap(int(startW),int(endW),int(startB),int(endB))
            genomeWinDict[(chrW,startW,endW)][id] += op

# filter samples for each window
logging.info(f'filter samples...')
os.system(f'rm -r *chr*')
out = open('windows.csv', 'w')
for (chrW,startW,endW),winDict in genomeWinDict.items():

    logging.info(f'handling window: {chrW}:{startW}-{endW}')

    # obtain cut threshold for this window and the filtered haplotype list
    stats = [ v for k,v in winDict.items() ]
    cut = np.percentile(np.array(stats), 50) # lower 50% coverage cut for this window
    filteredHaps = []
    for i,id in enumerate(ids):
        if winDict[id] > cut:
            if id in dipDict1: filteredHaps += [id, dipDict1[id]]
            if id in dipDict2: filteredHaps += [id, dipDict2[id]]

    for id in ids:
        batch,sample = id.split('__')
        path = f'{bamStatsDir}/{batch}/bam_mm2_2.23_GRCh38/{sample}/{sample}_bamCoveredRegions.bed'
        with open(path) as f:
            for line in f:
                chrS,startS,endS,cover = line.strip().split()[:4]
                if chrS != chrW: continue
                if cover != 'uncovered': continue
                op = overlap(int(startW),int(endW),int(startS),int(endS))
                if op > 0 and int(endS) - int(startS) > 100000:
                    logging.info(f'genome of large missing segment: {id}: {chrS}_{startS}_{endS}')
                    if id in dipDict1: filteredHaps += [id, dipDict1[id]]
                    if id in dipDict2: filteredHaps += [id, dipDict2[id]]

    # output
    key = f'{chrW}_{startW}_{endW}'
    winDir = f'{chrW}/{key}'
    os.system(f'mkdir -p {winDir}')
    outSampleList = winDir +'/'+ sampleList.replace('.list',f'.sample.{key}.tsv')
    outBed = winDir +'/'+ sampleList.replace('.list',f'.mask.{key}.bed')
    sortedBed = winDir +'/'+ sampleList.replace('.list',f'.mask.{key}.sorted.bed')
    mergedBed = winDir +'/'+ sampleList.replace('.list',f'.mask.{key}.merged.bed')
    out1 = open(outSampleList, 'w')
    out2 = open(outBed, 'w')
    out.write(','.join([chrW,str(startW),str(endW),f'../../s01_mask/{mergedBed}',f'../../s01_mask/{outSampleList}']) +'\n')

    # write the tr.centro masking bed first
    with open(mask1) as f:
        for line in f:
            chrT,startT,endT = line.strip().split()[:3]
            if chrT == chrW:
                op = overlap(int(startW),int(endW),int(startT),int(endT))
                if op > 0:
                    temp = str( int(endT) - int(startT) )
                    out2.write('\t'.join([chrT,startT,endT,temp,'TRCentro'])+'\n')

    # filter sample and write uncovered regions in kept sample
    for id in ids:
        batch,sample = id.split('__')
        path = f'{bamStatsDir}/{batch}/bam_mm2_2.23_GRCh38/{sample}/{sample}_bamCoveredRegions.bed'

        # this haploid id is filtered
        if id in filteredHaps: continue

        out1.write('\t'.join([f'{batch}__{sample}',str(cut),str(winDict[id])])+'\n')
        with open(path) as f:
            for line in f:
                chrS,startS,endS,cover = line.strip().split()[:4]
                if chrS != chrW: continue
                if cover != 'uncovered': continue
                op = overlap(int(startW),int(endW),int(startS),int(endS))
                if op > 0:
                    temp = str( int(endT) - int(startT) )
                    if int(temp) >= 10000:
                        out2.write('\t'.join([chrS,startS,endS,temp,'uncover-large'])+'\n')
                    else:
                        out2.write('\t'.join([chrS,startS,endS,temp,'uncover-small'])+'\n')
    out1.close()
    out2.close()

    # sort and merge
    os.system(f'bedtools sort -g /project/mchaisso_100/cmb-16/bidagu/databases/references/name.txt -header -i {outBed} > {sortedBed}')
    os.system(f'bedtools merge -i {sortedBed} > {mergedBed}')

out.close()

