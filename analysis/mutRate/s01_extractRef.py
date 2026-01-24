import sys
import logging
from Bio import SeqIO

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

#ref = sys.argv[1]
#inBed = sys.argv[2]
#outFasta = sys.argv[3]

ref = '/project2/mchaisso_100/bidagu/databases/references/GRCh38/fasta/GCA_000001405.15_GRCh38_no_alt.fasta'
inBed = '../s06_hotspotWindow/all.sorted.hotspot.tsv'
inSegdup = 'grch38_superdups.max_identity.bed'
#inSegdup = 'grch38_superdups.merge.bed'
outFastaHotspot = 'extractRefHotspot.fa'
outFastaSegdup = 'extractRefSegdup.fa'

# read in ref fasta
logging.info(f'reading reference fasta...')
seqs = { str(fa.id):str(fa.seq) for fa in SeqIO.parse(open(ref),'fasta') }

# read in segdup data
out = open(outFastaSegdup, 'w')
segDupdict = {}
with open(inSegdup) as f:
    for line in f:
        chr,start,end,identity = line.strip().split()
        startNew = int(start)
        endNew = int(end)
        startExpand = int(startNew)-50000
        endExpand = int(endNew)+50000
        if startExpand < 0: startExpand = '0'
        startExpand, endExpand = str(startExpand), str(endExpand)
        for id,seq in seqs.items():
            if id == chr:
                extract = seq[int(start)-1:int(end)].upper()
                out.write(f'>{chr}__{start}__{end}--{chr}__{startNew}__{endNew}--{chr}__{startExpand}__{endExpand}\n{extract}\n')
                break

        if chr not in segDupdict: segDupdict[chr] = []
        segDupdict[chr].append((start,end))
out.close()

def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)

# pick the overall coordinates for every TR mutRate hotspot and segdup overlapping regions
# at least the TR mutRate hotspot will be picked so to analyze every TR mutRate hotspot
out = open(outFastaHotspot, 'w')
with open(inBed) as f:
    for line in f:
        chr,start,end,_,_,tag,_ = line.strip().split('\t')
        if tag == 'normal': continue
        opStarts, opEnds = [], []
        for s,e in segDupdict[chr]:
            op = overlap(int(start),int(end),int(s),int(e))
            if op <= 0: continue
            logging.info(f'segdup overlap detected for hotspot {chr}:{start}-{end} with segdup {chr}:{s}-{e}')
            opStarts.append(int(s))
            opEnds.append(int(e))
        # add in also the raw start and end
        opStarts.append(int(start))
        opEnds.append(int(end))
        # pick overall minimum start and maximum end
        # at least the raw start and end will be picked
        startNew = min(opStarts)
        endNew = max(opEnds)
        startExpand = int(startNew)-50000
        endExpand = int(endNew)+50000
        if startExpand < 0: startExpand = '0'
        startExpand, endExpand = str(startExpand), str(endExpand)
        for id,seq in seqs.items():
            if id == chr:
                extract = seq[int(start)-1:int(end)].upper()
                out.write(f'>{chr}__{start}__{end}--{chr}__{startNew}__{endNew}--{chr}__{startExpand}__{endExpand}\n{extract}\n')
                break
out.close()

