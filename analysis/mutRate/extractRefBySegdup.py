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
inSegdup = 'grch38_superdups.max_identity.bed'
outFasta = 'extractRefSegdup.fa'

# read in ref fasta
logging.info(f'reading reference fasta...')
seqs = { str(fa.id):str(fa.seq) for fa in SeqIO.parse(open(ref),'fasta') }

out = open(outFasta, 'w')
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
out.close()

