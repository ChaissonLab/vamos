import os

out = open('genes.tsv', 'w')
written = []
with open('EnsembleID.BioMart.output.txt') as f:
    for i,line in enumerate(f):
        if i == 0: continue
        id,_,_,_,chrom,start,end,_,_ = line.strip().split('\t')
        if chrom == 'MT': chrom = 'M'
        if (chrom,start,end) not in written:
            out.write('\t'.join([f'chr{chrom}',start,end]) +'\n')
            written.append((chrom,start,end))
out.close()

os.system(f'bedtools sort -g /project2/mchaisso_100/bidagu/databases/references/name.txt -i genes.tsv > temp')
os.system(f'mv temp genes.tsv')

os.system(f'bedtools sort -g /project2/mchaisso_100/bidagu/databases/references/name.txt -i selectedATs.tsv > temp')
os.system(f'mv temp selectedATs.tsv')
