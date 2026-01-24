import os
import bedlib

# for oriMotif
vc = 'vamosVariantCluster_v3.0_oriMotifs.tsv'
genomic = 'vamosGenomicTR_v3.0_oriMotifs.tsv'
outFile = 'vamosExpanded_v3.0_oriMotifs.tsv'
# for effMotif-0.1
#vc = 'vamosVariantCluster_v3.0_effMotifs-0.1.tsv'
#genomic = 'vamosGenomicTR_v3.0_effMotifs-0.1.tsv'
#outFile = 'vamosExpanded_v3.0_effMotifs-0.1.tsv'
# for effMotif-0.1 entropy & purity
vc = 'vamosVariantCluster_v3.0_effMotifs-0.1_entropy_purity.tsv'
genomic = 'vamosGenomicTR_v3.0_effMotifs-0.1_entropy_purity.tsv'
outFile = 'vamosExpanded_v3.0_effMotifs-0.1_entropy_purity.tsv'


vcDict = {}
with open(vc) as f:
    for line in f:
        chrom,start,end = line.strip().split()[:3]
        vcDict[(chrom,start,end)] = [line]

genomicDict = {}
with open(genomic) as f:
    for line in f:
        chrom,start,end = line.strip().split()[:3]
        genomicDict[(chrom,start,end)] = [line,f'{chrom}_{start}_{end}']

tempDict = { k:coor for k,(_,coor) in genomicDict.items() }
vcDict = bedlib.bedIntersectTag(vcDict, tempDict)

out = open(outFile, 'w')
for (chrom,start,end),(line,coor) in genomicDict.items():
    out.write(line)
for (chrom,start,end),(line,tags) in vcDict.items():
    if tags == 'NULL':
        out.write(line)
out.close()

# sort
os.system(f'bedtools sort -g /project2/mchaisso_100/bidagu/databases/references/name.txt -i {outFile} > tmp')
os.system(f'mv tmp {outFile}')
