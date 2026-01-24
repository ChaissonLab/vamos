import os

inFile = 'all.mutRate.filterHotspot.tagBlock.closestRecomb.tsv'

def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)

boundaryDict = {}
with open(inFile) as f:
    for line in f:
        chr,start,end,numPattern,aveBlock,chrR,startR,endR = line.strip().split()[:8]
        if chr not in boundaryDict: boundaryDict[chr] = []
        boundaryDict[chr] += [startR,endR]
for chr,v in boundaryDict.items():
    boundaryDict[chr] = [int(v[0]), int(v[-1])]

out = open('temp', 'w')
out.write('\t'.join(['chr','start','end','numPattern','aveBlock','chrR','startR','endR','dist','telomere']) +'\n')
with open(inFile) as f:
    for line in f:
        chr,start,end,numPattern,aveBlock,chrR,startR,endR = line.strip().split()[:8]
        start,end,startR,endR = int(start),int(end),int(startR),int(endR)
        # get distance
        if overlap(start,end,startR,endR) == 0:
            if end <= startR:
                dist = startR - end
            else:
                dist = start - endR
        else:
            dist = 0
        # tag telomere TRs
        if end <= boundaryDict[chr][0] or start >= boundaryDict[chr][1]:
            telomere = 'telo-yes'
        else:
            telomere = 'telo-no'
        out.write(line.strip()+f'\t{dist}\t{telomere}\n')
out.close()
os.system(f'mv temp {inFile}')
