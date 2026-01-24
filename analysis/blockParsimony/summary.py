
allDict = {}
with open('all.parsimonyBlock.tsv') as f:
    for line in f:
        chr,start,end,pattern,count = line.strip().split()[:5]
        if (chr,start,end) not in allDict: allDict[(chr,start,end)] = []
        allDict[(chr,start,end)].append([pattern,count])

out = open('all.parsimonyBlock.merged.tsv', 'w')
out.write('\t'.join(['chr','start','end','totalRecurrence','patterns','counts']) +'\n')
for (chr,start,end),v in allDict.items():
    patterns,counts,total = '','',0
    patterns = ','.join([ p for p,c in v ])
    counts = ','.join([ c for p,c in v ])
    for p,c in v:
        total += int(c)
    out.write('\t'.join([chr,start,end,str(total),patterns,counts]) +'\n')
out.close()
