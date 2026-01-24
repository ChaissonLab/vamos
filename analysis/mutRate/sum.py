total = 0

with open('GSE180714_HARs.bed') as f:
    for i,line in enumerate(f):
        if i == 0: continue
        chrom,start,end = line.strip().split('\t')[:3]
        total += int(end) - int(start)
print(total)
