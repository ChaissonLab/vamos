out = open('liftoverUCSC.input.bed', 'w')
with open('Gymrek_etal_SupplementalData1_v2.bed') as f:
    for i,line in enumerate(f):
        if i == 0: continue
        chr,start,end = line.strip().split()[:3]
        out.write(f'chr{chr}\t{start}\t{end}\n')
out.close()
