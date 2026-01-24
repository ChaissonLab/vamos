outFiles = ['windows.0001-0500.csv','windows.0501-1000.csv','windows.1001-1505.csv']

for outFile in outFiles:
    inFile = '../s02_alleles/'+outFile
    out = open(outFile, 'w')
    with open(inFile) as f:
        for line in f:
            chr,start,end = line.strip().split(',')[:3]
            coor = f'{chr}_{start}_{end}'
            if chr in ['chrX','chrY']: continue
            path1 = f'../../s03_trees/{chr}/{coor}/{coor}.map.newick'
            path2 = f'../../s02_alleles/{chr}/{coor}/{coor}.tr.vcf'
            out.write(','.join([chr,start,end,path1,path2]) +'\n')
    out.close()
