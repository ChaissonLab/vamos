import logging

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

maniFiles = ['windows.0001-0500.csv','windows.0501-1000.csv','windows.1001-1505.csv']

def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)


out = open('all.mutRate.tsv','w')
# write header
with open(maniFiles[0]) as f:
    for i,line in enumerate(f):
        if i == 0:
            chrW,startW,endW = line.strip().split(',')[:3]
            coor = f'{chrW}_{startW}_{endW}'
            winFile = f'{chrW}/{coor}/{coor}.tsv'
            with open(winFile) as f2:
                for i,line in enumerate(f2):
                    if i == 0:
                        out.write(line)
                        break
            break

# write data
treeDict, trDict = {}, {}
for manifest in maniFiles:
    with open(manifest) as f:
        for line in f:
            chrW,startW,endW = line.strip().split(',')[:3]
            coor = f'{chrW}_{startW}_{endW}'
            winFile = f'{chrW}/{coor}/{coor}.tsv'
            with open(winFile) as f2:
                for i,line in enumerate(f2):
                    if i == 0:
                        continue
                    chr,start,end,treeS,treeE,TR2Tree,nSNP,treeLen,winS,winE,period,type,\
                        annotated,totalAllele,spanned,totalEd,totalDist,totalLen,rate,rateNormBefore,rateNormAfter = line.strip().split('\t')
                    # skip entries that not picked up by TR2Tree (this entry uses tree that is not the closest tree)
                    if spanned == 'not-unpick':
                        continue
                    # skip constant entries
                    #if totalEd == '0':
                    #    continue
                    if totalEd != '0':
                        # skip trees larger than
                        if int(treeLen) > 10000:
                            logging.info(f'large tree {chr}:{start}-{end}')
                            continue
                        # skip trees with low SNP density
                        if int(nSNP) != 0:
                            if int(treeLen) / int(nSNP) > 500:
                                logging.info(f'low SNP density {chr}:{start}-{end}')
                                continue
                    if winS != '1':
                        winS = int(winS) + 50000
                    winE = int(winE) - 50000
                    # skip TR locus that is not in the window (window was extended by 50kb on both sides)
                    if overlap(int(start),int(end),int(winS),int(winE)) == 0:
                        continue
                    # skip TR locus that is already recorded (occasionally one TR may span two adjacent windows)
                    # this ensures only the first record is taken
                    if (chr,start,end) in trDict:
                        continue
                    if (chr,treeS,treeE) not in treeDict:
                        treeDict[(chr,treeS,treeE)] = []

                    denominator = int(annotated) * (int(annotated)-1) / 2
                    averageEd = int(annotated) * (int(totalEd) / int(totalLen)) / denominator
                    averageDist = int(totalDist.split('.')[0]) / denominator
                    treeDict[(chr,treeS,treeE)] += [(averageEd,averageDist,int(nSNP),int(treeLen))]
                    if (chr,start,end) not in trDict:
                        trDict[(chr,start,end)] = []
                    out.write(line)
out.close()

out = open('all.tree.tsv','w')
for (chr,start,end),averageList in treeDict.items():
    averageEd = sum([average[0] for average in averageList]) / len(averageList)
    averageDist = sum([average[1] for average in averageList]) / len(averageList)
    averageNSNP = sum([average[2] for average in averageList]) / len(averageList)
    averageTreeLen = sum([average[3] for average in averageList]) / len(averageList)
    out.write(f'{chr}\t{start}\t{end}\t{averageNSNP}\t{averageTreeLen}\t{averageEd}\t{averageDist}\n')
out.close()
