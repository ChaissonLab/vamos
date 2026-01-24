import os
import sys
import logging

import edlib
import seqAlign 

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

inBlock = sys.argv[1]
inMutRate = sys.argv[2]
outFile = sys.argv[3]

def get_all_rotations(s1):
    rotations = []
    for i in range(len(s1)):
        rotations.append(s1[i:] + s1[:i])
    return rotations

def group_similar_strings(strings, threshold = 0.4):
    groups = {}
    processed = {}
    maxNumMotifs = 0
    for s in strings:
        for m in s:
            if int(m) > maxNumMotifs:
                maxNumMotifs = int(m)
    maxNumMotifs += 1
    logging.info(f'max number of motifs: {maxNumMotifs}')

    for i,s1 in enumerate(strings,start=1):
        if i % 5000 == 0:
            logging.info(f'{i} out of {len(strings)} patterns checked')
        if s1 in processed:
            continue

        group = [s1]
        processed[s1] = ''

        # get all cyclic rotations of s1
        rotations = get_all_rotations(s1)
        
        for s2 in strings:
            if s2 not in processed:
                # Check against all rotations of s1
                for rot in rotations:
                    if rot == s2:
                        ed = 0
                    elif abs(len(s1)-len(s2)) / len(s1) > threshold:
                        ed = 2*len(s1)
                    elif maxNumMotifs <= 256:
                        anno1 = [chr(int(a)) for a in rot]
                        anno2 = [chr(int(a)) for a in s2]
                        ed = edlib.align(anno1, anno2)['editDistance']
                    else:
                        ed = seqAlign.alignGlobal(rot, s2, distance=True, match=0, mismatch=1, indel=1)
                    if ed / len(s1) <= threshold:
                        group.append(s2)
                        processed[s2] = ''
                        break

        groups[s1] = group

    return dict(groups)

mutRateDict = {}
with open(inMutRate) as f:
    for i,line in enumerate(f):
        if i == 0: continue
        chrom,start,end,treeS,treeE,TR2Tree,nSNP,treeLen,winS,winE,period,type,annotated,totalAllele,spanned,totalEd,\
            totalDist,totalLen,rate,rateNormBefore,rateNormAfter = line.strip().split('\t')
        if int(period) < 7: continue
        if spanned != 'full': continue
        mutRateDict[(chrom,start,end,winS,winE)] = treeS,treeE,period,type,annotated,totalAllele,spanned,totalEd,totalDist,\
            totalLen,rate,rateNormBefore,rateNormAfter

out = open(outFile, 'w')

_,CHROM,winS,winE = outFile.split('_')[:4]
logging.info(f'handling window: {CHROM}:{winS}-{winE}')
with open(inBlock) as f:
    for i,line in enumerate(f):

        chrom,start,end,annotated,block,patterns,counts,nAffected = line.strip().split('\t')[:8]
        alleles = line.strip().split('\t')[8:]

        if block == 'False': continue
        if not (chrom,start,end,winS,winE) in mutRateDict: continue
        treeS,treeE,period,type,annotated,totalAllele,spanned,totalEd,totalDist,totalLen,rate,rateNormBefore,\
            rateNormAfter = mutRateDict[(chrom,start,end,winS,winE)]

        logging.info(f'handling locus {chrom}:{start}-{end}')
        # group patterns by similarity (cyclic rotation also considered)
        patterns = [ tuple(p.split('-')) for p in patterns.split(',') ]
        logging.info(f'total number of patterns to check: {len(patterns)}')
        logging.info(f'total number of unique patterns to check: {len(set(patterns))}')
        grouppedPatterns = group_similar_strings(patterns, threshold = 0.4)


        for i,group in enumerate(grouppedPatterns.keys()):
            patternsOfThisGroup = { '-'.join(p):i for p in grouppedPatterns[group] }
            alleleByBlocks = []
            for allele in alleles:
                blocks = 0
                for a in allele.split(','):
                    if a in patternsOfThisGroup: blocks += 1
                alleleByBlocks.append(blocks)

            cons = '-'.join(group)
            patternsOfThisGroup = ','.join(list(patternsOfThisGroup.keys()))
        
            out.write('\t'.join(map(str,[chrom,start,end,treeS,treeE,period,type,annotated,spanned,cons,\
                                         patternsOfThisGroup]+alleleByBlocks)) +'\n')

out.close()

