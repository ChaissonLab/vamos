import sys
from scipy.stats import entropy
import math
from collections import Counter
import edlib
import logging

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

'''def shannon_entropy(sequence):
    """
    Calculates the Shannon entropy of a DNA sequence.

    Args:
        sequence (str): The DNA sequence (e.g., "ATGC").

    Returns:
        float: The Shannon entropy of the sequence.
    """
    base_counts = Counter(sequence)
    probabilities = [count / len(sequence) for count in base_counts.values()]
    print(probabilities)
    return entropy(probabilities)
'''

def purity(motifs, counts, consensus):
    purityList = []
    counts = [ int(c) for c in counts ]
    for i,m in enumerate(motifs):
        if m == consensus:
            ed = 0
        else:
            ed = edlib.align(m, consensus, mode='NW')['editDistance']
        if ed > len(consensus):
            purity = 0
        else:
            purity = 1 - ed / len(consensus)
        purityList.append(purity * counts[i])

    purity = sum(purityList) / sum(counts)

    return purity

inFile = sys.argv[1] # input processed catalog file
outFile = inFile.replace('.tsv','_entropy_purity.tsv')
inCSV = sys.argv[2]
#inOriMotifFile = sys.argv[2] # input original motif file
#inEffMotifFile = sys.argv[3] # input efficient motif file

with open(inCSV) as f:
    for i,line in enumerate(f):
        if i == 0: inOriMotifFile = line.strip().split(',')[1]
        if i == 1: inEffMotifFile = line.strip().split(',')[1]

oriDict, effDict, consDict = {}, {}, {}
with open(inOriMotifFile) as f:
    for i,line in enumerate(f):
        if i == 0: continue
        chr,start,end,motifs,counts,consensus = line.strip().split()[:6]
        #if start != '203810128': continue
        if i % 1000000 == 0:
            logging.info(f'oriMotif: {i} loci handled')
        consDict[(chr,start,end)] = consensus
        motifs = motifs.split(',')
        counts = counts.split(',')
        oriDict[(chr,start,end)] = purity(motifs, counts, consensus)

with open(inEffMotifFile) as f:
    for i,line in enumerate(f):
        if i == 0: continue
        chr,start,end,motifs,counts = line.strip().split()[:5]
        consensus = consDict[(chr,start,end)]
        #if start != '203810128': continue
        if i % 1000000 == 0:
            logging.info(f'effMotif: {i} loci handled')
        motifs = motifs.split(',')
        counts = counts.split(',')
        effDict[(chr,start,end)] = purity(motifs, counts, consensus)

out = open(outFile, 'w')
with open(inFile) as f:
    for line in f:
        chr,start,end,motifs,_,ty,period,cons = line.strip().split()[:8]
        probs = [ cons.count(b) for b in ['A','C','G','T'] ]
        #entr = shannon_entropy(cons)
        entr = entropy(probs)
        temp = f'\t{str(entr)}\t{oriDict[(chr,start,end)]}\t{effDict[(chr,start,end)]}'
        out.write(line.strip()+temp+'\n')
out.close()




