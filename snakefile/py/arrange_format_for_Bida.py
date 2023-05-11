import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 
import re 
from Bio import SeqIO


parser = argparse.ArgumentParser('Arrange assembly fa into the format Bida wants')
parser.add_argument('-f', '--input_fa', nargs='?', required=True)
parser.add_argument('-b', '--input_bed', nargs='?', required=True)
parser.add_argument('-o', '--output', nargs='?', required=True)

args = parser.parse_args()
input_fa = args.input_fa
input_bed = args.input_bed
output = args.output

def readBed(input_bed):
    contig_info = defaultdict(str)

    with open(input_bed, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.rstrip('\n')
            fields = line.split('\t')
            contig, contig_cpname, strand = fields[0], fields[1], int(fields[2])

            if strand == 0: # forward strand
              contig_cpname = contig_cpname[:-1] + '1'
            else: # reverse strand
              contig_cpname = contig_cpname[:-1] + '-1'
            contig_info[contig] = contig_cpname

    return contig_info

def readFa(input_fa):
    asm_fa = defaultdict(str) 
    fasta_sequences = SeqIO.parse(open(input_fa), 'fasta')
    for fasta in fasta_sequences:
        contig, sequence = fasta.id, str(fasta.seq)
        contig_cpname = contig_info[contig]
        asm_fa[contig_cpname] = sequence

    return asm_fa

def outputFa(asm_fa, output):

    with open(output, 'w') as outfile:
        for contig_cpname, seq in asm_fa.items():
            outfile.write(f'>{contig_cpname}' + '\n') 
            outfile.write(seq + '\n')       
    return 


if __name__ == "__main__":
    contig_info = readBed(input_bed) 
    asm_fa = readFa(input_fa)
    outputFa(asm_fa, output)

















