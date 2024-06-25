import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 
import re 
from Bio import SeqIO


parser = argparse.ArgumentParser('Arrange VNTR ref file, contains all assemblies fa')
parser.add_argument('-i', '--input_fa_dir', nargs='?', required=True)
parser.add_argument('-b', '--input_bed_dir', nargs='?', required=True)
parser.add_argument('-o', '--output_dir', nargs='?', required=True)

args = parser.parse_args()
input_fa_dir = args.input_fa_dir
input_bed_dir = args.input_bed_dir
output_dir = args.output_dir
asms = [asm.rstrip(".fa") for asm in os.listdir(input_fa_dir) if asm.endswith(".fa")]

def read(asm):
    print(asm)
    contig_region_2_ref_region = defaultdict()
    with open(input_bed_dir + '/' + asm + '.full.bed', 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.rstrip('\n')
            contig, contig_cpname, strand, ref_region = line.split('\t')
            contig_region_2_ref_region[contig] = (ref_region, strand)


    fasta_sequences = SeqIO.parse(open(input_fa_dir + '/' + asm + '.fa'), 'fasta')
    for fasta in fasta_sequences:
        contig_region, sequence = fasta.id, str(fasta.seq)
        ref_region, strand = contig_region_2_ref_region[contig_region]
        ref_to_asm_fa[ref_region].append((contig_region, sequence, ';'.join([contig_region, ref_region, strand])))

    return 

def output(ref_region):
    print(f'outputing region - {ref_region}')
    fields = re.split(':|-', ref_region)
    replace_ref_region = fields[0] + '_' + fields[1] + '-' + fields[2]

    with open(f'{output_dir}/{replace_ref_region}.fa', 'w') as outfile:
        for contig_region, seq, contig_cpname in ref_to_asm_fa[ref_region]:
            outfile.write(f'>{contig_cpname}' + '\n') 
            outfile.write(seq + '\n')       
    return 


if __name__ == "__main__":
    ref_to_asm_fa = defaultdict(list)
    for asm in asms:
        read(asm)

    for ref_region in ref_to_asm_fa:
        output(ref_region)














