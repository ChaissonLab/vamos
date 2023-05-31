import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 
import re 

parser = argparse.ArgumentParser('Extract intervals for asm and ref, and arrange them into a bed file.')
parser.add_argument('-i', '--input', nargs='?', required=True)
parser.add_argument('-o', '--output_bed', nargs='?', required=True)
parser.add_argument('-f', '--output_full_bed', nargs='?', required=True)
parser.add_argument('-a', '--asm_sample', nargs='?', required=True)
parser.add_argument('-d', '--dataset', nargs='?', required=True)

args = parser.parse_args()
input = args.input
output_bed = args.output_bed
output_full_bed = args.output_full_bed
input_dataset = args.dataset

asm_sample = args.asm_sample
asm_name, asm_hap = asm_sample.split('_')[0], asm_sample.split('_')[1]
if input_dataset == 'hprc':
  replace_asm_name = input_dataset + '_' + asm_name + '#' + ('1' if asm_hap == 'paternal' else '2')
elif input_dataset == 'hgsvc':
  replace_asm_name = input_dataset + '_' + asm_name + '#' + ('2' if asm_hap == 'paternal' else '1')


class Asm_Intv_Info(object):
    def __init__(self, contig, start, end, ref_region, strand, asm):
        self.contig = contig
        self.strand = strand
        self.start = start 
        self.end = end 
        self.ref_region = ref_region 
        self.contig_region = self.contig + ":" + str(start) + "-" + str(end);
        self.contig_cpname = ';'.join([self.contig_region, self.ref_region, str(self.strand)])
        self.asm = asm 
        return 


def readBed(input):
    asm_to_ref_intv_dict = defaultdict(list)  # asm -> [asm_intv_info]

    with open(input, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.rstrip('\n')
            fields = line.split('\t')
            chrom, start, end = fields[0], fields[1], fields[2]
            asm_region_list = list(fields[3].split(','))
            asm_strand = list(fields[4].split(','))
            ref_region = f'{chrom}:{start}-{end}'

            assert(len(asm_region_list) == len(asm_strand)), "ERROR: different lengths!"

            for i in range(len(asm_region_list)):
                asm_region = asm_region_list[i]
                strand = asm_strand[i]
                coordinate = asm_region.split(':')[1]
                dataset = asm_region.split('_')[0]
                contig = '_'.join((asm_region.split(':')[0]).split('_')[1:])
                start, end = coordinate.split('-')
                asm = '#'.join(contig.split('#')[:2])

                if dataset == 'hgsvc':
                    contig = contig.split('#')[-1]

                asm_intv_info = Asm_Intv_Info(contig, start, end, ref_region, asm_strand[i], asm)
                asm_to_ref_intv_dict[dataset + '_' + asm].append(asm_intv_info)

    return asm_to_ref_intv_dict


def touch(path):
    with open(path, 'a'):
        os.utime(path, None)


def output(asm_to_ref_intv_dict):

    if not asm_to_ref_intv_dict[replace_asm_name]:
        touch(output_full_bed)
        touch(output_bed)
        return 


    print(f'outputing asm - {asm_sample}')

    with open(output_full_bed, 'w') as file:
        
        for asm_intv_info in asm_to_ref_intv_dict[replace_asm_name]:
            file.write('\t'.join([asm_intv_info.contig_region, \
                                  asm_intv_info.contig_cpname, \
                                  asm_intv_info.strand, \
                                  asm_intv_info.ref_region]) + '\n')


    print(f'outputing asm - {asm_sample}')

    with open(output_bed, 'w') as file:
        
        for asm_intv_info in asm_to_ref_intv_dict[replace_asm_name]:
            file.write(asm_intv_info.contig_region + '\n')

    return 


if __name__ == "__main__":
    asm_to_ref_intv_dict = readBed(input)
    output(asm_to_ref_intv_dict)














