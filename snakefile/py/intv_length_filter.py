import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 


parser = argparse.ArgumentParser('Filter VNTR entries by the comparison between ref_len and avg_asm_len')
parser.add_argument('-i', '--input', nargs='?', required=True)
parser.add_argument('-o', '--output', nargs='?', required=True)

args = parser.parse_args()
input = args.input
output = args.output

def readBed(input):
    intv_list = []
    with open(input, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.rstrip('\n')
            fields = line.split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            asm_info_list = [asm.split(':')[-1] for asm in fields[3].split(',')]
            asm_len_list = [int(asm.split('-')[1]) - int(asm.split('-')[0]) for asm in asm_info_list]
            avg = sum(asm_len_list) / len(asm_len_list)
            intv_list.append((end - start, avg, np.var(asm_len_list), (end - start) / avg, line))
    return intv_list

def output_intv(intv_list):
    with open(output, 'w') as file:
        for ref_len, avg_asm_len, len_var, len_per, line in intv_list:
            if len_per <= 0.5:
              file.write(line + '\n')


if __name__ == '__main__':
  intv_list = readBed(input)
  output_intv(intv_list)














