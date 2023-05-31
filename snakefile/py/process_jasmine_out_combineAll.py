import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 
import re 

"""
This script:
1. convert genotype column to "GT"
2. convert "/" to "|: in genotype
3. convert "./.:NA:NA:NA:NA" to 0|0

eg. GT:IS:OT:DV:DR  ./.:NA:NA:NA:NA 1|.:.:INS:0:.
"""

parser = argparse.ArgumentParser('Process the jasmine merged VCF to form a real multi-sample VCF')
# parser.add_argument('-i', '--input_vcf_dir', nargs='?', required=True)
parser.add_argument('-v', '--input_merged_vcf', nargs='?', required=True)
parser.add_argument('-f', '--input_file_list', nargs='?', required=True)
parser.add_argument('-o', '--output_merged_vcf', nargs='?', required=True)

args = parser.parse_args()
# input_vcf_dir = args.input_vcf_dir
input_merged_vcf = args.input_merged_vcf
output_merged_vcf = args.output_merged_vcf
input_file_list = args.input_file_list


def readFileList(input):
  single_vcfs = []
  samples = []
  with open(f"{input}", 'r') as file:
    lines = file.readlines()
    for line in lines:
      line = line.rstrip('\n')
      single_vcfs.append(line)

      name = line.split('.')[0]
      samples.append(name)

  return single_vcfs, samples


def readMergedVCF(input, sample2svs):

  SVlines = []
  header = []
  with open(f"{input}", 'r') as file:
    lines = file.readlines()
    for line in lines:
      line = line.rstrip('\n')
      if line.startswith('##'):
        header.append(line)
        continue 
      if line.startswith('#CHROM'):
        continue

      fields = line.split('\t')

      INFO = fields[7]

      fields[8] = 'GT'
      
      for i in range(9, len(fields)):
        gt = fields[i]
        if ':' in gt:
          fields[i] = gt.split(":")[0]
        if '/' in gt:
          fields[i] = fields[i].replace('/', '|')

        if fields[i] == '.|.' or fields[i] == '.|0' or fields[i] == '0|.':
          fields[i] = '0|0'
        elif fields[i] == '.|1':
          fields[i] = '0|1'
        elif fields[i] == '1|.':
          fields[i] = '1|0'
        elif fields[i] == '.|2':
          fields[i] = '0|2'
        elif fields[i] == '2|.':
          fields[i] = '2|0'          

      SVlines.append('\t'.join(fields))

  return SVlines, header 


def outputMergedVCF(output_merged_vcf, SVlines, header):

  with open(f"{output_merged_vcf}", 'w') as file:

    # write header
    header.pop() # remove the final line 

    for h in header:
      file.write(h + '\n')

    line = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
    file.write('\t'.join(line) + '\n')

    # write body 
    for line in SVlines:
      file.write(line + '\n')

  return 


if __name__ == "__main__":
  single_vcfs, samples = readFileList(input_file_list)
  num_samples = len(samples)

  sample2svs = defaultdict(lambda: defaultdict()) # sample -> {ID -> SV info (GT)}
  SVlines, header  = readMergedVCF(input_merged_vcf, sample2svs)

  outputMergedVCF(output_merged_vcf, SVlines, header)



















