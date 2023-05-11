import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 

parser = argparse.ArgumentParser('Filter SNP')
parser.add_argument('-i', '--input', nargs='?', required=True)
parser.add_argument('-o', '--out', nargs='?', required=True)
parser.add_argument('-a', '--asm', nargs='?', required=True)
parser.add_argument('-l', '--indel_length', nargs='?', required=True)

args = parser.parse_args()

input = args.input
out = args.out
asm = args.asm
indel_length = int(args.indel_length)


def transformInputToOutput(input, out):
  with open(input, 'r') as infile, open(out, 'w') as outfile:

    lines = infile.readlines()
    for i, line in enumerate(lines):

      if line.startswith('##'):
        outfile.write(line)
        continue 

      line = line.rstrip('\n')
      fields = line.split('\t')

      if line.startswith('#CHROM'):
        fields[-1] = asm
        outfile.write('\t'.join(fields) + '\n')
        continue 

      REF = fields[3]
      ALT = fields[4]

      if abs(len(REF) - len(ALT)) < indel_length:
        continue 

      # INFO = fields[7]
      # SVLEN = ''
      # for s in INFO.split(';'):
      #   if s.startswith('SVLEN='):
      #     SVLEN = s
      #     break 

      # if SVLEN != '':
      #   SVLEN = int(SVLEN.strip('SVLEN='))
      #   if (0 < SVLEN < 10) or (-10 < SVLEN < 0):
      #     continue 

      # outfile.write('\t'.join(fields) + '\n')
      fields[2] = str(i)
      fields[8] = 'GT'

      if ":" in fields[9]:
        fields[9] = fields[9].split(':')[0]
      if "/" in fields[9]:
        fields[9] = fields[9].replace('/', '|')
        
      outfile.write('\t'.join(fields) + '\n')

  return 

if __name__ == '__main__':
  transformInputToOutput(input, out)



















