import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 

parser = argparse.ArgumentParser('Filter SNP')
parser.add_argument('-i', '--input', nargs='?', required=True)
parser.add_argument('-o', '--out', nargs='?', required=True)

args = parser.parse_args()

input = args.input
out = args.out


def transformInputToOutput(input, out):
  with open(input, 'r') as infile, open(out, 'w') as outfile:

    lines = infile.readlines()
    for i, line in enumerate(lines):
      if line.startswith('#'):
        continue 
      line = line.rstrip('\n')
      fields = line.split('\t')
      INFO = fields[7]

      # if 'SVTYPE=SNV' in INFO:
      #   continue 

      SVLEN = ''
      for s in INFO.split(';'):
        if s.startswith('SVLEN='):
          SVLEN = s
          break 
      SVLEN = int(SVLEN.strip('SVLEN='))

      # if (0 < SVLEN < 10) or (-10 < SVLEN < 0):
      #   continue 

      if 'INS' in INFO or SVLEN > 0:
        new_line = [fields[0], fields[1], fields[1], "INS", str(i), str(SVLEN)] + fields[9:]
        outfile.write('\t'.join(new_line) + '\n')
      elif 'DEL' in INFO or SVLEN < 0:
        new_line = [fields[0], fields[1], str(int(fields[1]) + abs(SVLEN)), "DEL", str(i), str(SVLEN)] + fields[9:]
        outfile.write('\t'.join(new_line) + '\n')
  return 

if __name__ == '__main__':
  transformInputToOutput(input, out)



















