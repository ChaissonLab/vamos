import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 

parser = argparse.ArgumentParser('Extract flanking coordinate from TRF bed file')
parser.add_argument('-i', '--input', nargs='?', required=True)
parser.add_argument('-f', '--input_fa_path', nargs='?', required=True)
parser.add_argument('-o', '--output', nargs='?', required=True)

args = parser.parse_args()
input = args.input
output = args.output
input_fa_path = args.input_fa_path


class Intv(object):
    def __init__(self, contig, start, end, id):
        self.contig = contig
        self.start = start 
        self.end = end 
        self.id = id
        return 

def readBed(input):
    intv_list = []
    with open(input, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.rstrip('\n')
            fields = line.split('\t')
            contig, start, end, id = fields[0], int(fields[1]), int(fields[2]), fields[3]
            intv_list.append(Intv(contig, start, end, id))
    return intv_list


def readFai(input):
  contig2len = defaultdict(int)
  asm, hap = input.split("/")[-1].split('.')[0].split('_')
  asm_hap_fai = f'{input_fa_path}.fai' 

  with open(asm_hap_fai, 'r') as file:
    lines = file.readlines()
    for line in lines:
      line = line.rstrip('\n')
      fields = line.split('\t')
      contig, length = fields[0], int(fields[1])
      contig2len[contig] = length
  return contig2len


def extract_flank_coordinate(intv, contig2len):
  start, end = intv.start, intv.end 
  length = contig2len[intv.contig]
  start = max(0, start - 5000)
  end = min(length - 1, end + 5000)

  return (intv.contig, intv.id, (start, intv.start), (intv.end, end))

def output_intv(flank_coordinate):
    with open(output, 'w') as file:
        for i, intv in enumerate(flank_coordinate):
            contig, id = intv[0], intv[1]
            upstream_s, upstream_e = intv[2]
            downstream_s, downstream_e = intv[3]
            if upstream_s + 1000 <= upstream_e and downstream_s + 1000 < downstream_e:
              file.write(f'{contig}:{upstream_s}-{upstream_e}' + '\n')
              file.write(f'{contig}:{downstream_s}-{downstream_e}' +'\n')

    with open(output + '.full', 'w') as file:
        for i, intv in enumerate(flank_coordinate):
            contig, id = intv[0], intv[1]
            upstream_s, upstream_e = intv[2]
            downstream_s, downstream_e = intv[3]
            if upstream_s + 1000 <= upstream_e and downstream_s + 1000 < downstream_e:
              file.write(f'{contig}:{upstream_s}-{upstream_e}' + '\t' + 'upstream' + '\t' + f'{id}' + '\n')
              file.write(f'{contig}:{downstream_s}-{downstream_e}' + '\t' + 'downstream' + '\t' +  f'{id}' +'\n')

if __name__ == "__main__":
  hap2int = {"maternal" : "0", "paternal" : "1"}
  contig2len = defaultdict(int)
  flank = 5000 

  intv_list = readBed(input)
  contig2len = readFai(input)
  flank_coordinate = []

  for intv in intv_list:
    flank_coordinate.append(extract_flank_coordinate(intv, contig2len))

  output_intv(flank_coordinate)














