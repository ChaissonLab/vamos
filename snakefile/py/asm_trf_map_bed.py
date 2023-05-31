import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 
import re 

parser = argparse.ArgumentParser('Generate the mapped coordinate for each trf interval in asm')
parser.add_argument('-r', '--ref_region_bed_flank_seq', nargs='?', required=True)
parser.add_argument('-a', '--asm_region_bed_flank', nargs='?', required=True)
parser.add_argument('-o', '--output', nargs='?', required=True)
parser.add_argument('-t', '--dataset', nargs='?', required=True)
parser.add_argument('-p', '--asm_hap', nargs='?', required=True)


args = parser.parse_args()
ref_region_bed_flank_seq = args.ref_region_bed_flank_seq
asm_region_bed_flank = args.asm_region_bed_flank
output = args.output
dataset = args.dataset
asm_hap = args.asm_hap

class Intv(object):
    def __init__(self, contig, start, end, id):
        self.contig = contig
        self.start = start 
        self.end = end 
        self.id = id
        return 

def read_asm_bed(input):
  intv_dict = defaultdict()
  flank_dict = defaultdict()

  try:

    with open(input + ".full", 'r') as file:
      for line1 in file:
        line1 = line1.rstrip('\n')
        fields1 = line1.split('\t')
        contig1, s_up, e_up = re.split(":|-", fields1[0])
        direction1, id_up = fields1[1], fields1[2]
        assert(direction1 == 'upstream'), 'ERROR: wrong direction!'

        line2 = next(file)
        line2 = line2.rstrip('\n')
        fields2 = line2.split('\t')
        contig2, s_down, e_down = re.split(":|-", fields2[0])
        direction2, id_down = fields2[1], fields2[2]
        assert(direction2 == 'downstream'), 'ERROR: wrong direction!'

        assert(contig1 == contig2), "ERROR: inconsistent contig name!"
        assert(id_up == id_down), "ERROR: inconsistent direction!"

        contig = f'{contig1}:{e_up}-{s_down}'
        upcoord = fields1[0]
        downcoord = fields2[0]
        intv_dict[contig] = {'up': [upcoord, 0, False, ''], 'down': [downcoord, 0, False, '']} 

        flank_dict[upcoord] = ['up', contig]
        flank_dict[downcoord] = ['down', contig]

  except StopIteration:
    print("odd # of rows!")

  return intv_dict, flank_dict


def read_ref_bed(input, flank_dict, intv_dict):
  # print(flank_dict)

  with open(input, 'r') as file:
    lines = file.readlines()
    for line in lines:
      line = line.rstrip('\n')
      fields = line.split('\t')
      contig, flag, chrom, ref_start, ref_end, mapq, cigar = fields 
      flag, mapq = int(flag), int(mapq)

      assert(contig in flank_dict), 'ERROR: contig not in intv_dict!'
      direction, intv = flank_dict[contig]
      assert(intv_dict[intv][direction][0] == contig), "ERROR: inconsistent!"

      intv_dict[intv][direction][1] = 0 if flag == 0 else 1 
      intv_dict[intv][direction][2] = True
      intv_dict[intv][direction][3] = f'{chrom}:{ref_start}-{ref_end}'

  return 

def output_intv_and_mapped_intv(output, intv_dict):
    with open(output, 'w') as file:

      for intv in intv_dict:
        upcoord, up_di, up_state, up_ref_region = intv_dict[intv]['up']
        downcoord, down_di, down_state, down_ref_region = intv_dict[intv]['down']

        if (up_state and down_state) == False:
          continue 

        up_chrom, up_ref_start, up_ref_end = re.split(":|-", up_ref_region) 
        up_ref_start, up_ref_end = int(up_ref_start), int(up_ref_end)

        down_chrom, down_ref_start, down_ref_end = re.split(":|-", down_ref_region)
        down_ref_start, down_ref_end = int(down_ref_start), int(down_ref_end)  

        if dataset == "hprc":
          if ((up_di or down_di) == 0) and (up_chrom == down_chrom) and (up_ref_end < down_ref_start):
            file.write('\t'.join([up_chrom, str(up_ref_end), str(down_ref_start), intv, '0']) + '\n')
          elif   ((up_di & down_di) == 1) and (up_chrom == down_chrom) and (up_ref_start > down_ref_end):
            file.write('\t'.join([up_chrom, str(down_ref_end), str(up_ref_start), intv, '1']) + '\n')
        elif dataset == "hgsvc":
          # add sample and hap name to intv
          sample, hap = asm_hap.split('_')

          new_intv = f"{sample}#1#{intv}" if hap == 'maternal' else f"{sample}#2#{intv}"

          if ((up_di or down_di) == 0) and (up_chrom == down_chrom) and (up_ref_end < down_ref_start):
            file.write('\t'.join([up_chrom, str(up_ref_end), str(down_ref_start), new_intv, '0']) + '\n')
          elif   ((up_di & down_di) == 1) and (up_chrom == down_chrom) and (up_ref_start > down_ref_end):
            file.write('\t'.join([up_chrom, str(down_ref_end), str(up_ref_start), new_intv, '1']) + '\n')


if __name__ == "__main__":
  intv_dict, flank_dict = read_asm_bed(asm_region_bed_flank)
  print("finish read_asm_bed")
  read_ref_bed(ref_region_bed_flank_seq, flank_dict, intv_dict)
  print("finish read_ref_bed")
  output_intv_and_mapped_intv(output, intv_dict)












