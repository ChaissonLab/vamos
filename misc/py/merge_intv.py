import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 

parser = argparse.ArgumentParser('Merge intervals in a TRF Bed file')
parser.add_argument('-i', '--input', nargs='?', required=True)
parser.add_argument('-o', '--output', nargs='?', required=True)
# parser.add_argument('-p', '--out_path', nargs='?', required=True)

args = parser.parse_args()
input = args.input
output = args.output
# out_path = args.out_path


class Intv(object):
    def __init__(self, contig, start, end):
        self.contig = contig
        self.start = start 
        self.end = end 
        return 

def readBed(input):
    intv_list = []
    with open(input, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.rstrip('\n')
            fields = line.split('\t')
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            intv_list.append(Intv(contig, start, end))
    return intv_list


def merge(intervals):

    res = []
    for intv in intervals:
        if not res:
            res.append(intv)
            continue 

        if res[-1].contig != intv.contig:
            res.append(intv)
        elif res[-1].end >= intv.start:
            res[-1].end = max(res[-1].end, intv.end)
        else:
            res.append(intv)

    return res 


def output_intv(intervals):
    with open(output, 'w') as file:
        for i, intv in enumerate(intervals):
            if intv.end - intv.start <= 10000:
              file.write('\t'.join([intv.contig.split('/')[0], str(intv.start), str(intv.end), "ID=" + str(i)]) + '\n')

# def output_distributed_intv(intervals):
#     for i, intv in enumerate(intervals):
#         with open(out_path + "/interval_" + str(i) + ".bed", 'w') as file:
#             file.write('\t'.join([intv.contig, str(intv.start), str(intv.end)]) + '\n')
    
intv_list = readBed(input)
merged_intv_list = merge(intv_list)
output_intv(merged_intv_list)
# output_distributed_intv(merged_intv_list)














