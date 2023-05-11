import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 

parser = argparse.ArgumentParser('Merge intervals in a TRF Bed file')
parser.add_argument('-i', '--input', nargs='?', required=True)
parser.add_argument('-o', '--output', nargs='?', required=True)

args = parser.parse_args()
input = args.input
output = args.output


class Intv(object):
    def __init__(self, chrom, start, end, contig_region_list, strand_list):
        self.contig_region_list = contig_region_list
        self.strand_list = strand_list
        self.start = start 
        self.end = end 
        self.chrom = chrom
        return 


def readBed(input):
    intv_list = []

    with open(input, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.rstrip('\n')
            fields = line.split('\t')
            chrom, start, end, contig_region_list, strand_list, dataset = fields[0], \
                                                                          int(fields[1]), \
                                                                          int(fields[2]), \
                                                                          fields[3], \
                                                                          fields[4], \
                                                                          fields[9]
            intv_list.append(Intv(chrom, start, end, \
                                  [f"{dataset}_{contig}" for contig in contig_region_list.split(',')], \
                                  [strand for strand in strand_list.split(',')]))
    return intv_list


def merge(intervals):

    res = []
    for intv in intervals:
        if not res:
            res.append(intv)
            continue 

        if res[-1].chrom != intv.chrom:
            res.append(intv)
        elif res[-1].end >= intv.start:
            res[-1].end = max(res[-1].end, intv.end)
            res[-1].contig_region_list.extend(intv.contig_region_list)
            res[-1].strand_list.extend(intv.strand_list)
        else:
            res.append(intv)

    return res 


def output_intv(intervals):

    with open(output, 'w') as file:

        for i, intv in enumerate(intervals):

            haps = defaultdict(int)
            for contig in intv.contig_region_list:
                haps[contig.split(':')[0]] += 1 

            haps_str, cnt_str = [], []
            for contig in sorted(haps.keys()):
                cnt = haps[contig]
                haps_str.append(contig)
                cnt_str.append(str(cnt))

            file.write('\t'.join([intv.chrom, \
                                 str(intv.start), \
                                 str(intv.end), \
                                 ','.join(intv.contig_region_list), \
                                 ','.join([str(strand) for strand in intv.strand_list]), \
                                 ','.join(haps_str), \
                                 ','.join(cnt_str), str(len(haps_str)), str(intv.end - intv.start)]) + '\n')


if __name__ == "__main__":
    intv_list = readBed(input)
    merged_intv_list = merge(intv_list)
    output_intv(merged_intv_list)














