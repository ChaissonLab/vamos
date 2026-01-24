import os
import sys
import logging
import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

inBed1 = sys.argv[1]
inBed2 = sys.argv[2]
outBed = sys.argv[3]


bed1Dict, bed2Dict = {}, {}
with open(inBed1) as f:
    for line in f:
        fields = line.strip().split('\t')
        chrom,start,end = fields[:3]
        bed1Dict[(chrom,start,end)] = []

with open(inBed2) as f:
    for line in f:
        fields = line.strip().split('\t')
        chrom,start,end = fields[:3]
        bed2Dict[(chrom,start,end)] = f'{chrom}_{start}_{end}'

bed1Dict = bedlib.bedIntersectTag(bed1Dict, bed2Dict)
out = open(outBed, 'w')
for (chrom,start,end),tags in bed1Dict.items():
    if tags == 'NULL': continue
    out.write(f'{chrom}\t{start}\t{end}\t{tags}\n')
out.close()

