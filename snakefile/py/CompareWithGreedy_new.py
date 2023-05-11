#!~/.conda/envs/nb/bin/python3

import os
import sys
import re
from collections import defaultdict
import Levenshtein as lv
import pandas as pd
import argparse

parser = argparse.ArgumentParser('Compare greedy annotation and vamos annotation for one assembly sample')
parser.add_argument('-i', '--input_asm_name', nargs='?', required=True)
parser.add_argument('-o', '--output_bed_file', nargs='?', required=True)
parser.add_argument('-m', '--q_mode', nargs='?', required=True)
parser.add_argument('-v', '--vamos_anno_path', nargs='?', required=True)
parser.add_argument('-g', '--greedy_anno_path', nargs='?', required=True)
parser.add_argument('-l', '--asm_liftover_path', nargs='?', required=True)
parser.add_argument('-f', '--original_vntr_motif_path', nargs='?', required=True)
parser.add_argument('-p', '--original_vntr_motif_filter', nargs='?', required=True)

args = parser.parse_args()

input_asm_name = args.input_asm_name
output_bed_file = args.output_bed_file
q_mode = args.q_mode
vamos_anno_path = args.vamos_anno_path
greedy_anno_path = args.greedy_anno_path 
asm_liftover_path = args.asm_liftover_path
original_vntr_motif_filter = int(args.original_vntr_motif_filter)
original_vntr_motif_path = args.original_vntr_motif_path

def readBed(original_vntr_motif_path, vntr_set):

	with open(original_vntr_motif_path, 'r') as fin:

		lines = fin.readlines()

		for line in lines:

			if line.startswith('#'): continue 

			line = line.rstrip('\n')
			fields = re.split('\t', line)
			chrom, start, end = fields[0], fields[1], fields[2]
			motif_list = fields[3]
			vntr = f"{chrom}_{start}-{end}"
			if len(motif_list.split(',')) >= original_vntr_motif_filter:
				vntr_set.add(vntr)

	return 

def readVCF(vcf_in_file, asm_dict, asm, asm_liftover_path):

	# read asm annotation vcf
	with open(vcf_in_file, 'r') as fin:

		lines = fin.readlines()

		for idx, line in enumerate(lines):

			if line.startswith('#'): continue

			line = line.rstrip("\n")
			fields = re.split('\t', line)
			chrom, start = fields[0], fields[1]
			infos = re.split(';', fields[7])

			assert infos[0][:4] == 'END='
			assert infos[1][:3] == 'RU='
			assert infos[3][:11] == 'ALTANNO_H1='

			end = infos[0][4:]
			motif_list = list(re.split(',', infos[1][3:]))
			anno_index = [int(i) for i in (infos[3][11:]).split(',')]
			anno_str = ''.join([motif_list[i] for i in anno_index])
			vntr = f"{chrom}_{start}-{end}"

			if vntr not in vntr_set: continue 
			if (len(motif_list) == 1): continue
			asm_dict[vntr]['motif'] = motif_list
			asm_dict[vntr]['annoIndex'] = anno_index
			asm_dict[vntr]['annoStr'] = anno_str

	print(f"finish reading VCF")

	# read liftover seq fasta file
	asm_fa_file = f"{asm_liftover_path}"
	with open(asm_fa_file, 'r') as fin:

		lines = fin.readlines()
		nlines = len(lines)

		# read 2 lines at a time
		for idx in range(0, nlines, 2):

			line = lines[idx]
			line = line.rstrip('\n')
			assert line.startswith('>')
			read_name = line[1:]
			vntr = read_name.split(';')[1]

			if vntr not in asm_dict: continue
			if vntr not in vntr_set: continue 


			line = lines[idx + 1]
			line = line.rstrip('\n')
			read_seq = line 
			anno_str = asm_dict[vntr]['annoStr']

			if (len(read_seq) > 20000 and len(anno_str) > 20000): 
				print(vntr, len(read_seq), len(anno_str))
			if len(read_seq) > 2000000 and len(anno_str) > 2000000: 
				del asm_dict[vntr]
				continue

			d = lv.distance(read_seq, anno_str)
			asm_dict[vntr]['seq'] = read_seq
			asm_dict[vntr]['dist'] = d

	print(f"finish reading liftover fa")

	return 


if __name__ == "__main__":

	# read orignal vntr motif bed 
	vntr_set = set()
	readBed(original_vntr_motif_path, vntr_set)

	# asm -> {vntr -> {"motif": motif_list, "annoIndex": [annotation index], "anno_str": annotation nucl string, "seq": read seq, "dist": the edit distance between read_seq and annotation string}}
	asm_emotif = defaultdict(lambda: defaultdict())
	asm_greedy = defaultdict(lambda: defaultdict()) 

	readVCF(vamos_anno_path, asm_emotif, input_asm_name, asm_liftover_path)
	readVCF(greedy_anno_path, asm_greedy, input_asm_name, asm_liftover_path)

	gt_dist = defaultdict(list) # {'greedy': [], "greedy_ratio": greedy_dist / read_seq_len, q_mode: [], 'vntr': []}

	for vntr in asm_emotif:
		if 'dist' not in asm_emotif[vntr]:
			continue
		gt_dist['vntr'].append(vntr)
		gt_dist[q_mode].append(asm_emotif[vntr]['dist'])
		gt_dist['greedy'].append(asm_greedy[vntr]['dist'])
		r1 = asm_emotif[vntr]['dist'] / len(asm_emotif[vntr]['seq'] )
		gt_dist[q_mode + '_ratio'].append(min(r1, 1.0))
		r2 = asm_greedy[vntr]['dist'] / len(asm_greedy[vntr]['seq'])
		gt_dist['greedy_ratio'].append(min(r2, 1.0))
		gt_dist['seq_len'].append(len(asm_emotif[vntr]['seq']))


	df = pd.DataFrame.from_dict(gt_dist)

	"""
	emotif annotation is better
	"""
	emotif_better = list(df.loc[(df[q_mode] < df['greedy'])]['vntr'])
	emotif_better_10 = list(df.loc[(df[q_mode] < 0.8 * df['greedy'])]['vntr'])
	"""
	greedy annotation is better
	"""
	greedy_better = list(df.loc[(df[q_mode] > df['greedy'])]['vntr'])
	greedy_better_10 = list(df.loc[(0.8 * df[q_mode] > df['greedy'])]['vntr'])

	with open(output_bed_file, 'w') as f:
		f.write('\t'.join([input_asm_name, \
						   str(len(df[q_mode])), \
						   str(len(emotif_better)), \
						   str(len(emotif_better_10)), \
						   str(len(greedy_better)), \
						   str(len(greedy_better_10))]) + "\n")
