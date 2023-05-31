import argparse 
from collections import defaultdict, Counter
from os import listdir 
import os 

parser = argparse.ArgumentParser('Summarize the alleles for samples at given VNTR loci')
parser.add_argument('-v', '--vntr_loci_path', nargs='?', required=True)
parser.add_argument('-s', '--vntr_sv_ovp_path', nargs='?', required=True)
parser.add_argument('-a', '--asm_anno_vcf_folder', nargs='?', required=True)
parser.add_argument('-o', '--out_path', nargs='?', required=True)


args = parser.parse_args()

vntr_loci_path = args.vntr_loci_path
asm_anno_vcf_folder = args.asm_anno_vcf_folder.rstrip("/")
vntr_sv_ovp_path = args.vntr_sv_ovp_path
out_path = args.out_path


vntr2allele = defaultdict()
vntr2sv = defaultdict(lambda: defaultdict()) # vntr -> sv -> allels for hap
vntr2svalleles = defaultdict() # vntr -> # of SV alleles

asm_hap_samples = [asm.rstrip(".anno.vcf") for asm in listdir(asm_anno_vcf_folder) if asm.endswith(".anno.vcf")]
print(len(asm_hap_samples))

asm_samples = list(set(['_'.join(asm.split('_')[:2]) for asm in listdir(asm_anno_vcf_folder) if asm.endswith(".anno.vcf")]))
print(asm_samples)
print(len(asm_samples))


def readBED(vntr_loci_path):
    with open(vntr_loci_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            chrom, start, end = line.rstrip('\n').split('\t')
            vntr2allele[f"{chrom}_{start}-{end}"] = {'len': Counter(), 'string': Counter()}


def readAnnoFromAsm(asm_hap_sample):
    asm_file = f'{asm_anno_vcf_folder}/{asm_hap_sample}.anno.vcf'
    print("start: ", asm_file)

    with open(asm_file, 'r') as file:
        lines = file.readlines()
        for line in lines:

            if line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            chrom = fields[0]
            start = fields[1]
            info = fields[7]
            end = info.split(';')[0].strip('END=')

            identifier = f"{chrom}_{start}-{end}"

            anno = ''
            for a in info.split(';'):
                if a.startswith('ALTANNO_H1='):
                    anno = a

            assert(identifier in vntr2allele), "identifier not in vntr2allele dict"

            vntr2allele[identifier]['len'][len(anno[11:].split(','))] += 1
            vntr2allele[identifier]['string'][anno[11:]] += 1

    return 


def readVCF2SV(vntr_sv_ovp_path):

    with open(vntr_sv_ovp_path) as file:
        lines = file.readlines()
        for line in lines:
            line_sp = line.rstrip('\n').split('\t')
            vntr_chrom, vntr_start, vntr_end, sv_chrom, sv_start, sv_end, sv_type, sv_rid, sv_len = line_sp[:9]

            sv_alleles = line_sp[9:-1] # the last pos is ovp
            if len(sv_alleles) != len(asm_samples):
                print(len(sv_alleles), len(asm_samples))
                raise AssertionError()

            identifier = vntr_chrom + '_' + vntr_start + '-' + vntr_end
            sv = sv_chrom + '_' + sv_start + '-' + sv_end + ',sv_type=' + sv_type + ',sv_rid=' + sv_rid
            vntr2sv[identifier][sv] = {'sv_alleles' : sv_alleles, 'sv_len' : int(sv_len)}


    # count SV alleles, since there can be multiple SVs at same VNTR locus
    for vntr in vntr2sv:
        sv_seen = Counter()
        len_seen = Counter()
        sv_seen_sample = defaultdict(list)

        for i in range(len(asm_samples)):
            allele_sv = [[], []]
            allele_len = [0, 0]

            for sv in vntr2sv[vntr]:
                sv_len = vntr2sv[vntr][sv]['sv_len']
                if i >= len(vntr2sv[vntr][sv]['sv_alleles']) or '|' not in vntr2sv[vntr][sv]['sv_alleles'][i]:
                    print(i, len(vntr2sv[vntr][sv]['sv_alleles']), vntr, sv, vntr2sv[vntr][sv]['sv_alleles'])
                    raise AssertionError()
                hap = vntr2sv[vntr][sv]['sv_alleles'][i].split('|')

                allele_sv[0].append(hap[0])
                allele_sv[1].append(hap[1])

                for k in range(2):
                    if hap[k] == '1':
                        allele_len[k] += sv_len

            for k in range(2):
                if '.' in allele_sv[k]:
                    continue
                sv_seen[','.join(allele_sv[k])] += 1
                sv_seen_sample[','.join(allele_sv[k])].append(asm_samples[i])
                len_seen[allele_len[k]] += 1

        sv_allele_string = []
        sv_allele_string_cnt = []
        sv_seen_sample_string = []
        for k, v in sv_seen.items():
            sv_allele_string.append(k)
            sv_allele_string_cnt.append(str(v))
            sv_seen_sample_string.append(','.join(sv_seen_sample[k]))

        sv_len_allele_string = []
        sv_len_allele_string_cnt = []
        for k, v in len_seen.items():
            sv_len_allele_string.append(str(k))
            sv_len_allele_string_cnt.append(str(v))   

        vntr2svalleles[vntr] = {'sv_allele_cnt' : len(sv_seen), \
                                'sv_len_cnt' : len(len_seen), \
                                "sv_allele_string" : ';'.join(sv_allele_string), \
                                "sv_allele_string_cnt" : ';'.join(sv_allele_string_cnt), \
                                "sv_len_allele" : ';'.join(sv_len_allele_string), \
                                "sv_len_allele_cnt" : ';'.join(sv_len_allele_string_cnt), \
                                "sv_allele_string_sample" : ';'.join(sv_seen_sample_string)}

    return 


def outputVCF_highvamosanno_lowsvalleles(out_path):
    folder_path = '/' + '/'.join(out_path.split('/')[1:-1])
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    with open(out_path, 'w') as file:

        for vntr in vntr2svalleles:

            assert(vntr in vntr2allele), "vntr not in vntr2allele dict"

            numSValleles = vntr2svalleles[vntr]['sv_allele_cnt'] 
            numSVLenalleles = vntr2svalleles[vntr]['sv_len_cnt']

            numLenAnno = len(vntr2allele[vntr]['len']) 
            numStringAnno = len(vntr2allele[vntr]['string']) 
                
            vntr_anno_len = []
            vntr_anno_len_cnt = []
            for k, v in vntr2allele[vntr]['len'].items():
                vntr_anno_len.append(str(k))
                vntr_anno_len_cnt.append(str(v))


            vntr_anno_string = []
            vntr_anno_string_cnt = []
            for k, v in vntr2allele[vntr]['string'].items():
                vntr_anno_string.append(k)
                vntr_anno_string_cnt.append(str(v))

            file.write('\t'.join([vntr, \
                                     str(numSValleles), str(numSVLenalleles), \
                                     str(numLenAnno), str(numStringAnno), \
                                     ';'.join(vntr2sv[vntr]), \
                                     'anno_len=' + ';'.join(vntr_anno_len), \
                                     'anno_cnt=' + ';'.join(vntr_anno_len_cnt), \
                                     'anno_string=' + ';'.join(vntr_anno_string), \
                                     'anno_string_cnt=' + ';'.join(vntr_anno_string_cnt), \
                                     'sv_allele_string=' + vntr2svalleles[vntr]['sv_allele_string'], \
                                     'sv_allele_string_cnt=' + vntr2svalleles[vntr]['sv_allele_string_cnt'], \
                                     'sv_len_allele=' + vntr2svalleles[vntr]['sv_len_allele'], \
                                     'sv_len_allele_cnt=' + vntr2svalleles[vntr]['sv_len_allele_cnt'], \
                                     'sv_allele_samples=' + vntr2svalleles[vntr]['sv_allele_string_sample']]) + '\n')

    return 


if __name__ == '__main__':
    readBED(vntr_loci_path)

    for asm_hap_sample in asm_hap_samples:
        print(asm_hap_sample)
        readAnnoFromAsm(asm_hap_sample)

    readVCF2SV(vntr_sv_ovp_path)
    outputVCF_highvamosanno_lowsvalleles(out_path)


























