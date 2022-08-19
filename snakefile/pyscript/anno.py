# -*- coding: UTF-8 -*-

import os, sys, getopt
import pyabpoa as pa
import edlib as ed
from multiprocessing import Pool

def parse(argv):

    try:
        opts, args = getopt.getopt(argv, 'hs:v:f:', ['sample=', 'vamos=', 'filter='])
    except getopt.GetoptError as err:
        print(err)
        sys.stdout.write('\nanno_cons_first.py -s <sample> -v <vamos> -f <filter>\n'); sys.stdout.flush()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            sys.stdout.write('\nanno_cons_first.py -s <sample> -v <vamos> -f <filter>\n'); sys.stdout.flush()
            sys.exit()

    for opt, arg in opts:
        if opt in ('-s', '--sample'):
            sample = arg
        elif opt in ('-v', '--vamos'):
            vamos = arg
        elif opt in ('-f', '--filter'):
            filter = arg.split(',')

    return(sample, vamos, filter)


def readvcf(vcf):
    with open(vcf) as f:
        for line in f:
            if line.startswith('#'): continue

            chr, pos, id, ref, alt, qual, filter, info, format, genotype = line.strip().split('\t')
            end, ru, svtype, altanno, len = info.split(';')[:-1]
            annoDict[(chr,pos,end)] = ['\t'.join([chr, pos, id, ref, alt, qual, filter]), ';'.join([end, ru, svtype]), ';'.join([altanno, len]), ';'.join([altanno, len])]


def anno(key):

    if not os.path.exists('analysis/%s/anno/%s' %(sample, key)):
        os.makedirs('analysis/%s/anno/%s' %(sample, key))

    vntr = 'analysis/split/vntrs.%s.bed' %(key)
    bam = 'analysis/%s/phase/%s/%s.phase.bam' %(sample, key, key)
    anno = 'analysis/%s/anno/%s/anno.bed' %(sample, key)

    fa = 'analysis/%s/anno/%s/tmp.anno.fa' %(sample, key)
    bed = 'analysis/%s/anno/%s/tmp.anno.bed' %(sample, key)
    csv = 'analysis/%s/anno/%s/tmp.anno.csv' %(sample, key)
    vcf = 'analysis/%s/anno/%s/tmp.anno.vcf' %(sample, key)

    # vamos
    os.system('%s --output_read --per_read -b %s -r %s -s %s -o %s' %(vamos, bam, vntr, sample, anno))

    annoDict = {}

    with open(anno) as f:

        for index,line in enumerate(f):

            try:
                chr,start,end,temp,motifs,fields = line.strip().split('\t')
            except:
                continue

            reads = fields.split(';')[:-1]
            ids = [r.split(':')[0] for r in reads]
            groups = [r.split(':')[1] for r in reads]
            annos = [r.split(':')[3] for r in reads]
            seqs = [r.split(':')[4] for r in reads]

            seqs1 = [ seqs[i] for i,g in enumerate(groups) if g == '1' ]
            seqs2 = [ seqs[i] for i,g in enumerate(groups) if g == '2' ]

            # filter by count of seqs
            if len(seqs1) < filterDepth or len(seqs2) < filterDepth:
                sys.stdout.write('Warnings: skipping locus of poor coverage, %s:%s-%s\n' %(chr,start,end)); sys.stdout.flush()
                continue

            if sum([len(s) for s in seqs1]) > 1000000 or sum([len(s) for s in seqs2]) > 1000000:
                sys.stdout.write('Warnings: skipping sequences unhandable by MSA, %s:%s-%s\n' %(chr,start,end)); sys.stdout.flush()
                continue

            cons1 = pa.msa_aligner().msa(seqs1, out_cons=True, out_msa=False).cons_seq[0]
            cons2 = pa.msa_aligner().msa(seqs2, out_cons=True, out_msa=False).cons_seq[0]

            # filter by sequence similarities for consensus
            #distances1 = [ed.align(cons1, seq)['editDistance'] for seq in seqs1]
            #distances2 = [ed.align(cons2, seq)['editDistance'] for seq in seqs2]
            #if ( sum(distances1) / len(distances1) ) / len(cons1) >= filterDist or ( sum(distances2) / len(distances2) ) / len(cons2) >= filterDist: continue

            # annotation of the consensus sequence
            out = open(fa, 'w')
            out.write('>cons\n' + cons1 +'\n')
            out.close()
            out = open(bed, 'w')
            out.write('\t'.join([chr,start,end]) +'\n')
            out.close()
            out = open(csv, 'w')
            out.write(motifs +'\n')
            out.close()
            os.system('%s --conseq_anno -i %s -v %s -m %s -s %s -o %s' %(vamos, fa, bed, csv, sample, vcf))

            with open(vcf) as f:
                for l in f:
                    if l.startswith('#'): continue

                    chr, pos, id, ref, alt, qual, filter, info, format, genotype = l.strip().split('\t')
                    end, ru, svtype, altanno_h1, len_h1 = info.split(';')[:-1]
                    annoDict[(chr,pos,end)] = ['\t'.join([chr, pos, id, ref, alt, qual, filter]), ';'.join([end, ru, svtype]), ';'.join([altanno_h1, len_h1]), ';'.join([altanno_h1, len_h1])]

            out = open(fa, 'w')
            out.write('>cons\n' + cons2 +'\n')
            out.close()
            os.system('%s --conseq_anno -i %s -v %s -m %s -s %s -o %s' %(vamos, fa, bed, csv, sample, vcf))

            with open(vcf) as f:
                for l in f:
                    if l.startswith('#'): continue

                    chr, pos, id, ref, alt, qual, filter, info, format, genotype = l.strip().split('\t')
                    end, ru, svtype, altanno_h2, len_h2 = info.split(';')[:-1]
                    if not (chr,pos,end) in annoDict:
                        annoDict[(chr,pos,end)] = ['\t'.join([chr, pos, id, ref, alt, qual, filter]), ';'.join([end, ru, svtype]), ';'.join([altanno_h2, len_h2]), ';'.join([altanno_h2, len_h2])]
                    annoDict[(chr,pos,end)][3] = ';'.join([altanno_h2, len_h2])

    if not os.path.exists(vcf): return

    with open(vcf) as f: header = [l for l in f if l.startswith('#')]
    out = open(vcf, 'w')
    for h in header: out.write(h)
    for coor,text in annoDict.items():
        #if text[2] == text[3]: # homozygous by anno&length
        if text[2].split(';')[0] == text[3].split(';')[0]: # homozygous by anno
            out.write('\t'.join([text[0], text[1]+';'+text[2]+';', 'GT', '1/1']) +'\n')
        else: # heterozygous
            out.write('\t'.join([text[0], text[1]+';'+text[2]+';'+text[3].replace('ALTANNO_H1','ALTANNO_H2')+';', 'GT', '0/1']) +'\n')

    out.close()
    os.system('rm %s %s %s' %(sample, fa, bed, csv))
    os.system('touch analysis/%s/anno/%s/done' %(sample, key))

def multi():

    keys = [f.split('.')[1] for f in os.listdir('analysis/split') if 'vntrs.' in f]
    if not os.path.exists('analysis/%s/anno' %(sample)):
        os.makedirs('analysis/%s/anno' %(sample))

    with Pool(16) as p:
        p.map(anno, keys)

    texts = []
    for key in keys:
        v = 'analysis/%s/anno/%s/tmp.anno.vcf' %(sample, key)

        if not os.path.exists(v): continue

        with open(v) as f:
            for line in f:
                if not line.startswith('#'): texts.append(line)

    with open(v) as f: header = [l for l in f if l.startswith('#')]

    out = open('analysis/%s/anno/%s.anno.vcf' %(sample, sample), 'w')
    for h in header: out.write(h)
    for t in texts: out.write(t)
    out.close()

if __name__ == "__main__":

    sample, vamos, filter = parse(sys.argv[1:])
    filterDepth, filterDist = int(filter[0]), float(filter[1])

    multi()


