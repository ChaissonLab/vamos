# -*- coding: UTF-8 -*-

import os, sys, getopt
import pyabpoa as pa
from multiprocessing import Pool

def parse(argv):

    try:
        opts, args = getopt.getopt(argv, 'hs:v:t:f:', ['sample=', 'vamos=', 'thread=', 'filter='])
    except getopt.GetoptError as err:
        print(err)
        sys.stdout.write('\nanno_cons_first.py -s <sample> -v <vamos> -f <filter>\n'); sys.stdout.flush()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            sys.stdout.write('\nanno_cons_first.py -s <sample> -v <vamos> -f <filter>\n'); sys.stdout.flush()
            sys.exit()

    thread = 1

    for opt, arg in opts:
        if opt in ('-s', '--sample'):
            sample = arg
        elif opt in ('-v', '--vamos'):
            vamos = arg
        elif opt in ('-t', '--thread'):
            thread = int(arg)
        elif opt in ('-f', '--filter'):
            minDepth = int(arg)

    return(sample, vamos, thread, minDepth)


def readvcf(vcf):
    with open(vcf) as f:
        for line in f:
            if line.startswith('#'): continue

            chr, pos, id, ref, alt, qual, filter, info, format, genotype = line.strip().split('\t')
            end, ru, svtype, altanno, len = info.split(';')[:-1]
            annoDict[(chr,pos,end)] = ['\t'.join([chr, pos, id, ref, alt, qual, filter]), ';'.join([end, ru, svtype]), ';'.join([altanno, len]), ';'.join([altanno, len])]

'''
def medaka(fa, seqs):

    out = open(fa+'_seqs', 'w')
    for i,s in enumerate(seqs): out.write('>%s\n%s\n' %(str(i),s))
    out.close()
    dir = os.path.dirname(fa)

    if os.path.exists(fa+'.fai'): os.remove(fa+'.fai')
    if os.path.exists(fa+'.map-ont.mmi'): os.remove(fa+'.map-ont.mmi')
    if os.path.exists(fa+'.polished.fa.gaps_in_draft_coords.bed'): os.remove(fa+'.polished.fa.gaps_in_draft_coords.bed')
    if os.path.exists(fa+'bam.bai'): os.remove(fa+'bam.bai')
    #os.system('mini_align -i %s -r %s -m -p %s' %(bam.replace('bam','fasta'), fa, bam))
    os.system('medaka_consensus -f -x -i %s_seqs -d %s -o %s -m r941_prom_sup_g507' %(fa, fa, dir))
    os.system('mv %s/consensus.fasta %s' %(dir, fa))
'''

def anno(key):

    if not os.path.exists('analysis/%s/anno/%s' %(sample, key)):
        os.makedirs('analysis/%s/anno/%s' %(sample, key))

    vntr = 'analysis/split/vntrs.%s.bed' %(key)
    bam = 'analysis/%s/phase/%s/%s.phase.bam' %(sample, key, key)
    anno = 'analysis/%s/anno/%s/anno.bed' %(sample, key)

    fa = 'analysis/%s/anno/%s/tmp.anno.fa' %(sample, key)
    bed = 'analysis/%s/anno/%s/tmp.anno.bed' %(sample, key)
    vcf = 'analysis/%s/anno/%s/tmp.anno.vcf' %(sample, key)

    # vamos
    if not os.path.exists(bam): return
    os.system('%s --readwise -b %s -r %s -s %s -o %s' %(vamos, bam, vntr, sample, anno))

    annoDict = {}

    with open(anno) as file:

        for line in file:

            if line.startswith('#'): continue

            try:
                chr,start,end,motifs,fields = line.strip().split('\t')
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
            if len(seqs1) < minDepth or len(seqs2) < minDepth:
                sys.stdout.write('Warnings: skipping locus of poor coverage, %s:%s-%s\n' %(chr,start,end)); sys.stdout.flush()
                continue

            if sum([len(s) for s in seqs1]) > 1000000 or sum([len(s) for s in seqs2]) > 1000000:
                sys.stdout.write('Warnings: skipping sequences unhandable by MSA, %s:%s-%s\n' %(chr,start,end)); sys.stdout.flush()
                continue

            cons1 = pa.msa_aligner().msa(seqs1, out_cons=True, out_msa=False).cons_seq[0]
            cons2 = pa.msa_aligner().msa(seqs2, out_cons=True, out_msa=False).cons_seq[0]

            out = open(bed, 'w')
            out.write('\t'.join([chr,start,end,motifs]) +'\n')
            out.close()

            # annotation of the consensus sequence
            out = open(fa, 'w')
            out.write('>cons\n' + cons1 +'\n')
            out.close()

            '''# polish consensus for ont data
            medaka(fa, seqs1)'''

            os.system('%s --single_seq -b %s -r %s -s %s -o %s' %(vamos, fa, bed, sample, vcf))

            with open(vcf) as f:
                for l in f:
                    if l.startswith('#'): continue

                    chr, pos, id, ref, alt, qual, filter, info, format, genotype = l.strip().split('\t')
                    end, ru, svtype, altanno_h1, len_h1 = info.split(';')[:-1]
                    annoDict[(chr,pos,end)] = ['\t'.join([chr, pos, id, ref, alt, qual, filter]), ';'.join([end, ru, svtype]), ';'.join([altanno_h1, len_h1]), ';'.join([altanno_h1, len_h1])]

            out = open(fa, 'w')
            out.write('>cons\n' + cons2 +'\n')
            out.close()

            '''# polish consensus for ont data
            medaka(fa, seqs2)'''

            os.system('%s --single_seq -b %s -r %s -s %s -o %s' %(vamos, fa, bed, sample, vcf))

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
    for f in [fa, bed]: os.remove(f)


def multi():

    keys = [f.split('.')[1] for f in os.listdir('analysis/split') if 'vntrs.chr1_10627-10009901' in f]
    if not os.path.exists('analysis/%s/anno' %(sample)):
        os.makedirs('analysis/%s/anno' %(sample))

    with Pool(thread) as p:
        p.map(anno, keys)

    texts, header = [], []
    for key in keys:
        v = 'analysis/%s/anno/%s/tmp.anno.vcf' %(sample, key)

        if not os.path.exists(v): continue

        with open(v) as f: header = [ l for l in f if l.startswith('#') ]

        with open(v) as f:
            for line in f:
                if not line.startswith('#'): texts.append(line)

    out = open('analysis/%s/anno/%s.anno.vcf' %(sample, sample), 'w')
    for h in header: out.write(h)
    for t in texts: out.write(t)
    out.close()

if __name__ == "__main__":

    sample, vamos, thread, minDepth = parse(sys.argv[1:])

    multi()


