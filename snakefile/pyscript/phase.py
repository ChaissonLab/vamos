# -*- coding: UTF-8 -*-

import os, sys, getopt
from multiprocessing import Pool

def parse(argv):

    try:
        opts, args = getopt.getopt(argv, 'hb:w:l:r:s:t:o:', ['inBam=', 'winDir=', 'longshot=', 'ref=', 'samtools=', 'thread=', 'outDir='])
    except getopt.GetoptError as err:
        print(err)
        sys.stdout.write('\nphase.py -b <inBam> -w <winDir> -l <longshot> -r <ref> -s <samtools> -t <thread> -o <outDir>\n'); sys.stdout.flush()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            sys.stdout.write('\nphase.py -b <inBam> -w <winDir> -l <longshot> -r <ref> -s <samtools> -t <thread> -o <outDir>\n'); sys.stdout.flush()
            sys.exit()

    thread = 1

    for opt, arg in opts:
        if opt in ('-b', '--inBam'):
            inBam = arg
        elif opt in ('-w', '--winDir'):
            winDir = arg
        elif opt in ('-l', '--longshot'):
            longshot = arg
        elif opt in ('-r', '--ref'):
            ref = arg
        elif opt in ('-s', '--samtools'):
            samtools = arg
        elif opt in ('-t', '--thread'):
            thread = int(arg)
        elif opt in ('-o', '--outDir'):
            outDir = arg

    return(inBam, winDir, longshot, ref, samtools, thread, outDir)


def phase(key):

    chr,start,end = key
    os.makedirs(outDir +'/%s_%s-%s' %(chr, start, end))
    prefix = outDir +'/%s_%s-%s/%s_%s-%s' %(chr, start, end, chr, start, end)
    extractBam = prefix +'.extract.bam'
    phaseVCF = prefix +'.phase.vcf'
    phaseBam = prefix +'.phase.bam'
    phaseBam1 = prefix +'.phase.h1.bam'
    phaseBam2 = prefix +'.phase.h2.bam'

    os.system('%s view -h -S -b %s "%s:%s-%s" > %s' %(samtools, inBam, chr, start, end, extractBam))
    os.system('%s index %s' %(samtools, extractBam))
    os.system('%s --bam %s --ref %s --out %s -O %s' %(longshot, extractBam, ref, phaseVCF, phaseBam))
    os.system('%s index %s' %(samtools, phaseBam))

    '''os.system('%s view -H %s > %s.h.tmp' %(samtools, phaseBam, phaseBam))
    os.system('%s view %s | grep "HP:i:1" > %s.sam.tmp' %(samtools, phaseBam, phaseBam))
    os.system('cat %s.h.tmp %s.sam.tmp | %s view -bS > %s' %(phaseBam, phaseBam, samtools, phaseBam1))
    os.system('%s index %s' %(samtools, phaseBam1))
    os.system('%s view %s | grep "HP:i:2" > %s.sam.tmp' %(samtools, phaseBam, phaseBam))
    os.system('cat %s.h.tmp %s.sam.tmp | %s view -bS > %s' %(phaseBam, phaseBam, samtools, phaseBam2))
    os.system('%s index %s' %(samtools, phaseBam2))
    os.system('rm %s* %s.h.tmp %s.sam.tmp' %(extractBam, phaseBam, phaseBam))'''


def multi(inBam, winDir, longshot, ref, samtools, thread, outDir):

    if not os.path.exists(outDir): os.makedirs(outDir)

    files, keys = [f for f in os.listdir(winDir) if 'vntrs.' in f], []

    for file in files:

        key = os.path.basename(file).split('.')[1]
        chr,coor = key.split('_')

        keys.append( [chr] + coor.split('-') )

    with Pool(thread) as p:
        p.map(phase, keys)


if __name__ == "__main__":

    inBam, winDir, longshot, ref, samtools, thread, outDir = parse(sys.argv[1:])

    multi(inBam, winDir, longshot, ref, samtools, thread, outDir)
    

