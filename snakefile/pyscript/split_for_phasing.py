# -*- coding: UTF-8 -*-

import os, sys, getopt, datetime

def parse(argv):

    try:
        opts, args = getopt.getopt(argv, 'hr:w:o:', ['inVNTRs=', 'winSize=', 'outDir='])
    except getopt.GetoptError as err:
        print(err)
        sys.stdout.write('\nsplit.py -r <inVNTRs> -w <winSize> -o <outDir>\n'); sys.stdout.flush()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            sys.stdout.write('\nsplit.py -r <inVNTRs>  -w <winSize> -o <outDir>\n'); sys.stdout.flush()
            sys.exit()

    for opt, arg in opts:
        if opt in ('-r', '--inVNTRs'):
            inVNTRs = arg
        elif opt in ('-w', '--winSize'):
            winSize = int(arg)
        elif opt in ('-o', '--outDir'):
            outDir = arg

    return(inVNTRs, winSize, outDir)


def split(inVNTRs, winSize, outDir):

    vntrsDict, motifsDict, vntrsBedDict, moitfsCSVDict = {}, {}, {}, {}

    # read vntrs and motifs
    with open(inVNTRs) as f:
        for line in f:

            chr, start, end, motifs = line.strip().split('\t')
            if chr not in vntrsDict: vntrsDict[chr] = []

            vntrsDict[chr].append([int(start), int(end)])
            motifsDict[chr+'_'+start+'-'+end] = motifs

    for chr,coor in vntrsDict.items(): vntrsDict[chr] = sorted(coor, key=lambda x: x[0])

    # set up genomic windows
    if not os.path.exists(outDir): os.makedirs(outDir)

    for chr,coor in vntrsDict.items():

        head, tail, bedO, csvO = coor[0][0], 1, [], []

        for start,end in coor:
            if start > head + winSize:

                tag = '%s_%s-%s' %(chr,str(head),str(tail))
                bed = open(outDir + '/vntrs.'+ tag +'.bed', 'w')
                for o in bedO: bed.write(o +'\n')
                bed.close()

                head, bedO, csvO = start, [], []

            tail = end
            bedO.append('\t'.join([chr,str(start),str(end),motifsDict[chr+'_'+str(start)+'-'+str(end)]]))

        tag = '%s_%s-%s' %(chr,str(head),str(tail))
        bed = open(outDir + '/vntrs.'+ tag +'.bed', 'w')
        for o in bedO: bed.write(o +'\n')
        bed.close()

if __name__ == "__main__":

    inVNTRs, winSize, outDir = parse(sys.argv[1:])

    split(inVNTRs, winSize, outDir)
