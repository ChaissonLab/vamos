import sys
import os
import pysam
import logging

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)

bam = sys.argv[1]
outFile = sys.argv[2]

allDict = {}

samfile = pysam.AlignmentFile(bam, "rb")

# Iterate through reads and get their aligned ending points
for read in samfile.fetch():
    # Check if the read is mapped (not unmapped)
    if not (read.is_unmapped or read.is_secondary or read.is_supplementary):
    #if not (read.is_unmapped or read.is_secondary):
        supp = read.is_supplementary
        readID = read.query_name
        readRefChrom = read.reference_name
        readRefStart = read.reference_start
        readRefEnd = read.reference_end
        if readID not in allDict:
             allDict[readID] = []
        allDict[readID].append((readRefChrom,readRefStart,readRefEnd,supp))
# Close the BAM file
samfile.close()

out = open(outFile, 'w')
for seq,aligns in allDict.items():
    coor,coorNew,coorExpand = seq.split('--')
    chr,start,end = coor.split('__')
    chrNew,startNew,endNew = coorNew.split('__')
    chrExpand,startExpand,endExpand = coorExpand.split('__')
    lengthNew = int(endNew) - int(startNew)
    startDiff = int(startNew) - int(startExpand)
    endDiff = int(endExpand) - int(endNew)
    contigs = [ a[0] for a in aligns ]
    starts = [ int(a[1]) for a in aligns ]
    ends = [ int(a[2]) for a in aligns ]
    supps = [ int(a[3]) for a in aligns ]
    if len(set(contigs)) != 1:
        logging.info('Warning! Divergent alignment on different contigs!')
        tag = 'divergent'
        print(seq)
        print(aligns)
        # filter alignment not on the primary contig
        primaryIndex = [ j for j,s in enumerate(supps) if not s ][0]
        primaryContig = contigs[primaryIndex]
        pickIndex = [ j for j,c in enumerate(contigs) if c == primaryContig ]
        contigs = [ contigs[j] for j in pickIndex ]
        starts = [ starts[j] for j in pickIndex ]
        ends = [ ends[j] for j in pickIndex ]
        supps = [ supps[j] for j in pickIndex ]
        contig = primaryContig
        start = min(starts)
        end = max(ends)
    else:
        if len(aligns) != 1:
            logging.info('Warning! Multiple alignment on the same contig!')
            tag = 'multiple'
        else:
            tag = 'single'
        contig = contigs[0]
        start = min(starts)
        end = max(ends)
    diff = end-start - lengthNew
    prop = diff / lengthNew
    out.write('\t'.join(map(str, [coor,coorNew,coorExpand,tag,contig,start,end,lengthNew,end-start,diff,prop])) +'\n')

out.close()


