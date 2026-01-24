import sys
import datetime
import logging

import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

chrList = [f'chr{i}' for i in range(1, 23)]

Evan = 'Evan.hg38.bed'
EvanLiftover = 'liftoverEvanUCSC.output.bed'
vamos = '../s08_hotspotSegdup/all.mutRate.filterHotspot.tagSegdup.tsv'

EvanLiftoverDict, tempDict = {}, {}
with open(Evan) as f:
    for i,line in enumerate(f):
        #if i == 0: continue
        fields = line.strip().split()
        chr,start,end,motifs = fields[:4]
        EvanLiftoverDict[(chr,start,end)] = [motifs]


'''EvanDict = {}
with open(Evan) as f:
    for i,line in enumerate(f):
        #if i == 0: continue
        fields = line.strip().split()
        chrT2T,startT2T,endT2T = fields[:3]
        EvanDict[(chrT2T,startT2T,endT2T)] = ['']

# homogenize EvanListover results of multiple mapping regions
EvanLiftoverDict, tempDict = {}, {}
with open(EvanLiftover) as f:
    for line in f:
        chr,start,end,coorT2T,match = line.strip().split()
        chrT2T,temp = coorT2T.split(':')
        startT2T,endT2T = temp.split('-')
        startT2T = str(int(startT2T) - 1)
        if chr not in chrList: continue
        if (chrT2T,startT2T,endT2T) not in tempDict:
            tempDict[(chrT2T,startT2T,endT2T)] = []
        tempDict[(chrT2T,startT2T,endT2T)] += [(chr,start,end,match)]
for (chrT2T,startT2T,endT2T),maps in tempDict.items():
    starts = [int(start) for chr,start,end,match in maps]
    ends = [int(end) for chr,start,end,match in maps]
    start = str(min(starts))
    end = str(max(ends))
    EvanLiftoverDict[(chrT2T,startT2T,endT2T)] = EvanDict[(chrT2T,startT2T,endT2T)]'''

vamosDict,tempDict = {},{}
with open(vamos) as f:
    for i,line in enumerate(f):
        if i == 0: continue
        fields = line.strip().split()
        chr,start,end,treeS,treeE,TR2Tree,nSNP,treeLen,winS,winE,period = fields[:11]
        rateV,rateNormBefore,rateNormAfter = fields[18:21]
        vamosDict[(chr,start,end)] = [rateV,period,rateNormBefore,rateNormAfter]
        tempDict[(chr,start,end)] = f'{chr}_{start}_{end}'

# intersect EvanLiftoverDict by vamosDict
EvanLiftoverDict = bedlib.bedIntersectTag(EvanLiftoverDict,tempDict)

# check intersections of G by V and output results
out = open('intersectEvanbyVamos.tsv','w')
out.write('chr\tstart\tend\tstartV\tendV\tperiod\trateV\trateNormBefore\trateNormAfter\n')
for (chr,start,end),(motifs,tags) in EvanLiftoverDict.items():
    if tags == 'NULL':
        outList1 = [chr,start,end,'nan','nan','nan']
        outList2 = ['nan','nan','nan']
        out.write('\t'.join(outList1+outList2) + '\n')
    else:
        selectedTags = []
        for tag in tags.split(','):
            chrV,startV,endV = tag.split('_')
            rateV,period,rateNormBefore,rateNormAfter = vamosDict[(chrV,startV,endV)]
            selectedTags.append([rateV,period,rateNormBefore,rateNormAfter])
        # select the tag with the highest rateV
        selectedTags.sort(key=lambda x: x[0], reverse=True)
        rateV,period,rateNormBefore,rateNormAfter = selectedTags[0]
        #if abs(int(startV) - int(start)) < int(period) and abs(int(endV) - int(end)) < int(period):
        outList1 = [chr,start,end,startV,endV,period]
        outList2 = [rateV,rateNormBefore,rateNormAfter]
        out.write('\t'.join(outList1+outList2) + '\n')
out.close()

# intersect EvanLiftoverDict by vamosDict
EvanLiftoverDict = {(chr,start,end):f'{chr}_{start}_{end}' for chr,start,end in EvanLiftoverDict.keys()}
vamosDict = bedlib.bedIntersectTag(vamosDict,EvanLiftoverDict)

# check intersections of V by G and output results
out = open('intersectVamosbyEvan.tsv','w')
out.write('chr\tstart\tend\tperiod\ttype\tlength\ttag\n')
for (chr,start,end),(rateV,period,rateNormBefore,rateNormAfter,tags) in vamosDict.items():
    ty = 'STR'
    if int(period) > 6: ty = 'VNTR'
    size = int(end) - int(start)
    length = 'short'
    if size >= 20: length = 'medium'
    if size >= 50: length = 'long'
    if size >= 300: length = 'ultra'
    if tags == 'NULL':
        outList = [chr,start,end,period,ty,length,'nan']
        out.write('\t'.join(outList) + '\n')
    else:
        outList = [chr,start,end,period,ty,length,tags]
        out.write('\t'.join(outList) + '\n')
out.close()

logging.info('End of Program\n')
