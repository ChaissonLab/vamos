import sys
import datetime
import logging

import bedlib

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)


Gymrek = 'Gymrek_etal_SupplementalData1_v2.bed'
GymrekLiftover = 'liftoverGymrekUCSC.output.bed'
vamos = '../s08_hotspotSegdup/all.mutRate.filterHotspot.tagSegdup.tsv'

GymrekDict = {}
with open(Gymrek) as f:
    for i,line in enumerate(f):
        if i == 0: continue
        fields = line.strip().split()
        chr19,start19,end19,rate = fields[:4]
        consensus = fields[13]
        GymrekDict[(f'chr{chr19}',start19,end19)] = [rate,consensus]

GymrekLiftoverDict = {}
with open(GymrekLiftover) as f:
    for line in f:
        chr,start,end,coor19,match = line.strip().split()
        chr19,temp = coor19.split(':')
        start19,end19 = temp.split('-')
        start19 = str(int(start19) - 1)
        GymrekLiftoverDict[(chr,start,end)] = GymrekDict[(chr19,start19,end19)]

vamosDict,tempDict = {},{}
with open(vamos) as f:
    for i,line in enumerate(f):
        if i == 0: continue
        fields = line.strip().split()
        chr,start,end,treeS,treeE,TR2Tree,nSNP,treeLen,winS,winE,period = fields[:11]
        rateV,rateNormBefore,rateNormAfter = fields[18:21]
        vamosDict[(chr,start,end)] = [rateV,period,rateNormBefore,rateNormAfter]
        tempDict[(chr,start,end)] = f'{chr}_{start}_{end}'

# intersect GymrekLiftoverDict by vamosDict
GymrekLiftoverDict = bedlib.bedIntersectTag(GymrekLiftoverDict,tempDict)

# check intersections of G by V and output results
out = open('intersectGymrekbyVamos.tsv','w')
out.write('chr\tstart\tend\tstartV\tendV\tperiod\tconsensus\trateG\trateV\trateNormBefore\trateNormAfter\n')
for (chr,start,end),(rate,consensus,tags) in GymrekLiftoverDict.items():
    if tags == 'NULL':
        outList1 = [chr,start,end,'nan','nan','nan',consensus]
        outList2 = [rate,'nan','nan','nan']
        out.write('\t'.join(outList1+outList2) + '\n')
    else:
        for tag in tags.split(','):
            chrV,startV,endV = tag.split('_')
            rateV,period,rateNormBefore,rateNormAfter = vamosDict[(chrV,startV,endV)]
            #if abs(int(startV) - int(start)) < int(period) and abs(int(endV) - int(end)) < int(period):
            outList1 = [chr,start,end,startV,endV,period,consensus]
            outList2 = [rate,rateV,rateNormBefore,rateNormAfter]
            out.write('\t'.join(outList1+outList2) + '\n')
out.close()

# intersect GymrekLiftoverDict by vamosDict
GymrekLiftoverDict = {(chr,start,end):f'{chr}_{start}_{end}' for chr,start,end in GymrekLiftoverDict.keys()}
vamosDict = bedlib.bedIntersectTag(vamosDict,GymrekLiftoverDict)

# check intersections of V by G and output results
out = open('intersectVamosbyGymrek.tsv','w')
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

