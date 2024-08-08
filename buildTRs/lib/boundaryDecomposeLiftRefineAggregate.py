#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
import logging
import statistics
import seqShift
import seqDecompose
import pyabpoa

import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inFiles', 'outFile']
optList = ['outCluster', 'sampleCut']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
Refined TR entries from "boundaryDecomposeLiftRefine.py" of different samples
may not align perfectly by the coordinates and/or consensuses. The program
aggregates TR calls from different samples, grouping calls with any overlaps
and clustering by coordinates and consensus length using kmeans. Clusters are
regrouped by overlaps to give the final TR boundary. Note this program assumes
all input calls are on the same chromosome.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput files csv (sample,file),  e.g. /in/Files.csv')
parser.add_argument(posList[1], type=str, \
    help='string\toutput file,  e.g. /out/File')
# optional arguments
parser.add_argument('-c', '--'+optList[0], type=str, metavar='', \
    default=None, help='string\toutput cluster file,  default None')
parser.add_argument('-s', '--'+optList[1], type=int, metavar='', \
    default=3, help='integer\tmin num of samples to support a TR,  default 3')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
logging.info('Parsing Input Arguements...')
for key, value in argsDict.items():
    if key in posList: logging.info('Required Argument - %s: %s' %(key,value))
    if key in optList: logging.info('Optional Argument - %s: %s' %(key,value))
    vars()[key] = value # assign values of arguments into global variables
logging.info('Parsing Input Arguements Completed\n')


#--------------------------------------------------------
# global variables and user-defined functions
#--------------------------------------------------------

def overlap(s1,e1,s2,e2,base=1):
    if base==0: return max(0, min(e1, e2) - max(s1, s2))
    if base==1: return max(0, min(e1, e2) - max(s1, s2) + 1)

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    with open(inFiles) as f:
        files = { l.strip().split(',')[0]:l.strip().split(',')[1] for l in f }

    logging.info('Reading input TRs...')
    TRs = []
    for sample,file in files.items():
        # skip empty chromosome file (some sample may not have chrX/chrY)
        #if not os.path.exists(file): continue
        with open(file) as f:
            for line in f:
                fields = line.strip().split('\t')
                chr,start,end,_,_,_,_,_,_,consensus,_,motifs = fields
                TRs.append([sample,chr,int(start),int(end),consensus,motifs])

    TRs = sorted(TRs, key=lambda x: x[2]) # sort all calls by start coordinate
    logging.info('Reading input TRs finish')

    # check for empty input TR list
    if TRs == []:
        out = open(outFile, 'w')
        out.close()
        if outCluster:
            outC = open(outCluster, 'w')
            outC.close()
        logging.info('Empty input TR list!')
        sys.exit()

    # grouping different entries
    logging.info('Grouping TRs...')
    groupAll, group, checks = [], [TRs[0]], 0
    for i,(sample1,chr1,s1,e1,cons1,motifs1) in enumerate(TRs,start=1):

        # skip TRs[0] since it's already initialized
        if i == 1: continue

        for sample2,chr2,s2,e2,cons2,motifs2 in group:
            op = overlap(s1,e1,s2,e2)
            #if op >= 0.6*len1 or op >= 0.6*len2: checks += 1
            if op > 0: checks += 1

        #if checks == len(group):
        if checks > 0:
            group.append([sample1,chr1,s1,e1,cons1,motifs1])
            checks = 0
        else:
            groupAll.append(group)
            group, checks = [[sample1,chr1,s1,e1,cons1,motifs1]], 0

        # append the last group
        if i == len(TRs): groupAll.append(group)
    logging.info('Grouping TRs finish')

    # clustering of each group
    logging.info('Clustering TRs...')
    out = open(outFile+'.tmp', 'w')
    if outCluster: outC = open(outCluster, 'w')
    for group in groupAll:
        n = len(group)
        samples = ','.join([ g[0] for g in group ])
        chr = group[0][1]
        start = min([g[2] for g in group])
        end = max([g[3] for g in group])
        cons = ','.join([ g[4] for g in group ])

        # data for k means
        starts = [ g[2] for g in group ]
        ends = [ g[3] for g in group ]
        periods = [ len(g[4]) for g in group ]
        data = list(zip(starts, ends, periods))

        # scale data to -1 to 1 (invariant data is scaled to -1)
        scaler = MinMaxScaler(feature_range=(-1,1))
        scaler.fit(data) # fit the scaler
        data = scaler.transform(data) # transform the data

        # determine range of k for analysis
        kmin = 1
        kmax = 5*(n // len(files) + 1)
        if kmax > len(files):
            logging.info(f'Enormous group: {chr}:{start}_{end}, {kmax}')
            kmin = kmax - len(files)
        uniqueData = len(np.unique(data, axis=0))
        if uniqueData < kmax: kmax = uniqueData
        if len(data) == kmax: kmax = len(data) - 1 # special case: n = kmax

        # k means clustering
        inertias = []
        silhouettes, clusterCounts = [0], [str(n)] # special case: n = 1
        for i in range(kmin,kmax+1):
            kmeans = KMeans(n_clusters=i, n_init=10, random_state=123)
            kmeans.fit(data)
            inertias.append(round(kmeans.inertia_))
            if i > 1:
                silhouette = round(silhouette_score(data, kmeans.labels_),3)
                silhouettes.append(silhouette)
                temp = list(np.unique(kmeans.labels_, return_counts=True)[1])
                clusterCounts.append('-'.join(map(str, temp)))

        # picking the bestK (number of clusters) by the highest silhouette score
        bestK = silhouettes.index(max(silhouettes)) + kmin
        bestCounts = clusterCounts[bestK-kmin]
        silhouettes = ','.join(map(str,silhouettes))

        # calculate the consensus,start,end of each cluster under bestK
        kmeans = KMeans(n_clusters=bestK, n_init=10, random_state=123)
        kmeans.fit(data)
        clusters = { l:{'group':[]} for l in list(np.unique(kmeans.labels_)) }

        # save all raw group data
        for i,label in enumerate(kmeans.labels_):
            clusters[label]['group'].append(group[i])

        consensuses = []
        for label,cluster in clusters.items():
            temp = cluster['group']
            clusters[label]['n'] = len(temp)
            clusters[label]['start'] = round( np.mean([d[2] for d in temp]) )
            clusters[label]['end'] = round( np.mean([d[3] for d in temp]) )

            # shift all sequences and get consensus within each cluster
            seqs = [ d[4] for d in temp ]
            _,shiftDict = seqShift.shiftFrame(seqs,use='longest',lexi=False)
            shifted = [ shiftDict[s][0] for s in seqs ]
            clusters[label]['shifted'] = shifted
            res = pyabpoa.msa_aligner().msa(shifted,out_cons=True,out_msa=False)
            clusters[label]['consensus'] = res.cons_seq[0]
            consensuses += [res.cons_seq[0]] * len(cluster)

        # shift all cluster consensuses to match that of the bigest cluster
        _,shiftDict = seqShift.shiftFrame(consensuses,use='longest',lexi=False)
        for l,c in clusters.items():
            clusters[l]['shiftedCons'] = shiftDict[c['consensus']][0]
        #shifted = [ shiftDict[c['consensus']][0] for l,c in clusters.items() ]

        # sort all clusters by consensuses first and then coordinate for fix
        # output (cluster may assign cluster labels randomly)
        clusters = sorted(clusters.items(), key=lambda x: x[1]['shiftedCons'])
        clusters = sorted(clusters, key=lambda x: x[1]['n'])
        clusters = sorted(clusters, key=lambda x: x[1]['end'])
        clusters = sorted(clusters, key=lambda x: x[1]['start'])
        clusters = { j:v for j,(k,v) in enumerate(clusters) }

        # output
        if outCluster:
            # check if the optimal number of clusters is global optimum
            optimal = (bestK > kmin and bestK <= kmax) if bestK != 1 else True
            # output general info of each group
            temp = map(str, [chr,start,end,n,bestK,optimal,bestCounts,silhouettes,cons])
            outC.write('\t'.join(temp) +'\n')
            # output general info of each cluster
            for label,c in clusters.items():
                temp = [label,c['consensus'],c['start'],c['end'],c['n']]
                outC.write('\t'+ '\t'.join(map(str,temp)) +'\n')
                # output each sample under the cluster
                for j,d in enumerate(c['group']):
                    temp = [label,c['shifted'][j],d[4],d[2],d[3],d[0]]
                    outC.write('\t'+ '\t'.join(map(str,temp)) +'\n')

        for l,c in clusters.items():
            tag = f"{chr}:{start}-{end}"
            temp = [chr,c['start'],c['end'],tag,c['n'],l,c['shiftedCons'],c['consensus']]
            out.write('\t'.join(map(str,temp)) +'\n')

    out.close()
    if outCluster: outC.close()
    logging.info('Clustering TRs finish')

    # regrouping to get the final TR boundaries
    logging.info('Regrouping TRs...')
    TRs = []
    with open(outFile+'.tmp') as f:
        for line in f:
            fields = line.strip().split('\t')
            chr,start,end,coor,n,label,consensus,_ = fields
            start, end = int(start), int(end)
            # remove poorly supported call (support < sampleCut)
            if int(n) < sampleCut: continue
            # remove small mono/di-nucleotide calls inside bigger calls
            coor = coor.split(':')[1]
            s,e = coor.split('-')
            if len(consensus) < 3 and end-start < 0.5*(int(e)-int(s)): continue

            TRs.append([chr,start,end,n,consensus])

    TRs = sorted(TRs, key=lambda x: x[1]) # sort all calls by start coordinate

    # check for empty TR list after refinement
    if TRs == []:
        out = open(outFile, 'w')
        out.close()
        logging.info('Error: empty TR list after refinement!')
        sys.exit()

    groupAll, group, checks = [], [TRs[0]], 0
    for i,(chr1,s1,e1,n1,cons1) in enumerate(TRs,start=1):

        # skip TRs[0] since it's already initialized
        if i == 1: continue

        for chr2,s2,e2,n2,cons2 in group:
            op = overlap(s1,e1,s2,e2)
            #if op >= 0.6*len1 or op >= 0.6*len2: checks += 1
            if op > 0: checks += 1

        #if checks == len(group):
        if checks > 0:
            group.append([chr1,s1,e1,n1,cons1])
            checks = 0
        else:
            groupAll.append(group)
            group, checks = [[chr1,s1,e1,n1,cons1]], 0

        # append the last group
        if i == len(TRs): groupAll.append(group)
    logging.info('Regrouping TRs finish')

    out = open(outFile, 'w')
    for group in groupAll:
        chr = group[0][0]
        start = min([g[1] for g in group])
        end = max([g[2] for g in group])
        ns = [int(g[3]) for g in group]
        n = sum(ns)
        seqs = [ g[4] for g in group ]
        '''res = pyabpoa.msa_aligner().msa(seqs,out_cons=True,out_msa=False)
        cons = res.cons_seq[0]'''
        cons = seqs[ns.index(max(ns))]
        seqs = ','.join(seqs)
        ns = ','.join(map(str,ns))

        out.write('\t'.join(map(str,[chr,start,end,n,cons,ns,seqs])) +'\n')
    out.close()


logging.info('End of Program\n')

