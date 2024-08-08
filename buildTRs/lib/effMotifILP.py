#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
import logging
import numpy as np
from ortools.sat.python import cp_model
import Levenshtein as lv

from statistics import mean
import time 

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inFile', 'outFile']
optList = ['q', 'timeLimit', 'threads']
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0'
description = '''\nDescription:
This program applies the vamos ILP solution to find TR efficient motif sets for
given TR original motif sets. Searches passing the search time limit would most
likely return FEASIBLE (may not be OPTIMAL) solutions.
'''

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author, version, description]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput TR oriMotifs,  e.g. /in/oriMotifs.bed')
parser.add_argument(posList[1], type=str, \
    help='string\toutput TR effMotifs,  e.g. /out/effMotifs.bed')
# optional arguments
parser.add_argument('-q', '--'+optList[0], type=float, metavar='', default=0.1,
    help='string\tq compression parameter, e.g. 0.1, default 0.1')
parser.add_argument('-l', '--'+optList[1], type=int, metavar='', default=600,
    help='string\ttime limit (sec) for solution search, e.g. 600, default 600')
parser.add_argument('-t', '--'+optList[2], type=int, metavar='', default=2,
    help='string\tnumber of threads (2 recommended), e.g. 2, default 2')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
logging.info('Parsing Input Arguements...')
for key, value in argsDict.items():
    if key in posList: logging.info(f'Required Argument - {key}: {value}')
    if key in optList: logging.info(f'Optional Argument - {key}: {value}')
    vars()[key] = value # assign values of arguments into global variables
logging.info('Parsing Input Arguements Completed\n')

#--------------------------------------------------------
# global variables and user-defined functions
#--------------------------------------------------------

class MotifsInfo:
    def __init__(self):
        self.coords = []
        self.motifs = []
        self.counts = []
        self.countsTotals = [] # list of total motif counts of each TR
        self.editdists = []
        self.numOfTR = 0
        self.deltas = []
      
    """
    each row is a TR site
    """
    def readFile(self, inFile):
        with open(inFile, 'r') as f:
            for i,line in enumerate(f):
                chr,start,end,motifs,counts = line.strip().split('\t')[:5]
                self.coords.append((chr,start,end))
                self.motifs.append(motifs.split(','))
                self.counts.append([int(c) for c in counts.split(',')])
                self.countsTotals.append(sum(self.counts[-1]))
        self.numOfTR = len(self.motifs)
        
    """
    compute the pairwise Levenshtein edit distance between any pair of motifs
    """ 
    def pairwiseLevenshteinDist(self):
        for i in range(self.numOfTR):
            motifs, distDict = self.motifs[i], {}
            dist = [[0 for m1 in motifs] for m2 in motifs]
            for i1,m1 in enumerate(motifs):
                for i2,m2 in enumerate(motifs):
                    if i1 <= i2: continue
                    dist[i1][i2] = lv.distance(m1, m2)
                    dist[i2][i1] = dist[i1][i2]
            self.editdists.append(dist)
    
    """
    compute the delta based on the mean pairwise Levenshtein edit distance
    """
    def decideDelta(self):
        self.deltas = [0] * self.numOfTR
        for i in range(self.numOfTR):
            motifs = self.motifs[i]
            editdists = self.editdists[i]
            counts = self.counts[i]
            numOfmotifs = len(self.editdists[i])
            if len(motifs) == 1:
                self.deltas[i] = 0
                continue

            temp = []
            for i1,m1 in enumerate(motifs):
                for i2,m2 in enumerate(motifs):
                    if i1 < i2: # k replace h
                        temp.extend([editdists[i1][i2]] * counts[i1])
                    elif i2 < i1:
                        temp.extend([editdists[i2][i1]] * counts[i2])
            # percentile threshold of edit distance sum
            temp = np.array(temp)
            temp.sort()
            percentile = np.quantile(temp, q)
            # delta strategy 1
            # self.deltas[i] = int( percentile * self.countsTotals[i] )
            # delta strategy 2
            # self.deltas[i] = int( percentile * int(self.countsTotals[i]*q) )
            # delta strategy 3
            self.deltas[i] = int( percentile * int(self.countsTotals[i]*q) )
        return

"""
build up the ILP solver
"""
def ILPsolver(oriMotifs, oriCounts, costs, idxTR, delta):
    startTime = time.time()
    logging.info(f'Processing TR: {idxTR}')
    model = cp_model.CpModel()
    oriNum = len(oriMotifs)

    """
    create the variables
    """
    # the replacement indicator
    # x[i][j] == 1: motif i is replaced by motif j
    x = []
    for i in range(oriNum):
        t = [ model.NewBoolVar(f'x[{i},{j}]') for j in range(oriNum) ]
        x.append(t)

    # the sum replacement indicator for ILP form
    # y[j] == 1: motif j 
    y = [ model.NewBoolVar(f'y[{j}]') for j in range(oriNum) ]

    """
    create the constraints
    """
    # Each motif is replaced by exactly one motif (can be itself)
    # => all x replacement indicator sum to 1
    for i in range(oriNum):
        model.Add( sum(x[i][j] for j in range(oriNum)) == 1 )

    # The ILP upper bound of "1-sum(x[i][j] for i in range(oriNum)"
    for j in range(oriNum):
        model.Add( 1-sum(x[i][j] for i in range(oriNum)) <= 1-y[j] )

    # The ILP lower bound of "1-sum(x[i][j] for i in range(oriNum)"
    for j in range(oriNum):
        model.Add( 1-sum(x[i][j] for i in range(oriNum)) >= -oriNum*y[j]+1 )

    # Declare our intermediate boolean variable for channeling constraint
    # Channeling is usually implemented using half-reified linear constraints:
    # one constraint implies another (a => b), but not necessarily the other
    # way around (a <= b)
    b = []
    for i in range(oriNum):
        t = [ model.NewBoolVar(f'b[{i},{j}]') for j in range(oriNum) ]
        b.append(t)

    # motif i can be only replaced by motif j with higher or equal count
    # b[i][j]       => (counts[i] > counts[j])
    # not(b[i][j])  => (counts[i] <= counts[j])
    for i in range(oriNum):
        for j in range(oriNum):
            model.Add(oriCounts[i] > oriCounts[j]).OnlyEnforceIf(b[i][j])
            model.Add(oriCounts[i] <= oriCounts[j]).OnlyEnforceIf(b[i][j].Not())
    # (counts[i] > counts[j])   => (x[i][j] == 0)
    # (counts[i] <= counts[j])  => (x[i][j] >= 0)
    for i in range(oriNum):
        for j in range(oriNum):
            model.Add(x[i][j] == 0).OnlyEnforceIf(b[i][j])
            model.Add(x[i][j] >= 0).OnlyEnforceIf(b[i][j].Not())

    # The total cost is less than or equal to \delta
    model.Add(sum( oriCounts[i] * sum(x[i][j] * int(costs[i][j]) \
                    for j in range(oriNum)) \
                    for i in range(oriNum) ) <= delta)

    """
    create objective function
    minimizing the sum of efficient motif set size and replacement cost
    """
    # model.Minimize(sum(y))

    model.Minimize(delta * sum(y) + \
                   sum( oriCounts[i] * sum(x[i][j] * int(costs[i][j]) \
                        for j in range(oriNum)) \
                        for i in range(oriNum)))

    """
    invoke the solver
    """
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeLimit # search time limit
    solver.parameters.num_search_workers = threads # number of search threads
    # solver.parameters.enumerate_all_solutions = True

    status = solver.Solve(model)

    """
    logging the solution
    """
    targetMotifs = [] # targetMotif is replacement counterpart of each oriMotif
    effMotifs = []
    effCounts = [0] * oriNum
    mapcosts = []
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        if status == cp_model.OPTIMAL:
            statusTag = 'OPTIMAL'
            logging.info(f'OPTIMAL solution found')
        else:
            statusTag = 'FEASIBLE'
            logging.info(f'A solution found, but may not be optimal')

        for i in range(oriNum):
            check = False
            for j in range(oriNum):
                # motif i is replaced by motif j
                if solver.BooleanValue(x[i][j]):
                    targetMotifs.append(oriMotifs[j])
                    effCounts[j] += oriCounts[i] # add motif i counts to motif j
                    mapcosts.append(costs[i][j])
                    check = True
                    break
            assert check == True, f'motif i - {i} is mapped to nothing!'

        indicator = ['0' if c == 0 else '1' for c in effCounts]
        effMotifs = [ oriMotifs[i] for i,c in enumerate(effCounts) if c > 0 ]
        effCounts = [ c for c in effCounts if c > 0 ] # effCounts is updated!
        #effNum = int(solver.ObjectiveValue())
        effNum = len(effMotifs)
        #assert len(effMotifs) == effNum, f'unmatched solution set size!'
        endTime = time.time()

        # statistics
        logging.info(f'Time: {round(endTime-startTime,4)} sec')
        logging.info(f'Efficient/Original set size: {effNum}/{oriNum}\n')
        '''print('\nStatistics')
        print(f'  VNTR idx : {idxTR}')
        print(f'  status   : {solver.StatusName(status)}')
        print(f'  conflicts: {solver.NumConflicts()}')
        print(f'  branches : {solver.NumBranches()}')
        print(f'  wall time: {solver.WallTime()} s')'''

        return list(zip(oriMotifs, targetMotifs)), mapcosts, oriNum, effNum, \
            effMotifs, effCounts, indicator, statusTag

    # return the original motifs for TRs that have no feasible solutions
    else:
        statusTag = 'INFEASIBLE'
        logging.info(f'No solution found, returning the original motif set!')
        mapcosts = [ '0' for c in range(oriNum) ]
        indicator = [ '1' for cnt in oriCounts ]
        return list(zip(oriMotifs, oriMotifs)), mapcosts, oriNum, oriNum, \
            oriMotifs, oriCounts, indicator, statusTag

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    # read input
    logging.info(f'Reading input...')
    motifsInfo = MotifsInfo() 
    motifsInfo.readFile(inFile)
    motifsInfo.pairwiseLevenshteinDist()
    motifsInfo.decideDelta()
    logging.info(f'Reading input finish\n')

    out = open(outFile, 'w')
    for idxTR in range(motifsInfo.numOfTR):
        oriMotifs = motifsInfo.motifs[idxTR]
        oriCounts = motifsInfo.counts[idxTR]
        dists = motifsInfo.editdists[idxTR]
        delta = motifsInfo.deltas[idxTR]
        pairs,mapcosts,oriNum,effNum,effMotifs,effCounts,indicator,status = \
            ILPsolver(oriMotifs, oriCounts, dists, idxTR, delta)

        if len(pairs) > 0:
            chr,start,end = motifsInfo.coords[idxTR]
            out.write(chr +'\t'+ start +'\t'+ end +'\t'+ 
                ','.join(effMotifs) +'\t'+ 
                ','.join([str(c) for c in effCounts]) +'\t'+ 
                status +'\t'+ 
                f'{effNum}/{oriNum}' +'\t'+ 
                ','.join([str(pairs[i][0]) for i in range(len(pairs))]) +'\t'+ 
                ','.join([str(pairs[i][1]) for i in range(len(pairs))]) +'\t'+ 
                ','.join([str(cost) for cost in mapcosts]) +'\t'+ 
                ','.join(indicator) +'\t'+ 
                str(delta) +'\n')
    out.close()

logging.info('End of Program\n')

