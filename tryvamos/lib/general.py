# -*- coding: UTF-8 -*-
from typing import Any

def parseBed(bedFile):

    bedDict = {}
    with open(bedFile) as f:
        for line in f:
            chr,start,end = line.strip().split('\t')[:3]
            if chr not in bedDict: bedDict[chr] = {}
            bedDict[chr][(start,end)] = 0

    return bedDict


def alignGlobal(seqA:str, seqB:str, distance:bool=False, \
                match:int=1, mismatch:int=-2, indel:int=-2) -> int:
    """
    function for global alignment (Wunch-Needleman) of seqA and seqB

    Args:
        seqA (str): sequence A for alignment
        seqB (str): sequence B for alignment
        distance (bool, optional): use distance configuration instead of
                                   scoring configuration. Defaults to False.
        match (int, optional): score for match. Defaults to 1.
        mismatch (int, optional): score for mismatch. Defaults to -2.
        indel (int, optional): score for indel. Defaults to -2.

    Returns:
        int: score of the final global alignment
    """

    m, n, update, results, scoreTable = len(seqA), len(seqB), [0,0,0], 0, []

    # initialize the DP table (n+1 rows, m+1 cols, propagate from top-left)
    for i in range(n+1):
        if i == 0: # first row as dels from seqA
            scoreTable.append( [indel*j for j in range(m+1)] )
        else: # first col as ins into seqA
            scoreTable.append( [indel*i if j == 0 else 0 for j in range(m+1)] )

    # initialize the traceback table (n+1 rows, m+1 cols)
    # traceTable = [ [[0,0] for j in range(m+1)] for i in range(n+1) ]

    for i in range(1, n+1): # n rows to propagate
        for j in range(1, m+1): # m cols to propagate

            # compute all of the recurrence scores
            update[0] = scoreTable[i][j-1] + indel # left-to-right
            update[1] = scoreTable[i-1][j-1] + \
                    [mismatch, match][ int(seqA[j-1]==seqB[i-1]) ] # diagonal
            update[2] = scoreTable[i-1][j] + indel # upper-to-lower

            # update the DP table
            scoreTable[i][j] = min(update) if distance else max(update)

    results = scoreTable[n][m]

    return results


