from typing import Any

#@title function for the string decomposing algorithm
def stringDecomposer(vntr:str, motifs:list[str], \
                     match:int=0, mismatch:int=1, indel:int=1) -> list[str]:
    """
    function for the string decomposing algorithm

    Args:
        vntr (str): vntr string to be decomposed
        motifs (list[str]): motif set used for vntr decomposition
        match (int, optional): score for match. Defaults to 0.
        mismatch (int, optional): score for mismatch. Defaults to 1.
        indel (int, optional): score for indel. Defaults to 1.

    Returns:
        list[str]: decomposed vntr
    """

    optMotifs = [] # final optimal motifs
    vntrJunctions = [] # final optimal junction points on vntr
    N = len(vntr)    
    right = 0 # del-right: deletion from vntr to make the motif
    down = 1 # ins-down: insertion into vntr to make the motif
    diagonal = 2 # diag: match or mismatch

    # each motif has one DP sheet (motif indexed by row, vntr indexed by col)
    # the DP book is to glue all DP sheets of individual motifs by first row
    scoreBook, traceBook = [], []
    for m in motifs:
        scoreBook.append([ [0]*(N+1) for i in range(len(m)) ])
        traceBook.append([ [0]*(N+1) for i in range(len(m)) ])

    # the general DP book first row is handled separately
    # (traceRow0 records wrap around motifs not directions)
    scoreRow0 = [ 0 for n in range(N+1) ]
    traceRow0 = [ 0 for n in range(N+1) ]

    # for each DP sheet, initialize first column as gap
    for m,motif in enumerate(motifs):
        for i in range(len(motif)):
            scoreBook[m][i][0] = - i * indel
            traceBook[m][i][0] = 1

    '''
    dynamic programming propogation
    propogate col after col instead of row after row (hence j on outer loop)
    propogate from top-left to bottom-right
    first row of each DP sheet is propagated from the general DP book first row
    first row of the general DP book is handled seperately (next loop block)
    '''
    for j in range(N): # j indexes vntr position (current j+1)
        for m,motif in enumerate(motifs): # m indexes motif id
            for i in range(len(motif)): # i indexes motif position

                # handle first row
                if i == 0:
                    ms = 0 if motif[0] == vntr[j] else -mismatch
                    diagScore = scoreRow0[j] + ms
                    insScore = scoreRow0[j] - indel
                # handle other rows
                else:
                    ms = 0 if motif[i] == vntr[j] else -mismatch
                    diagScore = scoreBook[m][i-1][j] + ms
                    insScore = scoreBook[m][i-1][j+1] - indel
                delScore = scoreBook[m][i][j] - indel

                # determine the optimal score and path (for other rows)
                if diagScore >= insScore and diagScore >= delScore:
                    optScore=diagScore
                    optPath=diagonal
                elif insScore >= diagScore and insScore >= delScore:
                    optScore=insScore
                    optPath=down
                else:
                    optScore=delScore
                    optPath=right

                scoreBook[m][i][j+1] = optScore
                traceBook[m][i][j+1] = optPath

        '''
        transition at the first row (switch between motifs)
        '''
        maxScore = scoreBook[0][len(motifs[0])-1][j+1]
        maxIndex = 0
        for m,motif in enumerate(motifs): # m indexes motif id
            if scoreBook[m][len(motifs[m])-1][j+1] > maxScore:
                maxScore = scoreBook[m][len(motifs[m])-1][j+1]
                maxIndex = m
        scoreRow0[j+1] = maxScore
        traceRow0[j+1] = maxIndex

    '''
    trace back
    om: index to find the optimal motif
    oi: index on motif to find the optimal ending position of the optimal motif
    oj: index on vntr to find the optimal junction point for each optimal motif
    '''
    # find the last optimal motif
    for m,motif in enumerate(motifs): # m indexes motif id
        if m==0 or scoreBook[m][len(motif)-1][N] > maxScore:
            maxIndex = m
            maxScore = scoreBook[m][len(motif)-1][N]

    # trace back
    om, oi, oj = maxIndex, len(motifs[maxIndex])-1, N
    optMotifs.append(om)
    vntrJunctions.append(oj)
    while oj > 0: # trace back starts from oj=N and ends at oj=0

        direction = traceBook[om][oi][oj]

        if oi == 0: # a wrap around may occur at first row
            if direction == right: # still at current motif
                oj -= 1
            # wrap around to another motif
            # (traceRow0 records wrap around motifs)
            else:
                om = traceRow0[oj-1]
                oi = len(motifs[om]) - 1
                if oj > 1:
                    optMotifs.append(om)
                    vntrJunctions.append(oj-1)
                oj -= 1
        else: # normal trace back at other rows
            if direction == diagonal:
                oj -= 1
                oi -= 1
            elif direction == down:
                oi -= 1
            else:
                oj -= 1

    # reverse the traceback order
    optMotifs.reverse()
    vntrJunctions.reverse()
    traceback = []
    for idx,m in enumerate(optMotifs):
        if idx == 0:
            traceback.append([m,0,vntrJunctions[idx]])
        else:
            traceback.append([m,vntrJunctions[idx-1],vntrJunctions[idx]])

    decomposition = []
    for t in traceback: decomposition.append(vntr[t[1]:t[2]])

    return decomposition

#stringDecomposer('ACTACTAGTAT',['AC','ACT'])


# define function to check for perfect concatenation
def checkConcat(elements, lexi=True):

    if isinstance(elements, list):
        elementSet = list(set(elements))
    elif isinstance(elements, dict):
        elementSet = list(elements.keys())
    if lexi: elementSet = sorted(elementSet)

    # find all spliting cases
    splitDict = {}
    for e1 in elementSet:

        checkSet = [e for e in elementSet if e != e1]

        # vector for junction points at each iteration
        # (0 position is always a junction point)
        junctions = [1 if i == 0 else 0 for i in range(len(e1)+1)]
        # vector for splitting traceback
        splits = [[] for i in range(len(e1)+1)]

        proceed = True
        while proceed:
            proceed = False
            update = [0 for i in range(len(e1)+1)]
            # try to extend each existing junction point to the next one
            for pos,status in enumerate(junctions):
                if status == 0 or pos == len(junctions): continue

                e3 = e1[pos:] # sequence for concatenation checking
                for e2 in checkSet:
                    if e3.startswith(e2): # head concatenation found
                        proceed = True
                        update[pos+len(e2)] = 1 # add new junction point
                        # update splitting traceback, keep more splits for tie
                        if len(splits[pos])+1 > len(splits[pos+len(e2)]):
                            splits[pos+len(e2)] = splits[pos] + [e2]

            if proceed: junctions = update

        splitDict[e1] = splits[-1]

    # create output
    if isinstance(elements, dict):
        new = dict(elements)
        for e,v in splitDict.items():
            if v == []: continue
            del new[e]
            for s in v: new[s] += elements[e]

    return(new)

