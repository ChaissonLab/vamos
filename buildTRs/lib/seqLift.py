import sys
import re
from typing import Any

def liftOnePass(cigar:str, startAln:int, refAsKey:bool) -> dict[int,list[int]]:
    """
    function to map all query positions to ref positions (ref/query as key)
    (does not handle case when startRef and endRef are in different aligned
    segment e.g., startRef in primary but endRef in supplementary)
    clipped bases are ignored in the output

    Args:
        cigar (str): cigar string of the query alignment to the reference
        startAln (int): start position of the query alignment on the reference
        refAsKey (bool): if output uses ref positions as key (otherwise query)

    Returns:
        dict[int,list[int]]: mapped ref-to-query positions, unique ref/query
            positions as key
    """

    # convert the cigar string
    cigar = re.findall(r'(\d+)(\w)', cigar.replace('=','M'))
    cigar = [ [int(c[0]),c[1]] for c in cigar ]

    mapDict, pointerRef, pointerQuery = {}, startAln+1, 1

    for c in cigar:

        # block of clips
        if c[1] in ['H', 'S']:
            # no need to move the pointerRef since it's clip
            pointerQuery += c[0]

        # block of matches and mismatches
        elif c[1] in ['M', 'X']:
            for j in range(c[0]):
                if refAsKey:
                    if pointerRef not in mapDict: mapDict[pointerRef] = []
                    mapDict[pointerRef].append(pointerQuery)
                else:
                    if pointerQuery not in mapDict: mapDict[pointerQuery] = []
                    mapDict[pointerQuery].append(pointerRef)
                pointerRef += 1
                pointerQuery += 1

        # block of deletions in query
        elif c[1] == 'D':
            # no need to move the pointerQuery since it's query delete
            for j in range(c[0]):
                if refAsKey:
                    if pointerRef not in mapDict: mapDict[pointerRef] = []
                    mapDict[pointerRef].append(pointerQuery)
                else:
                    if pointerQuery not in mapDict: mapDict[pointerQuery] = []
                    mapDict[pointerQuery].append(pointerRef)
                pointerRef += 1

        # block of insertions in query
        elif c[1] == 'I':
            # no need to move the pointerRef since it's ref delete
            for j in range(c[0]):
                if refAsKey:
                    if pointerRef not in mapDict: mapDict[pointerRef] = []
                    mapDict[pointerRef].append(pointerQuery)
                else:
                    if pointerQuery not in mapDict: mapDict[pointerQuery] = []
                    mapDict[pointerQuery].append(pointerRef)
                pointerQuery += 1

        else:
            sys.exit('Error: unknown cigar notation')

    return mapDict

'''print(liftOnePass(cigar='5M', startAln=5, refAsKey=True))'''


def liftQuery2Ref(query:str, cigar1Pass:dict[int,list[int]], reverse:bool, \
                  startQuery:int, endQuery:int) -> tuple[int,int]:
    """
    function to map given query segment to the reference start/end positions
    (does not handle case when startQuery and endQuery are in different aligned
    segment e.g., startQuery in primary but endQuery in supplementary)

    Args:
        query (str): full query sequence (including hard clips)
        cigar1Pass (dict[int,list[int]]): parsed cigar (query positions as key)
        reverse (bool): True if the query aligns to the reverse strand
        startQuery (int): start position of the query segment (natural)
        endQuery (int): end position of the query segment (natural)

    Returns:
        tuple[int,int]: mapped reference start/end positions (forward)
    """

    # reverse startQuery/endQuery if on reverse strand
    # For reverse aligned query, bam cigar is the forward cigar
    # So, startQuery/endQuery should be reversed to match the cigar
    if reverse:
        total = len(query)
        temp = total - endQuery + 1
        endQuery = total - startQuery + 1
        startQuery = temp

    # get the lifted sequence from parsed cigar (skip half-clipped case)
    if startQuery in cigar1Pass and endQuery in cigar1Pass:
        startRef = min(cigar1Pass[int(startQuery)])
        endRef = max(cigar1Pass[int(endQuery)])
    else:
        startRef, endRef = 0, 0

    return (startRef, endRef)

'''print(liftQuery2Ref(query='AAGTC', reverse=False, startQuery=2, endQuery=4, \
        cigar1Pass={1:[6], 2:[7], 3:[8], 4:[9], 5:[10]}))'''


def liftRef2Query(query:str, cigar1Pass:dict[int,list[int]], reverse:bool, \
                  startRef:int, endRef:int) -> tuple[int,int]:
    """
    function to map given ref segment to the query start/end positions
    (does not handle case when startRef and endRef are in different aligned
    segment e.g., startRef in primary but endRef in supplementary)

    Args:
        query (str): full query sequence (including hard clips)
        cigar1Pass (dict[int,list[int]]): parsed cigar (ref positions as key)
        reverse (bool): True if the query aligns to the reverse strand
        startRef (int): start position of the reference segment (forward)
        endRef (int): end position of the reference segment (forward)

    Returns:
        tuple[int,int]: mapped query start/end positions (natural)
    """

    # get the lifted sequence from parsed cigar (skip half-clipped case)
    if startRef in cigar1Pass and endRef in cigar1Pass:
        startQuery = min(cigar1Pass[int(startRef)])
        endQuery = max(cigar1Pass[int(endRef)])
    else:
        startQuery, endQuery = 0, 0

    # reverse startQuery/endQuery if on reverse strand
    # For reverse aligned query, bam cigar is the forward cigar
    # So, startQuery/endQuery should be reversed to give the natural start/end
    if reverse:
        total = len(query)
        temp = total - endQuery + 1
        endQuery = total - startQuery + 1
        startQuery = temp

    return (startQuery, endQuery)

'''print(liftRef2Query(query='AAGTC', reverse=False, startRef=7, endRef=9, \
        cigar1Pass={6:[1], 7:[2], 8:[3], 9:[4], 10:[5]}))'''

