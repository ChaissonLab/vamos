from typing import Any
import pyabpoa
import edlib

def shiftFrame(elements:list[str], use:str='cons', lexi:bool=True) -> \
        dict[str,tuple[int,str]]:
    """
    function to find number of bases to shift elements into the same frame

    Args:
        elements (list[str]): input sequences
        use (str, optional): which sequence should all sequences be shifted to
                             "cons": the consensus, "common":the most frequent,
                             "longest": the longest. Can also directly give the
                             target sequence for shift. Defaults to "cons".
        lexi (bool, optional): should sequences be sorted in lexicographical
                               order. Defaults to True.

    Returns:
        dict[str,tuple[int,str]]: final shifted sequence and number of bases
                                  (from head to tail) to be shifted for
                                  each sequence
    """

    elementSet = list(set(elements))
    if lexi: elementSet = sorted(elementSet)

    if use == 'cons':
        res = pyabpoa.msa_aligner().msa(sorted(elements), \
                out_cons=True, out_msa=False)
        consensus = res.cons_seq[0]
    elif use == 'common':
        consensus = max(sorted(elementSet), key = elements.count)
    elif use == 'longest':
        consensus = max(sorted(elementSet), key = len)
    else:
        consensus = use

    outDict = {}
    for seq in elementSet:
        seqShifts, dists = [], []
        # align outcome of each possible number of shifts to the consensus
        for shift in range(len(seq)):
            seqShift = seq[shift:] + seq[:shift]
            ed = edlib.align(consensus, seqShift, mode='NW')['editDistance']
            seqShifts.append(seqShift)
            dists.append(ed)
        # find the best number of shifts from all alignment results
        idx = dists.index(min(dists))
        outDict[seq] = ( seqShifts[idx], idx )

    return(consensus, outDict)

'''
print( shiftFrame(['AGTGGC', 'CAGTGG', 'AGTGGC', \
                   'AGTGGG', 'GGGAGT', 'GTGGCA']) )
'''


# define function to check for perfect cyclic shifts
def checkCyc(elements, lexi=True):

    if isinstance(elements, list):
        elementSet = list(set(elements))
    elif isinstance(elements, dict):
        elementSet = list(elements.keys())
    if lexi: elementSet = sorted(elementSet)

    checkDict = {}
    for m in elementSet:
        if len(m) not in checkDict: checkDict[len(m)] = []
        checkDict[len(m)].append(m)

    if isinstance(elements, list):
        new = []
        for m in elements:
            check = 0
            if m in new:
                new.append(m)
                continue

            for n in checkDict[len(m)]:
                if m != n and m in n+n:
                    if n in new:
                        new.append(n)
                    else:
                        new.append(m)
                    check = 1
                    break

            if check == 0: new.append(m)

    elif isinstance(elements, dict):
        new = {}
        for m in elementSet:
            check = 0

            for n in checkDict[len(m)]:
                if m != n and m in n+n:
                    if n in new:
                        new[n] += elements[m]
                    else:
                        new[m] = elements[m]
                    check = 1
                    break

            if check == 0: new[m] = elements[m]

    return(new)

