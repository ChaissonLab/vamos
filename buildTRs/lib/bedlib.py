from typing import Any

#@title function to check for number of overlaps for adjacent bed regions
def checkOverlap(bedFile:str, closeEnd:bool=False) -> int:
    """
    function to check for number of overlaps for adjacent bed regions

    Args:
        bedFile (str): input bed file, bed regions must be sorted
        closeEnd (bool, optional): should first and last bases be counted as
                                   overlapping. Defaults to False.

    Returns:
        int: number of overlaps of adjacent bed regions
    """

    counter = 0
    with open('bedFile') as f:
        for i,line in enumerate(f,start=1):

            chr,start,end = line.strip().split()[:3]

            if i != 1:
                if chr == chrL:
                    if int(start) <= int(endL): counter += 1

            chrL,startL,endL = chr,start,end

    return counter


#@title function to tag features in A by all the overlapping features in B
def bedIntersectTag(dictA:dict[tuple[str,str,str],list[Any]], \
                    dictB:dict[tuple[str,str,str],str]) \
                    -> dict[tuple[str,str,str],list[Any]]:
#def bedIntersectTag(dictA:dict, dictB:dict) -> dict:
    """
    function to tag features in A by all the overlapping features in B. Linear
    time complexity is achieved by double while loop. This requires features
    in set B to be non-overlapping

    Args:
        dictA (dict[tuple[str,str,str],list[Any]]):
            features of setA, (chr,start,end) as keys, values as list
        dictB (dict[tuple[str,str,str],str]):
            features of setB, (chr,start,end) as keys, values as str

    Returns:
        dict[tuple[str,str,str],list[Any]]:
            sorted features of setA, (chr,start,end) as keys, values as list
            where last item is the overlapping info from setB
    """

    chrs = []
    for key in list(dictA.keys()):
        if key[0] not in chrs: chrs.append(key[0])

    for chr in chrs:
        overlaps, overlapsKey, i, j = [], [], 0, 0
        keysA = [ k for k in dictA.keys() if k[0]==chr ]
        keysA.sort(key = lambda k: int(k[1]))
        keysB = [ k for k in dictB.keys() if k[0]==chr ]
        keysB.sort(key = lambda k: int(k[1]))

        # handle chrs that are not in B at all
        if len(keysB) == 0:
            for k in keysA: dictA[k].append('NULL')

        while i < len(keysA) and j < len(keysB):

            chrA, startA, endA = keysA[i]
            chrB, startB, endB = keysB[j]

            # B not reaches A
            if int(endB) < int(startA):
                # no more next B, conclude this A and go to next A
                # also back to first B intersected with last A to re-check
                if j == len(keysB) - 1:
                    temp = ','.join(overlaps) if overlaps else 'NULL'
                    dictA[keysA[i]].append(temp)
                    i += 1
                    if overlapsKey: j = overlapsKey[0]
                    overlaps, overlapsKey = [], []
                # there are more next B, go to next B
                else:
                    j += 1
            # B passes A, conclude this A and go to next A
            # also back to first B intersected with last A to re-check
            elif int(endA) < int(startB):
                temp = ','.join(overlaps) if overlaps else 'NULL'
                dictA[keysA[i]].append(temp)
                i += 1
                if overlapsKey: j = overlapsKey[0]
                overlaps, overlapsKey = [], []
            # an intersect, record this B
            else:
                overlaps.append(dictB[keysB[j]])
                overlapsKey.append(j)
                # no more next B, conclude this A and go to next A
                # also back to first B intersected with last A to re-check
                if j == len(keysB) - 1:
                    temp = ','.join(overlaps) if overlaps else 'NULL'
                    dictA[keysA[i]].append(temp)
                    i += 1
                    if overlapsKey: j = overlapsKey[0]
                    overlaps, overlapsKey = [], []
                # there are more next B, go to next B
                else:
                    j += 1

    dictA = {k:v for k,v in sorted(dictA.items(), key=lambda x: int(x[0][1]))}
    # final sort of results
    outDict = {}
    for chr in chrs:
        keysA = [ k for k in dictA.keys() if k[0]==chr ]
        keysA.sort(key = lambda k: int(k[1]))
        for k in keysA: outDict[k] = dictA[k]

    return outDict

'''
bedIntersect({('chr1','1','5'):[], ('chr2','11','15'):[]}, \
             {('chr1','1','5'):'chr1_1_5', ('chr2','11','15'):'chr2_11_15'})
'''


#@title function to remove features in A if it overlaps any features in B
def bedIntersectRemove(dictA:dict[tuple[str,str,str],list[Any]], \
                       dictB:dict[tuple[str,str,str],Any]) \
                       -> dict[tuple[str,str,str],list[Any]]:
#def bedIntersectRemove(dictA:dict, dictB:dict) -> dict:
    """
    function to remove any features in A if it overlaps any features in B.
    Linear time complexity is achieved by double while loop. This requires
    features in each set to be non-overlapping

    Args:
        dictA (dict[tuple[str,str,str],list[Any]]):
            features of setA, (chr,start,end) as keys, values as list
        dictB (dict[tuple[str,str,str],Any]):
            features of setB, (chr,start,end) as keys, values as any

    Returns:
        dict[tuple[str,str,str],list[Any]]:
            sorted features of setA, (chr,start,end) as keys, values as list
            where last item is the overlapping info from setB

    """

    chrs, newA = [], {}
    for key in list(dictA.keys()):
        if key[0] not in chrs: chrs.append(key[0])

    for chr in chrs:
        i, j = 0, 0
        keysA = [ k for k in dictA.keys() if k[0]==chr ]
        keysA.sort(key = lambda k: int(k[1]))
        keysB = [ k for k in dictB.keys() if k[0]==chr ]
        keysB.sort(key = lambda k: int(k[1]))

        while i < len(keysA) and j < len(keysB):

            chrA, startA, endA = keysA[i]
            chrB, startB, endB = keysB[j]

            # B not reaches A
            if int(endB) < int(startA):
                # no more next B, record this A and go to next A
                if j == len(keysB) - 1:
                    newA[keysA[i]] = dictA[keysA[i]]
                    i += 1
                # there are more next B, go to next B
                else:
                    j += 1
            # B passes A, record this A and go to next A
            elif int(endA) < int(startB):
                newA[keysA[i]] = dictA[keysA[i]]
                i += 1
            # an intersect, skip this A and go to next A
            else:
                i += 1

    return newA

'''
bedIntersectRemove({('chr1','1','5'):[], ('chr2','11','15'):[]}, \
                   {('chr1','1','5'):'', ('chr2','16','18'):''})
'''

