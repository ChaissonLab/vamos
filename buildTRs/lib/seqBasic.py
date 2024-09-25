from typing import Any

def reComDNA(seq:str) -> str:
    """function to get revserse complement

    Args:
        seq (str): input DNA sequence

    Returns:
        str: reverse complement of the input DNA sequence
    """

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rc = ''.join(complement.get(base, base) for base in reversed(seq))

    return(rc)

#print( alignGlobal('ACCGGT', 'ACCGCT') )

