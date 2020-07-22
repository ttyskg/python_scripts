import random

# content of test_sample.py

def random_seq(len=20):
    return ''.join(random.choices('ATGC', k=len))


def calc_tm(seq):
    '''calculate Tm value of inputed sequence

    :parm seq: DNA sequence for calucating its Tm.
    :type seq: string
    
    :rtype: numeric
    :return: Tm value of the DNA sequence
    '''

    seq = seq.upper()
    nA = seq.count('A')
    nT = seq.count('T')
    nC = seq.count('C')
    nG = seq.count('G')

    if len(seq) < 14:
        return 2 * (nA + nT) + 4 * (nG + nC)
    else:
        return round(64.9 + 41 * (nG + nC - 16.4) / (len(seq)), 1)

