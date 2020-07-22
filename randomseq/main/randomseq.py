#!/usr/bin/python3

import csv
import random
import sys

# content of test_sample.py

def random_seq(slen=20):
    return ''.join(random.choices('ATGC', k=slen))


def calc_tm(seq):
    '''calculate Tm value of inputed sequence

    :parm seq: DNA sequence for calucating its Tm.
    :type seq: string
    
    :rtype: float, r
    :return: Tm value of the DNA sequence (round at the first decimal place)
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


def calc_gc(seq):
    '''calculate GC content of inputed sequence.

    :parm seq: DNA sequence for calucating its Tm.
    :type seq: string
    
    :rtype: float, r
    :return: GC% of the DNA sequence (round at the first decimal place)
    '''

    seq = seq.upper()
    nA = seq.count('A')
    nT = seq.count('T')
    nC = seq.count('C')
    nG = seq.count('G')

    return round((nC + nG) / len(seq) * 100, 1)


def select_seq(n=3, slen=20, minTm=59, maxTm=62):
    '''Create DNA sequences that meet the condition.

    :parm n: Number of sequences to create
    :type n: int
    :parm slen: Sequence length
    :parm slen: int
    :parm minTm: minimum Tm value
    :type minTm: float
    :parm maxTm: maximum Tm value
    :type maxTm: float

    :rtype: set
    :return: Set of tuples which contained (seq, Tm)
    '''

    res = set()
    while len(res) < n:
        seq = random_seq()
        tm = calc_tm(seq)
        if minTm <= tm <= maxTm:
            res.add((seq, tm))
    
    return list(res)


if __name__ == '__main__':
    args = sys.argv
    n = int(args[1])
    slen = str(args[2])
    minTm = float(args[3])
    maxTm = float(args[4])

    print('Inputs:\n',\
          '  Number of seq: {}\n'.format(n),\
          '  Sequence length: {}\n'.format(slen),\
          '  minimum Tm: {}\n'.format(minTm),\
          '  maximum Tm: {}\n'.format(maxTm))
    
    res = select_seq(n=n, slen=slen, minTm=minTm, maxTm=maxTm)
    with open('./selected_seq.csv', mode='w') as f:
        writer = csv.writer(f)
        writer.writerow(['sequences', 'Tm'])

        for (seq, Tm) in res:
            writer.writerow([seq, Tm])

    print('Finished!')