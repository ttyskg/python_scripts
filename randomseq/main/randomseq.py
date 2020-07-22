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


print(''.join(sorted(random_seq())))random_se