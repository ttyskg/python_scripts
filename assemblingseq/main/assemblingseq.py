#!/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math
import random
import subprocess
import sys


def make_assembling_seqs(n=5, slen=20):
    """Create DNA sequences for assembling, such as for Gibson assembly.
    The default parameters are:
        Second structure thleshold (deltaG): -4
        Self dimers threshold: 4
        Melting temperature: 58 - 62 Celsius temperature
    
    :parm n: Number os sequences to create
    :type n: int
    :param slen: Sequence length
    :type slen: int

    :rtype: set
    :return: Set of strings which are for assembly
    """
    res = set()
    while len(res) < n:
        seq = random_seq(slen)
        if not check_secondary_structure(seq):
            continue
        if not check_self_dimers(seq):
            continue
        if not check_melting_temperature(seq):
            continue
        res.add(seq)
    return res


def random_seq(slen):
    return ''.join(random.choices("ATGC", k=slen))


def check_secondary_structure(seq, cutoff=-4):
    score = calc_secondary_structure(seq)
    return score > cutoff


def calc_secondary_structure(seq):
    cmd = "hybrid-ss-min --NA=DNA -q " + str(seq)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    p.wait()
    return float(p.stdout.read())


def check_self_dimers(seq, cutoff=4):
    score = primer_dimers(seq, seq)
    return score < cutoff


def primer_dimers(primer1, primer2):
    score = []
    for i in range(len(primer1)):
        score.append(complement_compare(primer1[len(primer1)-i-1:len(primer1)], primer2[len(primer1)-i-1:len(primer1)]))
        return max(score)


def complement_compare(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError('primer sequences in complement_compare are not equql')
    score = 0
    for i in range(len(seq1)):
        if i < 10:
            score += basecompare(seq1[i], seq2[(i*-1)-1])*2
        else:
            score += basecompare(seq1[i], seq2[((i*-1)-1)])
    return score


def basecompare(base1, base2):
    if base1 == 'G' and base2 == 'C':
        return 1
    elif base1 == 'C' and base2 == 'G':
        return 1
    elif base1 == 'A' and base2 == 'T':
        return 1
    elif base1 == 'T' and base2 == 'A':
        return 1
    else:
        return -1


def check_melting_temperature(seq, min_threshold=58, max_threshold=62):
    tm = oligoTm(seq)
    return min_threshold < tm < max_threshold


def oligoTm(seqobj):
    """
    Takes either a SeqRecord object, a Seq object, or a string
    and computes the melting temp based on the NN model (yes?).
    This is Kun's code
    
    CHECK THE NN PARAMETERS
    From Uri Laserson
    """
    
    if isinstance(seqobj,SeqRecord):
        seq = str(seqobj.seq).upper()
    elif isinstance(seqobj,Seq):
        seq = str(seqobj).upper()
    elif isinstance(seqobj,str):
        seq = seqobj.upper()
    
    # set the default tm parameters
    C_primer = 200 # 200nM in standard pcr
    C_Mg = 1.5 #Taq Buffer i
    C_MonovalentIon = 70 #20mM Tris-Cl + 50mM KCL in Taq Buffer  
    C_dNTP = 0.25 #mM
    percentage_DMSO = 0
    percentage_annealed = 50
    
    percentage_annealed = percentage_annealed/100.0
    percentage_DMSO = percentage_DMSO/100.0
    
    # Some constants
    R = 1.987
    deltaH =  {"AA": -7.6,  "TT": -7.6, "AT": -7.2, "TA": -7.2, "CA": -8.5, "TG": -8.5, "GT": -8.4, "AC": -8.4,"CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,"CG": -10.6,"GC": -9.8, "GG": -8.0, "CC": -8.0, "A": 2.2, "T": 2.2, "G": 0.0, "C": 0.0}
    deltaS = {"AA": -21.3, "TT": -21.3, "AT": -20.4, "TA": -21.3, "CA": -22.7, "TG": -22.7, "GT": -22.4, "AC": -22.4, "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,"CG": -27.2, "GC": -24.4, "GG": -19.9, "CC":-19.9, "A": 6.9, "T": 6.9, "G": 0.0, "C": 0.0}
    
    C_SodiumEquivalent = C_MonovalentIon + 120 * math.sqrt(C_Mg-C_dNTP)
    seqLength = len(seq)
    dH = 0.2 + deltaH[str(seq[0])] + deltaH[str(seq[len(seq)-1])]
    dS = -5.7 + deltaS[seq[0]] + deltaS[seq[len(seq)-1]]
    for i in range(0, seqLength - 1):
        dH += deltaH[str(seq[i:i+2])]
        dS +=  deltaS[seq[i:i+2]]
    dS = dS + 0.368 * seqLength * math.log(C_SodiumEquivalent/1000.0)
    # val = math.log(C_primer*(1-percentage_annealed)/percentage_annealed)
    Tm = (dH * 1000) / (dS + R * (math.log(C_primer*(1-percentage_annealed)/percentage_annealed)-21.4164)) - 273.15 - 0.75*percentage_DMSO
    return Tm


if __name__ == "__main__":
    args = sys.argv
    n = int(args[1])
    slen = int(args[2])

    print("Inputs:\n",\
          "    Number of seq: {}\n".format(n),\
          "    Sequence length: {}\n".format(slen))

    res = make_assembling_seqs(n, slen)
    with open("./assembling_seq.txt", mode="w") as f:
        for seq in res:
            f.write(seq + '\n')
