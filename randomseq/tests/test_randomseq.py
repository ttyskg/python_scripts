from main.randomseq import random_seq, calc_tm, select_seq
import re

def test_length_random_seq():
    seq = random_seq()
    seq10 = random_seq(slen=10)
    assert len(seq) == 20
    assert len(seq10) == 10


def test_content_random_seq():
    seq = random_seq()
    seq = ''.join(sorted(seq))
    assert re.fullmatch(r'A*?C*?G*?T*?', seq)


def test_calc_tm_short():
    seq = 'CTCTGCCTAGC'
    assert calc_tm(seq) == 36


def test_calc_tm_long():
    seq = 'CTCTATCTAGCTCTCT'
    assert calc_tm(seq) == 40.8


def test_select_seq():
    res = select_seq()
    assert type(res) == list
    assert len(res) == 3
    assert type(res[0][0]) == str
    assert type(res[1][1]) == float
    assert res[0][1] <= 62
    assert res[0][1] >= 59