from main.randomseq import random_seq
import re

def test_length_random_seq():
    seq = random_seq()
    seq10 = random_seq(len=10)
    assert len(seq) == 20
    assert len(seq10) == 10


def test_content_random_seq():
    seq = random_seq()
    seq = ''.join(sorted(seq))
    assert re.fullmatch(r'A*?C*?G*?T*?', seq)
