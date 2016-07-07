#!/usr/bin/env python
import enum
from collections import namedtuple
import re
from cigar import CIGAR

EQUAL = CIGAR.EQUAL
DIFF = CIGAR.DIFF
DEL = CIGAR.DEL

class State:
    MATCH = CIGAR.EQUAL
    MUT = CIGAR.DIFF
    DEL = CIGAR.DEL

MDop = namedtuple('MDop', ['type', 'length', 'seq'])
def parse_tag(tag):
    "MD string to MD tuple"
    #MD tag format: [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
    groups = re.split("(\d+)", tag)
    result = list()
    # ['', 'num', 'ATCG', 'num', 'ATCG', ..., 'num', '']
    for i, mark in enumerate(groups):
        if i % 2 == 0:
            if mark == '':
                continue
            elif mark[0] == '^':
                seq = mark[1:]
                elt = MDop(State.DEL, len(seq), seq)
            else:
                elt = MDop(State.MUT, len(mark), mark)
        else:
            length = int(mark)
            if length == 0:
                continue
            else:
                elt = MDop(State.MATCH, length, None)
        result.append(elt)
    return result
