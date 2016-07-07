#!/usr/bin/env python
# enum is about 20x slower than direct access in python3.4
# it is not used for performance reason
# import enum
import re
import itertools
import operator
class CIGAR:
    MATCH   = 0
    INS     = 1
    DEL     = 2
    REF_SKIP    = 3
    SOFT_CLIP   = 4
    HARD_CLIP   = 5
    PAD     = 6
    EQUAL   = 7
    DIFF    = 8

    _names = ['MATCH', 'INS', 'DEL', 'REF_SKIP', 
            'SOFT_CLIP', 'HARD_CLIP', 'PAD', 'EQUAL', 'DIFF']
    _chars = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
    _char2int = dict(zip(_chars, range(len(_chars)) ))
    _cigar_pattern = re.compile(r"(?P<length>\d+)(?P<operation>[%s])" % ''.join(_chars))
    op = operator.itemgetter(0)
    length = operator.itemgetter(1)

    @classmethod
    def get_name(cls, num):
        return cls._names[num]
    @classmethod
    def get_char(cls, num):
        return cls._chars[num]

    @classmethod
    def str2tuple(cls, cigarstring):
        result = list()
        for m in re.finditer(cls._cigar_pattern, cigarstring):
            if m:
                result.append((
                    cls._char2int[m.group('operation')], 
                    int(m.group('length'))
                    ))
        return tuple(result)

class cigartuple(tuple):
    @classmethod
    def from_str(cls, cigarstring):
        return cls(CIGAR.str2tuple(cigarstring))

    def infer_read_length(self):
        if hasattr(self, '_read_length'):
            return self._read_length
        total = 0
        for op, length in self:
            if op in set((CIGAR.MATCH, CIGAR.INS, CIGAR.SOFT_CLIP, CIGAR.EQUAL, CIGAR.DIFF)):
                total += length
        self._read_length = total
        return total
    @staticmethod
    def _is_clip(item):
        return item[0] in (CIGAR.SOFT_CLIP, CIGAR.HARD_CLIP)
    @property
    def query_start(self):
        left_trim = sum(length for op, length in 
                itertools.takewhile(self._is_clip, self) if op == CIGAR.SOFT_CLIP)
        return left_trim

    @property
    def query_end(self):
        right_trim = sum(length for op, length in 
                itertools.takewhile(self._is_clip, reversed(self)) if op == CIGAR.SOFT_CLIP)
        return self.infer_read_length() - right_trim

    @property
    def left_hard_clip(self):
        if CIGAR.op(self[0]) == CIGAR.HARD_CLIP:
            return CIGAR.length(self[0])
        else:
            return 0
    @property
    def right_hard_clip(self):
        if CIGAR.op(self[-1]) == CIGAR.HARD_CLIP:
            return CIGAR.length(self[-1])
        else:
            return 0

