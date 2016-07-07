#!/usr/bin/env python3
import sys
from BioUtil import xzopen
assert(len(sys.argv[1:]) == 2)
afile, bfile = sys.argv[1:]
fields = ['info', 'qseq', 'match_pattern', 'ref_seq']

with open(afile) as afh, open(bfile) as bfh:
    line = 0
    diff_line = 0
    while True:
        if diff_line > 5:
            break
        try:
            aline = next(afh)
            bline = next(bfh)
            line += 1
        except StopIteration:
            break
        a_elts = aline.rstrip().rsplit(None, 3)
        b_elts = bline.rstrip().rsplit(None,3)
        if aline != bline:
            print("======line %d differ====" % line)
            diff_line += 1
        for i in range(len(fields)):
            if a_elts[i] != b_elts[i]:
                print("a", fields[i], ":", a_elts[i])
                print("b", fields[i], ":", b_elts[i])

