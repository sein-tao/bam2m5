#!/usr/bin/env python3
"""Convert bam file to pacbio m5 format
Usage: python3 <script> <in.bam> <ref.fa> [--score <score_scheme>] <out.m5>

in.bam should be sorted by coordinate for effiency
score_scheme: scoring parameter used for alignment, 
    in format: match,mismatch,gap_open,gap_extend
example: python3 bam2m5.py align.sorted.bam ref.fa -5,6,0,5 align.sorted.m5
"""

from cbam2m5 import bam2m5, mapScore
from BioUtil import samFile, xzopen, cachedFasta
import argparse

def main():
    parser = argparse.ArgumentParser(description = "Convert bam file to pacbio m5 format")
    parser.add_argument("inbam", metavar = "in.bam", help = "input bam file")
    parser.add_argument("fasta", metavar = "ref.fa", help = "reference file")
    parser.add_argument("--score", metavar = "match,mismatch,gap_open,gap_extend", 
            default="-5,6,0,5", help="scoring parameter for alignment, default:-5,6,0,5. see README for details")
    parser.add_argument("outm5", metavar="out.m5", help = "output m5 file")
    args = parser.parse_args()

    fa = cachedFasta(args.fasta)
    match, mismatch, gap_open, gap_extend = map(float, args.score.split(','))
    score_scheme = mapScore(mismatch, gap_open, gap_extend, match)
    with samFile(args.inbam, 'r') as sam, xzopen(args.outm5,'w') as out:
        ref_lengths = dict(zip(sam.references, sam.lengths))
        for rec in sam:
            m5 = bam2m5(rec, fa, ref_lengths, score_scheme)
            if m5 is not None:
                print(*m5, file=out)
    fa.close()

if __name__ == '__main__':
    main()

