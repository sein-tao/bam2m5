#!/usr/bin/env python3
"""Convert bam file to pacbio m5 format
Usage: python3 <script> <in.bam> <ref.fa> <score_scheme> <out.m5>

in.bam should be sorted by coordinate for effiency
score_scheme: scoring parameter used for alignment, 
    in format: match,mismatch,gap_open,gap_extend
example: python3 bam2m5.py align.sorted.bam ref.fa -5,6,0,5 align.sorted.m5
"""

nargs = 4
from cbam2m5 import bam2m5, mapScore
from BioUtil import samFile, xzopen, cachedFasta

def main(inbam, fasta_file, score, out_file):
    fa = cachedFasta(fasta_file)
    with samFile(inbam, 'r') as sam, xzopen(out_file,'w') as out:
        ref_lengths = dict(zip(sam.references, sam.lengths))
        match, mismatch, gap_open, gap_extend = map(float, score.split(','))
        score_scheme = mapScore(mismatch, gap_open, gap_extend, match)
        for rec in sam:
            m5 = bam2m5(rec, fa, ref_lengths, score_scheme)
            if m5 is not None:
                print(*m5, file=out)
    fa.close()

if __name__ == '__main__':
    import sys
    if len(sys.argv) - 1 != nargs:
        sys.exit(__doc__)
    main(*sys.argv[1:])

