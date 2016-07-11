#!/usr/bin/env python3
# cython: profile=True
"""Convert bam file to pacbio m5 format
Usage: python3 <script> <in.bam> <ref.fa> <score_scheme> <out.m5>

in.bam should be sorted by coordinate for effiency
score_scheme: scoring parameter used for alignment, 
    in format: match,mismatch,gap_open,gap_extend
example: python3 bam2m5.py align.sorted.bam ref.fa -5,6,0,5 align.sorted.m5
"""

from __future__ import print_function
nargs = 4
from cigar import CIGAR, cigartuple
from BioUtil import samFile, xzopen, cachedFasta

from pysam.calignmentfile cimport AlignmentFile, AlignedSegment
# from libcpp.vector cimport vector
import warnings

class mapScore:
    def __init__(self, mismatch, gap_open, gap_extend, match = 1):
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_extend

    def score(self, cigar_len, cigar_count):
        score = 0
        score += self.match * cigar_len[CIGAR.EQUAL]
        score += self.mismatch * cigar_len[CIGAR.DIFF]
        score += self.gap_open * (cigar_count[CIGAR.INS] + cigar_count[CIGAR.DEL]) 
        score += self.gap_extend * (cigar_len[CIGAR.INS] + cigar_len[CIGAR.DEL])
        return int(score)

# m5 fileds
# qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq
def bam2m5(AlignedSegment rec, object fa, dict ref_lengths, object score_scheme):
    cdef str qseq, rseq
    if rec.is_unmapped:
        return None
    if rec.query_sequence is None :
        warnings.warn("%s -> %s is not parsed due to no query sequece" 
                %(rec.query_name, rec.reference_name) )
        return None
    cigar = cigartuple(rec.cigartuples)
    qseq = rec.query_sequence[rec.query_alignment_start:rec.query_alignment_end]
    rseq = fa[rec.reference_name][rec.reference_start:rec.reference_end]
    qseq, rseq = padding_seq(qseq, rseq, rec)

    if rec.is_reverse:
        qseq = reverse_complement(qseq)
        rseq = reverse_complement(rseq)
    mp = match_pattern(qseq, rseq)

    cigar_len, cigar_count = rec.get_cigar_stats()
    if cigar_len[CIGAR.EQUAL] == 0 and cigar_len[CIGAR.MATCH] > 0:
        cigar_len[CIGAR.DIFF] = mp.count('*') - (cigar_len[CIGAR.INS] + cigar_len[CIGAR.DEL])
        cigar_len[CIGAR.EQUAL] = cigar_len[CIGAR.MATCH] - cigar_len[CIGAR.DIFF]

    qstart = rec.query_alignment_start + cigar.left_hard_clip
    qend = rec.query_alignment_end + cigar.left_hard_clip
    qlen = rec.query_length + cigar_len[CIGAR.HARD_CLIP]
    if rec.is_reverse:
        qstart, qend = qlen - qend, qlen - qstart

    # qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq
    return (rec.query_name, qlen, qstart, qend, '+', '', # add a space to match m5 formating
            rec.reference_name, ref_lengths[rec.reference_name],
            rec.reference_start, rec.reference_end,
            '-' if rec.is_reverse else '+',
            score_scheme.score(cigar_len, cigar_count),
            *(cigar_len[op] for op in (CIGAR.EQUAL, CIGAR.DIFF, CIGAR.INS, CIGAR.DEL)), 
            # cigar_len[CIGAR.EQUAL], cigar_len[CIGAR.DIFF],
            # cigar_len[CIGAR.INS], cigar_len[CIGAR.DEL],
            rec.mapping_quality,
            qseq, mp, rseq,
            )


cdef padding_seq(str qseq, str rseq, AlignedSegment rec):
    cdef size_t qpos = 0
    cdef size_t rpos = 0
    cdef list qgap, rgap
    qgap, rgap = list(), list() # gap positions and length
    cdef size_t op, length
    for op, length in rec.cigartuples:
        if op in (CIGAR.EQUAL, CIGAR.DIFF, CIGAR.MATCH):
            qpos += length
            rpos += length
        elif op == CIGAR.INS:
            rgap.append((rpos, length))
            qpos += length
        elif  op == CIGAR.DEL:
            qgap.append((qpos, length))
            rpos += length
        elif op in (CIGAR.SOFT_CLIP, CIGAR.HARD_CLIP):
            continue
        else:
            warnings.warn("%s %s %s" %
                    (rec.query_name, CIGAR.get_name(op)))
    qseq = add_gap(qseq, qgap)
    rseq = add_gap(rseq, rgap)
    if (len(qseq) != len(rseq)):
        warnings.warn("sequence length not match, maybe error: %s -> %s"
                % (rec.query_name, rec.reference_name))
    return qseq, rseq

cdef str add_gap(str seq, list gaps):
    cdef size_t pos = 0
    cdef list frag = list()
    cdef size_t start, length
    for start, length in gaps:
        frag.append(seq[pos:start])
        frag.append('-' * length)
        pos = start
    else:
        frag.append(seq[pos:])
    return ''.join(frag)


cdef dict complement = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
cdef str reverse_complement(str seq):
    return "".join(complement.get(base, base) for base in reversed(seq) )

cdef str match_pattern(str qseq, str rseq):
    return "".join('|' if q == r else '*' for q, r in zip(qseq.upper(), rseq.upper()))



