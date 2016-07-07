m5 fileds
===========

- qName: query name
- qLength: the total length of query sequence, incluing clipped parts
- qStart: alignment start coordinate of query, see note 
- qEnd: alignment end coordnate of query, see note 
- qStrand: always '+'
- tName : reference name
- tLength: total length of reference
- tStart: alignment start coordinate on reference
- tEnd: alignment end coordinate on reference
- tStrand: alignment strand, '+' or '-'
- score: mapping score, calculated from mapping scheme
- numMatch: number of bases of equal match (exclude mismatch)
- numMismatch: number of bases of mismatch
- numIns: number of bases of insertion
- numDel: .. of deletion
- mapQV: mapping quality value
- qAlignedSeq: query sequece aligned, '-' for gap base
- matchPattern: relationship of each position in qAlignedSeq and tAlignedSeq, '|' for
  equal, '*' for diff
- tAlignedSeq: reference sequence aligned on, '-' for gap base.
  qAlignedSeq, matchPattern, tAlignedSeq should have the same length

Note
========

- *coordinate*: all coordinate are 0-based, start inclusive, end exclusive. like
  slice of python sequence.
- *reverse strand*: for reads mapped on reverse strand,
    + tAlignedSeq, qAlignedSeq is reverse complemented
    + qStart, qEnd on counted from the reverse order. (in the same order as
      sequence in tAlignedSeq and qAlignedSeq

