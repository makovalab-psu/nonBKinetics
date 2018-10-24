#!/usr/bin/env python
"""
input (on stdin) a whitespace-delimited tabular file, with columns motif,
chromosome, start, end.  Motif is an N-character microsatellite motif (e.g.
TAA).  Start and end are end-inclusive coordinates of a repeat of that motif
(e.g. a thrice-repeated TAAG could have coordinates 20..31).

The goal is to collect all instances of the same motif, including shifts (but
nor reverse-complements), determine a consistent midpoint so that the motif
"lines up" (e.g. for TAA and AAT, the Ts will have the same position relative
to the midpoint), and output a 100-bp window centered at the midpoint.

Input intervals are considered to be relative to counting from zero.  Output
intervals are relative to counting from one.

For example, this input
  CCGCT chr1 75605 75654
  CGCTC chr1 23666 23705
  CCGCT chr1 54649 54673
  CCGCT chr1 20969 20993
  CGCTC chr1 81379 81453
  CTCCG chr1 86059 86133
  GCTCC chr1 22334 22363
  CGCTC chr1 53454 53523
  CTCCG chr1 50117 50156
  CTCCG chr1 60018 60047
  
will produce this output (but only the first four columns):

  CCGCT chr1 75583 75682 .......................CCGCTCCGCTCCGCTCCGCTCCGCTCCgCTCCGCTCCGCTCCGCTCCGCT...........................
  CGCTC chr1 23638 23737 .............................CGCTCCGCTCCGCTCCGCTCCgCTCCGCTCCGCTCCGCTC...............................
  CCGCT chr1 54612 54711 ......................................CCGCTCCGCTCCgCTCCGCTCCGCT.....................................
  CCGCT chr1 20932 21031 ......................................CCGCTCCGCTCCgCTCCGCTCCGCT.....................................
  CGCTC chr1 81366 81465 ..............CGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCgCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTC...........
  CTCCG chr1 86049 86148 ...........CTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCgCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCG..............
  GCTCC chr1 22300 22399 ...................................GCTCCGCTCCGCTCCgCTCCGCTCCGCTCC...................................
  CGCTC chr1 53441 53540 ..............CGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCgCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTC................
  CTCCG chr1 50087 50186 ...............................CTCCGCTCCGCTCCGCTCCgCTCCGCTCCGCTCCGCTCCG.............................
  CTCCG chr1 59983 60082 ....................................CTCCGCTCCGCTCCgCTCCGCTCCGCTCCG..................................
"""

from sys import stdin
import sys

infile = open(sys.argv[1], 'r')

def main():

        windowSize = 100
        debug = True

        # collect instances by shift-equivalent motifs
        #
        # example of contents of motifToShift:
        #   CCGCT --> (0,CCGCT)  this is the "pureMotif"
        #   CGCTC --> (1,CCGCT)
        #   GCTCC --> (2,CCGCT)
        #   CTCCG --> (3,CCGCT)
        #   TCCGC --> (4,CCGCT)

        motifToShift    = {}
        motifToIntances = {}

        for (motif,chrom,start,end) in read_motif_locations(infile):
                if (motif in motifToShift):
                        (shift,pureMotif) = motifToShift[motif]
                else:
                        motifToShift[motif] = (0,motif)
                        for ix in xrange(1,len(motif)):
                                shiftedMotif = motif[ix:] + motif[:ix]
                                motifToShift[shiftedMotif] = (ix,motif)
                        (shift,pureMotif) = (0,motif)

                if (pureMotif not in motifToIntances):
                        motifToIntances[pureMotif] = []
                motifToIntances[pureMotif] += [(shift,end-start,chrom,start,motif)]

        # choose "best" midpoint for each group of motifs, by choosing the one
        # that satisifies the most observations
        #
        # any instance, in isolation, has a natural shift of its midpoint relative
        # to its start;  for example, in CCGCTCCGCTCCGCT the midpoint is the G in
        # the motif CCGCT, so if we did not apply any further shift the midpoint
        # would be aligned to the start of GCTCC, which is the shift=2 variant of
        # the motif;  we count how may instances there are for each natural shift,
        # and choose that as the "best" shift for the group

        motifToBestShift = {}
        for pureMotif in motifToIntances:
                observedShiftCount = [0] * len(pureMotif)
                for (shift,length,chrom,start,motif) in motifToIntances[pureMotif]:
                        mid = length // 2  # (integer division)
                        naturalShift = (mid + shift) % len(pureMotif)
                        observedShiftCount[naturalShift] += 1
                        #if (debug):
                        #       print "%s %d %s -> %d+%d == %d" % (chrom,start,motif,mid,shift,naturalShift)
                motifToBestShift[pureMotif] = argmax(observedShiftCount)

        # output the windows, shifting the centers accordingly
        #
        # each instance is shifted to realign its natural shift with the best
        # shift;  conceptually this is accomplished by shifting it right N bp,
        # where N = best-natural;  however, we want to keep the result centered as
        # much as possible, so we determine the equivalent N, modulo the motif
        # length L, and within the range -L/2..+L/2

        leftHalf = windowSize // 2

        for pureMotif in motifToIntances:
                if (debug):
                        print # (separating line between non-equivalent motifs)
                bestShift = motifToBestShift[pureMotif]
                for (shift,length,chrom,start,motif) in motifToIntances[pureMotif]:
                        mid = length // 2  # (integer division)
                        naturalShift = (mid + shift) % len(motif)
                        shiftNeeded = (bestShift - naturalShift) % len(motif)

                        shiftedMid = start + mid + shiftNeeded
                        leftSide  = shiftedMid+1 - start
                        rightSide = length-1 - leftSide
                        assert (leftSide >= rightSide), "internal error, tell Bob to fix it!"
                        if (leftSide - rightSide > len(pureMotif)):
                                shiftNeeded -= len(pureMotif)
                                shiftedMid = start + mid + shiftNeeded
                                leftSide  = shiftedMid+1 - start
                                rightSide = length-1 - leftSide
                        assert (abs(leftSide-rightSide) <= len(pureMotif)), "internal error, tell Bob to fix it!"

                        windowStart = shiftedMid - leftHalf
                        windowStart += 1  # effect change from origin-0 to origin-1
                        if (debug):
                                print "%s\t%s\t%d\t%d\t%d\t%d|%d\t%s" \
                                    % (motif,chrom,windowStart,windowStart+windowSize-1,length,
                                       leftSide,rightSide,
                                       window_string(windowStart,windowSize,start,length,motif))
                        else:
                                print "%s\t%s\t%d\t%d" % (motif,chrom,windowStart,windowStart+windowSize-1)


def read_motif_locations(f):
        lineNum = 0
        for line in f:
                lineNum += 1
                line = line.strip()
                fields = line.split()

                assert (len(fields) == 4), \
                       "unexpected number of fields in line %d (expected 4, observed %d)\n%s" \
             % (lineNum,len(fields),line)

                try:
                        motif = fields[0]
                        chrom = fields[1]
                        start = int(fields[2])
                        end   = int(fields[3])
                        if (end < start): raise ValueError
                except ValueError:
                        assert (False), "unable to parse line %d\n%s" \
                                      % (lineNum,line)

                if ((end+1 - start) % len(motif) != 0):
                        assert (False), \
                               "in line %d, interval length (%d) isn't a multiple of the motif length (%d)\n%s" \
                             % (lineNum,end+1-start,len(motif),line)

                end += 1  # internally, work with half-open intervals
                yield (motif,chrom,start,end)

def argmax(values):
        # note: ties are decided in favor of earliest index
        (bestIx,maxVal) = (values[0],0)
        for (ix,val) in enumerate(values):
                if (val > maxVal): (bestIx,maxVal) = (ix,val)
        return bestIx


def window_string(windowStart,windowSize,start,length,motif):
        window = ["."] * windowSize
        for ix in xrange(length):
                if (not 0 <= start-windowStart+ix < windowSize): continue
                window[start-windowStart+ix] = motif[ix % len(motif)]
        window[windowSize//2] = window[windowSize//2].lower()
        return "".join(window)


if __name__ == "__main__": main()
