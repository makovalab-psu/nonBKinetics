#!/usr/bin/env python
"""
Input 1 is a list of window intervals, one interval per line,
        chromosome start end
and must be in sorted order.  Intervals are origin-zero half-open.

Input 2 is a list of positions and values,
        chromsome position value
and must be in sorted order.  Positions are origin-zero.

Sorting can be accomplished like this:
    cat windows.dat | env LC_ALL=C sort -k 1,1d -k 2,2n > windows.sorted
"""

from sys import stdin,argv
import numpy

def main():
        assert (len(argv) == 3), \
            "you must give me exactly two file names, windows and values"

        windowsFilename = argv[1]
        valuesFilename  = argv[2]

        windowReader = read_windows(windowsFilename)
        (windowChrom,start,end,feature,length) = windowReader.next()
        assert (windowChrom != None), "%s contains no windows" % windowsFilename
        windowValues = new_window(windowChrom,start,end,feature,length)
        valueReader  = read_values(valuesFilename)
        (valChrom,valPos,strand,val) = valueReader.next()
        assert (valChrom != None), "%s contains no values" % valuesFilename

        while (valChrom != None):
                # is value before current window?

                if (valChrom < windowChrom) \
                  or ((valChrom == windowChrom) and (valPos < start)):
                        (valChrom,valPos,strand,val) = valueReader.next()
                        continue

                # is value within current window?
                if (valChrom == windowChrom) and (start <= valPos <= end):
                        windowValues[valPos-start] = val
                        (valChrom,valPos,strand,val) = valueReader.next()
                        continue

                # otherwise, value is after current window
                #print(windowValues)
                print_window(windowChrom,start,end,feature,length,windowValues)

                (windowChrom,start,end,feature,length) = windowReader.next()
                if (windowChrom == None):
                        break
                windowValues = new_window(windowChrom,start,end,feature,length)

        if (windowChrom != None):
                print_window(windowChrom,start,end,feature,length,windowValues)

def read_windows(filename):
        f = file(filename,"rt")

        prevChrom = prevEnd = None

        lineNumber = 0
        for line in f:
                lineNumber += 1
                line = line.strip()

                fields = line.split()

                assert (len(fields) == 5), \
                      "wrong number of fields in line %d of %s\n%s" \
                    % (lineNumber,filename,line)

                try:
                        chrom   = fields[0]
                        start   = int(fields[1])
                        end     = int(fields[2])
                        feature = fields[3]
                        length  = fields[4]

                except ValueError:
                        assert (False), \
                               "bad window (line %d in %s)\n%s" \
                            % (lineNumber,filename,line)

                assert (end != start), \
                       "empty window (line %d in %s)\n%s" \
                    % (lineNumber,filename,line)

                assert (end > start), \
                       "negative-length window (line %d in %s)\n%s" \
                    % (lineNumber,filename,line)

                if (prevChrom != None):
                        if (chrom < prevChrom) \
                          or ((chrom == prevChrom) and (start < prevEnd)):
                                assert (False), \
                                       "windows out of order or overlapping (at line %d in %s)\n%s" \
                                    % (lineNumber,filename,line)

                yield (chrom,start,end,feature,length)
                (prevChrom,prevEnd) = (chrom,end)

        yield (None,None,None,None,None)
        f.close()

def read_values(filename):
        f = file(filename,"rt")
        f.readline()
        prevChrom = prevPos = None

        lineNumber = 0
        for line in f:
                lineNumber += 1
                line = line.strip()

                fields = line.split()
                assert (len(fields) == 4), \
                      "wrong number of fields in line %d of %s\n%s" \
                    % (lineNumber,filename,line)

                try:
                        chrom  = fields[0]
                        pos    = int  (fields[1])
                        strand = int(fields[2])
                        val    = float(fields[3])
                except ValueError:
                        assert (False), \
                               "bad position/value (line %d in %s)\n%s" \
                            % (lineNumber,filename,line)
                if strand == 1:                                                                 #########STRAND SPECIFICITY IS HERE############
                        if (prevChrom != None):
                                if (chrom < prevChrom) \
                                  or ((chrom == prevChrom) and (pos < prevPos)):
                                        assert (False), \
                                               "positions out of order (at line %d in %s)\n%s" \
                                            % (lineNumber,filename,line)

                        yield (chrom,pos,strand,val)
                        (prevChrom,prevPos) = (chrom,pos)
        yield (None,None,None,None)
        f.close()


def new_window(chrom,start,end,feature,length):
        return [None] * (end-start+1)

def print_window(chrom,start,end,feature,length,values,defaultVal=numpy.nan):
        print chrom,start,end,feature,length,
        for pos in xrange(start,end+1):
                if (values[pos-start] != None):
                        print values[pos-start],
                else:
                        print defaultVal,
        print



if __name__ == "__main__": main()