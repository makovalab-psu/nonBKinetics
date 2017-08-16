#!/usr/bin/env python

import sys
import os
import re
import collections
import os.path
import itertools
import operator

from os.path import basename
from itertools import izip
from collections import defaultdict


######################################
#python reorder.py motifFileWithIPds motifFileWithErrors
motifFileWithIPds=str(sys.argv[1])
motifFileWithErrors=str(sys.argv[2])
######################################

print("motifFileWithIPds: " + motifFileWithIPds + "; " + "motifFileWithErrors: " + motifFileWithErrors)

errorDict=defaultdict()
outputEmpty=motifFileWithErrors.replace("merged_","")
outputEmpty=re.sub('$', '.collapsed', outputEmpty)

# Do not overwrite existing file
if (os.path.exists(outputEmpty)):
    print("File " + outputEmpty + " already exists. Quit.")
    sys.exit()

def load_unordered(motifFileWithErrors):
    """load the file that needs reordering

    Arguments:
        motifFileWithErrors input file handle for unordered file

    Returns:
        dictionary with the first three columns as a key and the rest as single value
    """
    with open(motifFileWithErrors) as f:
        for line_with_errors in f:
            array_with_errors = line_with_errors.rstrip().split(None, 3)
            key=str( array_with_errors[0:3] )
            values=array_with_errors[3:]
            errorDict[key]=values
    f.close()
    return errorDict

def walk_ordered(motifFileWithIPds, errorDict):
    """walk the ordered file and save the unordered values in order

    Arguments:
        motifFileWithIPds input file handle for ordered file
        errorDict dictionary with the unordered data

    Returns:
        ordered list with flattenned lines
    """
    ordered_errorDict = []
    with open(motifFileWithIPds) as f:
        i=0
        for line_with_ipds in f:
            i=i+1 #line_with_ipds number
            array_with_ipds=line_with_ipds.rstrip().split(None, 3)
            left = array_with_ipds[0:3]
            key=str(left)
            if errorDict.has_key(key):
                values=errorDict[key]
            if values:
                left.append(values[0])
            flattenned = '\t'.join( left )
            ordered_errorDict.append( flattenned )
    f.close()
    return ordered_errorDict


errorDict = load_unordered(motifFileWithErrors)
ordered_errorDict = walk_ordered(motifFileWithIPds, errorDict)


# Write results
fout = open(outputEmpty, 'w')
for item in ordered_errorDict:
    #fout.write(item + '\n')
    out_item=str(item)
    out_item=item.split(" ")[0:4] #remove last 4 columns as they are not needed in plotting script working with .collapsed files
    fout.write(str(' '.join(out_item)) + '\n')
fout.close()


