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
windowCount=defaultdict(int) #number of windows is integer
outputEmpty=motifFileWithErrors.replace("merged_","")
output_type=".collapsed"

if ("collapsed" in motifFileWithErrors):
    output_type=".reorder"

outputEmpty=re.sub('$', output_type, outputEmpty)

# Do not overwrite existing file
if (os.path.exists(outputEmpty)):
    print("File " + outputEmpty + " already exists. Quit.")
    sys.exit()

def get_window_counts(motifFileWithIPds):
    """walk the .mf file and save information about how many times each window occurs

    Arguments:
        motifFileWithIPds input .mf file

    Returns:
        dictionary with the number of times each window occurs
    """
    with open(motifFileWithIPds) as f:
        for line_with_ipds in f:
            array_with_ipds=line_with_ipds.rstrip().split(None, 3)
            left = array_with_ipds[0:3]
            key=str(left)
            windowCount[key] += 1
    f.close()
    return windowCount

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

def walk_ordered(motifFileWithIPds, errorDict, windowCount):
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
            try:
                if values:
                    number_of_window_occurencies=windowCount[key]
                    if (number_of_window_occurencies==1):
                        left.append(values[0])
                    else:
                        adjusted_rates=[] #the error rates need to be adjusted because control window is present multiple times
                        for rate in values[0].split('\t'):
                            numerical_rate=float(rate)
                            adjusted_rate=numerical_rate/float(number_of_window_occurencies) #adjustment
                            adjusted_rates.append(str(adjusted_rate))
                        adjusted_rates='\t'.join(adjusted_rates)
                        left.append(adjusted_rates)

            except UnboundLocalError:
                print("Looks like the coordinate lookup was unsuccessfull. Are you sure you are using correct input files? Aborting.")
                sys.exit()
            flattenned = '\t'.join( left )
            ordered_errorDict.append( flattenned )
    f.close()
    return ordered_errorDict


errorDict = load_unordered(motifFileWithErrors)
windowCount = get_window_counts(motifFileWithIPds)
ordered_errorDict = walk_ordered(motifFileWithIPds, errorDict, windowCount)


# Write results
fout = open(outputEmpty, 'w')
for item in ordered_errorDict:
    #fout.write(item + '\n')
    out_item=str(item)
    out_item=item.split(" ")
    if (output_type==".collapsed"):
        out_item=out_item[0:4] #remove last 4 columns as they are not needed in plotting script working with .collapsed files
    fout.write(str(' '.join(out_item)) + '\n')
fout.close()