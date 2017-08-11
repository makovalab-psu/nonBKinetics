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

def get_arg(index):
    try:
        sys.argv[index]
    except IndexError:
        return False
    else:
        return sys.argv[3]

######################################
#This script will convert filtered gff files back to mf files
#python gff2mf.py originalmfFile gffFile
originalMfFile=str(sys.argv[1])
gffFile=str(sys.argv[2])
mode=get_arg(3) #optional argument 'keep|remove' picks the subsetting mode
######################################

gffDict=defaultdict()
outputEmpty=gffFile.replace(".gff",".mf")
outputEmpty=basename(outputEmpty)

if(get_arg(3)==False):
    mode="keep"
else:
    mode=get_arg(3)

print("originalMfFile: " + originalMfFile + "; " + "gffFile: " + gffFile)
print("Note that the script will only work if two files have same formatting (100bp vs FeatureOnly)")

if (mode=="keep"):
    print("Mode is to subset to MATCHING LINES.")
    print("Output written to " + outputEmpty)
else:
    print("Mode is to subset to NON-MATCHING LINES.")
    outputEmpty=outputEmpty.replace(".mf","_complement.mf")
    print("Output written to " + outputEmpty)


# Do not overwrite existing file
if (os.path.exists(outputEmpty)):
    print("File " + outputEmpty + " already exists. Quit.")
    sys.exit()

def load_unoriginal(gffFile):
    """load the file that represents the needed subset

    Arguments:
        gffFile input file handle 

    Returns:
        dictionary with the first three columns as a key and the rest as single value
    """
    with open(gffFile) as f:
        for line_with_errors in f:
            array_with_errors = line_with_errors.rstrip().split("\t")
            #coordinates chr start end
            key=str(array_with_errors[0])+" "+str(array_with_errors[3])+" "+str(array_with_errors[4])

            values=True
            gffDict[key]=values
    f.close()
    return gffDict

def walk_original(originalMfFile, gffDict):
    """walk the original file and subset

    Arguments:
        originalMfFile input file handle
        dictionary with the subsetted data

    """
    original_gffDict = []
    fout = open(outputEmpty, 'w')
    with open(originalMfFile) as f:
        for line_with_ipds in f:
            array_with_ipds=line_with_ipds.rstrip().split("\t", 3)
            left = array_with_ipds[0:3]
            key=str(left)
            flattenned = ' '.join( left )
            if (gffDict.has_key(flattenned)==True):
                if (mode=="keep"):
                    fout.write(line_with_ipds)
            else:
                if (mode!="keep"):
                    fout.write(line_with_ipds)
    f.close()
    fout.close()

gffDict = load_unoriginal(gffFile)
walk_original(originalMfFile, gffDict)



