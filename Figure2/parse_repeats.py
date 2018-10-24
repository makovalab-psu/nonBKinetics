#!/usr/local/bin/python
# appends lexicographically smallest rotation at the end of the line
import sys
import os
import re
from Bio.Seq import Seq
import subprocess as sp

def RotateMe(text,mode=0,steps=1):
    # function from http://www.how2code.co.uk/2014/05/how-to-rotate-the-characters-in-a-text-string/
    # Takes a text string and rotates
    # the characters by the number of steps.
    # mode=0 rotate right
    # mode=1 rotate left
    length=len(text)

    for step in range(steps):
    # repeat for required steps
        if mode==0:
            # rotate right
            text=text[length-1] + text[0:length-1]
        else:
            # rotate left
            text=text[1:length] + text[0]
    return text
    
def SmallestRotation(seq):
    smallest=seq
    for i in range(0,len(seq)):
        actual=RotateMe(seq,0,i)
        #print ("*" + actual)
        if (actual<smallest):
            #found new minimum
            smallest=actual
    return smallest

def lexicographicallySmallestRotation(seq):     #Modified by Wil to remove reverse complement collapsing
    my_seq=Seq(seq)
    #reverse_complement=my_seq.reverse_complement()
    #reverse_complement=str(reverse_complement)

    smrt_seq=SmallestRotation(seq)
    #smrt_rev_compl_seq=SmallestRotation(reverse_complement)

    #lexicographically smallest rotation is either one of the rotations of the sequence or its reverse complement
    #if (smrt_seq < smrt_rev_compl_seq):
     #   return smrt_seq
    #else:
    #    return smrt_rev_compl_seq
    return smrt_seq

def parseTRF(file):
    number_of_sequences=0
    for line in open(file):
        li=line.strip()
        fields=li.split("\t")

        if (len(fields)==7):
            #append lexicographically smallest rotation at the end of the line
            repeat=fields[0]    #changed from 3
            #print (li + "\t" + lexicographicallySmallestRotation(repeat))
            f.write(li + "\t" + lexicographicallySmallestRotation(repeat) + "\n")

trf_file = sys.argv[1]
output = trf_file + "_output.txt"
f = open(output,"w")
parseTRF(trf_file)
f.close()
print ("Output printed into file " + output)

#split into multiple files based on the repeat representatives
args = ["awk", r'{print >  ""$8"n"}', output]
p = sp.Popen(args, stdin = sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE )
print(p.stdout.readline()) # will give you the first line of the awk output
print ("Splitting finished. Done.")