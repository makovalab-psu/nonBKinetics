#!/usr/bin/python2
import sys, os

file = open(sys.argv[1], 'r')

features = {}

for line in file:
        line = line.strip()
        array = line.split(' ')
        chrom, start, stop, feature, length = array[0], array[1], array[2], array[3], array[4]
        IPDs = array[5:]
        if feature in features:
                features[feature].append([chrom,start,stop,length,IPDs])
        else:
                features[feature] = [[chrom,start,stop,length,IPDs]]

for feature in features:
        out = open(feature, 'w')
        for window in features[feature]:
                out.write(str(window[0])+'\t'+str(window[1])+'\t'+str(window[2])+'\t'+str(window[3]))
                for IPD in window[4]:
                        out.write('\t'+str(IPD))
                out.write('\n')
