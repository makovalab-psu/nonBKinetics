#!/usr/bin/python2
import sys, os
from sys import stdin,argv
import numpy
import random
from random import randint
import scipy.stats as stats
import re
import time

#########  Storing Constant Variable  #########
window_size = 100

#########  Arguments  #########

gff_dir = sys.argv[1]
#pacbio_file = open(sys.argv[2], 'r')

#########  Main  #########
def main(gff_dir):
    data = read_files(gff_dir)
    cleaned_windows = clean_windows(data)
    emptied_windows = empty_windows(cleaned_windows)
    windows = cleaned_windows[0] + emptied_windows
    store_windows_in_physical_memory(windows)
    
#########  Library  #########
def read_files(gff_dir): # Loop through the files in the GFF directory

    gffs = os.listdir(gff_dir)
    dictionary = {}
    for gff in gffs:
        f=open(gff_dir+gff, 'rt')
        dictionary[gff] = get_coordinates(f)
        #print('Get '+str(gff)+' Coordinates Done')
        f.close()
    return dictionary

def get_coordinates(gff): # Use the GFF format
    gff.seek(0)
    coordinates = []
    for line in gff:
        line = line.strip()
        array = line.split('\t')
        if len(array) > 5:                      # Will recognize the GFFs from nonB-DB
            strand = array[6]
            if (strand == '+') or (strand == '-'):
                chrom = array[0][3:]
                start = array[3]
                end = array[4]
                length = int(end) - int(start) + 1
                coordinates.append([start, end, chrom, length])

        if len(array) ==5:                      # Will recognize our Microsatellite formats
            chrom = array[1][3:]
            start = array[2]
            end = array[3]
            length = array[4]
            coordinates.append([start, end, chrom, length])
    return coordinates

def clean_windows(data):        # Returns non-overlapping windows containing only 1 feature # Check the centers with Bob
    chr_dic = {}
    ranks = []
    for feature in data:
        for start, stop, chrom, length in data[feature]:
            center = (int(start) + int(stop)) / 2.0
            edge = (window_size - 1) / 2.0
            #length = int(stop) - int(start) + 1
            window_start = int(center - edge)
            window_stop = int(center + edge)
            rank_start = int(int(window_start)/100)
            rank_stop = int(int(window_stop)/100)
            ranks.append(str(rank_start) + '|' + str(chrom))
            ranks.append(str(rank_stop) + '|' + str(chrom))
            if chrom in chr_dic:
                chr_dic[chrom].append([window_start, window_stop, feature, length])
            else:
                chr_dic[chrom] = [[window_start, window_stop, feature, length]]

    clean_windows = []
    for chrom in chr_dic: # Here should take into account original start and stop if length > 100
        chr_dic[chrom].sort(key=lambda x:x[0])
        previous_stop = -1
        for i,(start, stop, feature, length) in enumerate(chr_dic[chrom]):
            if i == len(chr_dic[chrom])-1:
                if start > previous_stop:
                    clean_windows.append([chrom,start,stop,length,feature])
            elif start > previous_stop and stop < chr_dic[chrom][i+1][0]:
                clean_windows.append([chrom,start,stop,length,feature])
            previous_stop = stop
    return [clean_windows,ranks]
    
def empty_windows(clean_windows):       # Find all the ranks containing no feature
    ranks = clean_windows[1]
    positions = {}
    for rank in ranks:
        position , chrom = rank.split('|')
        if chrom in positions:
            positions[chrom].append(int(position))
        else:
            positions[chrom] = [int(position)]
    uniq_positions = {}
    for chrom in positions:
        uniq_positions[chrom] = sorted(set(positions[chrom]))
    empty_positions = {}
    for chrom in uniq_positions:
        empty_positions[chrom] = set(xrange(uniq_positions[chrom][0], uniq_positions[chrom][-1])) - set(positions[chrom])
    chromosomes = []
    for chrom in empty_positions: # This and...
        chromosomes.append(chrom)
    all_empty_positions = []
    for chrom in empty_positions: # ... this and ...
        for position in empty_positions[chrom]:
            all_empty_positions.append([position, chrom])
    windows = []
    for position in all_empty_positions: # ... this should be merged in one chunk
        windows.append([position[1],position[0]*100,position[0]*100+99,100,'Empty'])

    return windows

def store_windows_in_physical_memory(windows):
    for window in windows:
        print(str(window[0])+'\t'+str(window[1])+'\t'+str(window[2])+'\t'+str(window[4])+'\t'+str(window[3]))


main(gff_dir)