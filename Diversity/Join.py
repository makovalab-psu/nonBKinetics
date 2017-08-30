import sys
import re
#This script joins the stitched blocks that share extremities

infile = open(sys.argv[1], 'rt') #stitched file
infile2 = open(sys.argv[2], 'rt') #original intervals to map
outfile = open(sys.argv[1]+'.joined', 'w')

intervalmap = {}

for line in infile2:
        line = line.strip()
        chrom, start, stop, AF, ref, alt = line.split('\t')
        if AF > 0.05:
                key = chrom+'|'+start+'|'+stop
                intervalmap[key] = [ref,alt]


oldstart = 0
oldstop = 0
hom = ''
pan = ''
gor = ''
pon = ''
nom = ''

infile.readline() #remove header

for line in infile:
        line = line.strip()
        #chrom, newstart, newstop, strand, score, name, gor, hom, nom, pan, pon = line.split('\t')
        chrom, newstart, newstop = line.split('\t')[0:3]

        linelen = len(line.split('\t'))

        if linelen != 11: #missing species to the right end
                missing =  11 - linelen
                for i in range(0,missing):
                        line = line + '\t'


        if newstart != oldstop:
                key = chrom+'|'+str(oldstart)+'|'+str(oldstop) #print previous window
                if key in intervalmap:
                        ref, alt = intervalmap[key]
                        outfile.write(chrom+'\t'+oldstart+'\t'+oldstop+'\t'+hom+'\t'+pan+'\t'+gor+'\t'+pon+'\t'+nom+'\t'+ref+'\t'+alt+'\n')
                oldstart = newstart
                oldstop = newstop
                gor, hom, nom, pan, pon = line.split('\t')[-5:]
                
        else: #newstart == oldstop

                key = chrom+'|'+str(oldstart)+'|'+str(oldstop)
                if key in intervalmap:
                        ref, alt = intervalmap[key]
                        outfile.write(chrom+'\t'+oldstart+'\t'+oldstop+'\t'+hom+'\t'+pan+'\t'+gor+'\t'+pon+'\t'+nom+'\t'+ref+'\t'+alt+'\n')
                        oldstart = newstart
                        oldstop = newstop
                        gor, hom, nom, pan, pon = line.split('\t')[-5:]
                else: # merge windows
                        newgor, newhom, newnom, newpan, newpon = line.split('\t')[-5:]
                        oldstop = newstop
                        hom = str(hom)+str(newhom)
                        pan = str(pan)+str(newpan)
                        gor = str(gor)+str(newgor)
                        pon = str(pon)+str(newpon)
                        nom = str(nom)+str(newnom)



#print(oldstart,oldstop,hom,pan,gor,pon,nom)
outfile.write(chrom+'\t'+str(oldstart)+'\t'+str(oldstop)+'\t'+hom+'\t'+pan+'\t'+gor+'\t'+pon+'\t'+nom+'\t'+ref+'\t'+alt+'\n')