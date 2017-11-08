import sys
from string import maketrans
#This script polarize the indels from the .alleles


infile = open(sys.argv[1], 'rt')
nucl_mapping = maketrans("actg","ACTG")


def custom_translate(input_string):
   return input_string.translate(nucl_mapping)



for line in infile:
	line = line.strip()
#	print(line)
#	chrom, start, end, gor, hom, nom, pan, pon, ref, alt = line.split('\t')
        chrom, start, end, hom, pan, gor, pon, nom, ref, alt = line.split('\t')

	#print(gor, hom, nom, pan, pon)

	originalstart = str(int(start)+6)

        hom = custom_translate(hom)
        pan = custom_translate(pan)
        gor = custom_translate(gor)
        pon = custom_translate(pon)
	nom = custom_translate(nom)

	#species = [pan, gor, pon, nom]

	lenref = len(ref)
	#lenalt = len(alt) - 1
	lenalt = len(alt)

#	print(ref,lenref)
#	print(alt,lenalt)

	if lenref > lenalt:
		gap = lenref - lenalt
		#print(alt)
		alt += '-'*gap
		#print(alt)

	if lenref < lenalt:
		gap = lenalt - lenref
		#print(ref)
		ref += '-'*gap
		#print(ref)
		altgor = gor[:6] + '-'*gap + gor[6:]
		altnom = nom[:6] + '-'*gap + nom[6:]
		altpan = pan[:6] + '-'*gap + pan[6:]
		altpon = pon[:6] + '-'*gap + pon[6:]
	
	newref = hom[:5] + ref + hom[-5:]
	#newalt = hom[:5] + alt[1:] + hom[-5:]	
	newalt = hom[:5] + alt + hom[-5:]

	print(str('ref, alt, hom, newref, newalt, pan, gor, pon, nom'))
	print(ref, alt, hom, newref, newalt, pan, gor, pon, nom)
	#print(ref,hom,newref)
	#print(alt,hom,newalt)

	species = [pan, gor, pon, nom]

	count_ref = 0
	count_alt = 0

	for specie in species: #not grammatically correct but whatever
		if newref == specie:
			count_ref +=1
		if newalt == specie:
			count_alt += 1

        if lenref < lenalt:
		#print('TWO',ref, alt, hom, newref, newalt, altpan, altgor, altpon, altnom)
                altspecies = [altpan,altgor,altpon,altnom]
		for specie in altspecies:
	                if newref == specie:
       	 	               count_ref +=1
       	        	if newalt == specie:
      				count_alt += 1



        #if count_ref >= 2 and count_alt >= 2:
         #       print('STANDING VARIATION')
       	  #      print(newref,newalt,species,count_ref,count_alt,)

	#if lenref < lenalt:
	#	print(newref,newalt,species,count_ref,count_alt,)
	if count_ref >= 2 and count_alt < 2:
		if lenref > lenalt and count_ref > count_alt:
			print(chrom+'\t'+originalstart+'\tDEL\t in alt')
                if lenref < lenalt and count_ref > count_alt:
                        print(chrom+'\t'+originalstart+'\tINS\t in alt')

        if count_ref < 2 and count_alt >= 2:
		if lenref > lenalt and count_ref < count_alt:
			print(chrom+'\t'+originalstart+'\tINS\t in ref')
		if lenref < lenalt and count_ref < count_alt:
			print(chrom+'\t'+originalstart+'\tDEL\t in ref')

