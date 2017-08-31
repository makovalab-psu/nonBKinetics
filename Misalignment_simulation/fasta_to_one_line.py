#!/usr/bin/env python
"""
Convert a fasta file to a file with each sequence on a single line (including
header)
"""

from sys  import argv,stdin,stderr,exit
from math import ceil

def main():

	# parse the command line

	separator = "#"
	toUpper   = False
	noArrow   = False
	nameParse = None
	headLimit = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--noseparator") or (arg == "--nosep"):
			separator = None
		elif (arg.startswith("--separator=")):
			separator = argVal
		elif (arg == "--upper"):
			toUpper = True
		elif (arg == "--noarrow"):
			noArrow = True
		elif (arg == "--nameparse=darkspace"):
			nameParse = "darkspace"
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--")):
			assert (False), "unknown argument: %s" % arg
		else:
			assert (False), "unknown argument: %s" % arg

	# process the sequences

	if (noArrow): arrow = ""
	else:         arrow = ">"

	seqNum = 0
	for (name,seq) in fasta_sequences(stdin):
		seqNum += 1
		if (headLimit != None) and (seqNum > headLimit):
			print >>stderr, "limit of %d sequences reached" % headLimit
			break

		if (nameParse == "darkspace"):
			name = name.split()[0]
		if (separator != None) and (separator in name):
			name = name.replace(separator,"_")
		if (toUpper): seq = seq.upper()
		if (separator == None): print "%s%s %s"    % (arrow,name,seq)
		else:                   print "%s%s %s %s" % (arrow,name,separator,seq)


# fasta_sequences--
#	Read the fasta sequences from a file

def fasta_sequences(f):
	seqName = None
	seqNucs = None

	for line in f:
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				yield (seqName,"".join(seqNucs))
			seqName = line[1:].strip()
			seqNucs = []
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seqNucs += [line]

	if (seqName != None):
		yield (seqName,"".join(seqNucs))


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


if __name__ == "__main__": main()

