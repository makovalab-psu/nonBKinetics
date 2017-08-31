#!/usr/bin/env python
"""
Translate genomic positions in an intervals file.
"""

from sys    import argv,stdin,stdout,stderr,exit
from math   import ceil


def usage(s=None):
	message = """

usage: cat intervals | translate_positions [options]
  <translation_filename>  (required) file containing translation info (see note
                          below)
  --columns=<columns>     columns containing the chromosome, start and end; a
                          comma-separated list
                          (by default, this is column 1,2,3)
  --separator=<separator> specify separator for input and output fields
                          (default is tabs)
  --origin=zero           input and output intervals are origin-zero, half-open
                          (this is the default)
  --origin=one            input and output intervals are origin-one, closed
  --head=<number>         limit the number of input records
  --progress=<number>     periodically report how many records we've read

The input is intervals in tabular format. Output is in the same format, with
interval positions modified. Any additional columns are preserved.

The translation file describes how to translate positions.  The first two
columns indicate a position in the input genome (the genome of the intervals
coming into the program on stdin).  The third and fourth columns indicate the
corresponding position in the output genome (the genome of the intervals
coming output the program on stdout).  Note that these translation positions
are ALWAYS ORIGIN-ZER0, regardless of the --origin option.

  <-- from ---> <--- to ---->
  chr1fake 397  chr1  1785225
  chr1fake 794  chr1  3650510
  chr1fake 1191 chr1  3752149
  chr1fake 1588 chr1  3804422
   ...

Caveats:
(1) It is assumed (without checking) that each input interval lies completely
    on a single translation segment.
(2) Intervals from reverse complement strands probably will yield incorrectly
    translated intervals."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global requiredColumns,chromCol,startCol,endCol
	global debug

	# parse the command line

	translationFilename = None
	chromCol            = 0
	startCol            = 1
	endCol              = 2
	inputSeparator      = "\t"
	origin              = "zero"
	headLimit           = None
	reportProgress      = None
	debug               = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--translation=")) or (arg.startswith("--translate=")):
			translationFilename = argVal
		elif (arg.startswith("--columns=")):
			# $$$ ought to do some error checking on these
			(chromCol,startCol,endCol) = map(lambda s:int(s)-1,argVal.split(",",2))
		elif (arg.startswith("--separator=")) or (arg.startswith("--sep=")):
			inputSeparator = argVal
			if   (inputSeparator == "tab"):        inputSeparator = "\t"
			elif (inputSeparator == "space"):      inputSeparator = " "
			elif (inputSeparator == "whitespace"): inputSeparator = None
		elif (arg.startswith("--origin=")):
			origin = argVal
			if (origin == "0"): origin = "zero"
			if (origin == "1"): origin = "one"
			assert (origin in ["zero","one"]), "can't understand %s" % arg
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (translationFilename == None):
			translationFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (translationFilename == None):
		usage("you must provide a translation file")

	requiredColumns = 1 + max(chromCol,startCol,endCol)

	outputSeparator = inputSeparator
	if (inputSeparator == None): outputSeparator = " "

	# read the translation file

	translator = PositionTranslator(translationFilename)

	# process the intervals

	lineNumber = 0
	for line in stdin:
		lineNumber += 1
		if (reportProgress != None) and (lineNumber % reportProgress == 0):
			progressCount = commatize(lineNumber)
			print >>stderr, "progress: %s intervals processed" % commatize(lineNumber)

		if (headLimit != None) and (lineNumber > headLimit):
			print >>stderr, "limit of %d intervals reached" % headLimit
			break

		fields = line.split(inputSeparator)
		assert (len(fields) >= requiredColumns), \
		      "not enough fields, need %d (in input line %d): %s" \
		    % (requiredColumns,lineNumber,line)

		try:
			chrom =     fields[chromCol]
			start = int(fields[startCol])
			end   = int(fields[endCol  ])
			if (origin == "one"): start -= 1
			if (start >= end): raise ValueError
		except ValueError:
			assert (False), \
			      "bad interval in line %d: %s" % (lineNumber,line)

		try:
			(tChrom,tStart) = translator.lookup(chrom,start)
		except KeyError:
			assert (False), \
			       "%s doesn't appear in the translation file (line %d)\n%s" \
			     % (chrom,lineNumber,line)
		except ValueError:
			assert (False), \
			       "%s %d is not covered by the translation file (line %d)\n%s" \
			     % (chrom,start,lineNumber,line)

		tEnd = tStart + end-start
		if (origin == "one"): tStart += 1

		fields[chromCol] = tChrom
		fields[startCol] = str(tStart)
		fields[endCol  ] = str(tEnd)

		print outputSeparator.join(fields)


# PositionTranslator--
#	Class to translate positions from one realm to another.

class PositionTranslator(object):

	def __init__(self,filename):
		self.filename = filename

		f = file(filename,"rt")

		chromToSegments = {}

		lineNumber = 0
		for line in f:
			lineNumber += 1
			line = line.strip()

			fields = line.split()

			assert (len(fields) == 4), \
				  "wrong number of fields at line %d in %s (%d, expected %d)" \
				% (lineNumber,filename,len(fields),4)

			try:
				fromChrom =     fields[0]
				fromPos   = int(fields[1])
				toChrom   =     fields[2]
				toPos     = int(fields[3])
			except ValueError:
				assert (False), "bad translation record at line %d in %s\n%s" \
							  % (lineNumber,filename,line)

			if (fromChrom not in chromToSegments):
				chromToSegments[fromChrom] = []
			chromToSegments[fromChrom] += [(fromPos,toChrom,toPos)]

		f.close()

		for fromChrom in chromToSegments:
			chromToSegments[fromChrom].sort()
			segments = chromToSegments[fromChrom]
			for (ix,(fromPos,toChrom,toPos)) in enumerate(segments):
				if (ix == 0): continue
				(prevPos,_,_) = segments[ix-1]
				if (prevPos == fromPos):
					assert (prevPos != fromPos), \
					       "translation positions %s %d occurs more than once in %s" \
					     % (fromChrom,fromPos,filename)

		self.chromToSegments = chromToSegments


	def lookup(self,chrom,pos):
		# $$$ there are smarter, faster ways to do this

		if (chrom not in self.chromToSegments): raise KeyError
		segments = self.chromToSegments[chrom]

		prev = None
		for seg in segments:
			(fromPos,_,_) = seg
			if (fromPos > pos): break
			prev = seg
		if (prev == None): raise ValueError

		(fromPos,toChrom,toPos) = prev
		return (toChrom,toPos+pos-fromPos)


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


# commatize--
#	Convert a numeric string into one with commas.

def commatize(s):
	if (type(s) != str): s = str(s)
	(prefix,val,suffix) = ("",s,"")
	if (val.startswith("-")): (prefix,val) = ("-",val[1:])
	if ("." in val):
		(val,suffix) = val.split(".",1)
		suffix = "." + suffix

	try:    int(val)
	except: return s

	digits = len(val)
	if (digits > 3):
		leader = digits % 3
		chunks = []
		if (leader != 0):
			chunks += [val[:leader]]
		chunks += [val[ix:ix+3] for ix in xrange(leader,digits,3)]
		val = ",".join(chunks)

	return prefix + val + suffix


if __name__ == "__main__": main()
