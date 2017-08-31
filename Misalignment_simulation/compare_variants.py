#!/usr/bin/env python
"""
Compare variants from two variant-call files
"""

from sys  import argv,stdin,stdout,stderr,exit
from math import ceil
from sets import Set
from gzip import open as gzip_open


def usage(s=None):
	message = """

usage: compare_variants <true_variants> <observed_variants> [options]
  --model=1            compare variants as per 'error model 1' (see note below)
                       (this is the default)
  --model=12           compare variants as per 'error model 2' (see note below)
  <true_variants>      (required) file containing information about 'true'
                       variants (see note below); this can be a gzip file
  <observed_variants>  file containing information about 'observed' variants
                       (see note below); this can be a gzip file; if this is
                       absent, it is read from stdin
  --regions=<file>     (recommended) regions to restrict computation to; this
                       is also used to determine the denominator for reported
                       statistics; one interval per line, <chrom> <start> <end>
                       (by default, origin-zero half-open)
  --regions=origin1:<file> Same as --regions, but the intervals are origin-one
                       closed
  --identifier=<name>  identifying name for this data set (used in output
                       tables)
  --head=<number>      limit the number of observed variants we'll process
  --progress=<number>  periodically report how many observed variants we've
                       processed

The input is the chrom, pos, ref, and alt fields from vcf format, with an
optional header line may be included. As such, the positions are presumed to
be origin-1 (however, since both 'true' and 'observed' files has the same
coordiniates, it usually won't matter).

  #CHROM POS    REF ALT
  chr7   136988 G   A
  chr7   136999 GG  G
  chr7   137005 T   G,C
  chr7   150919 C   CT
  chr7   151225 A   G
   ...

The 'true' variants file is considered to be the ground truth.  The 'observed'
variants are compared against the truth.

Error Models:

Error models determine how certain differences between truth and observation
are treated.  In general we will report observations that differ from truth
as false positives (FP).  But the model determines what we consider to be
different.

Error model 1:
(note that this model was implemented for a specific project, and may not be
approriate for any other project)

This model reports how many observed variants are not truth.  If multiple
variants are observed for a single position, they are all reported (except any
of them that are truth).

Error model 2:
(note that this model was implemented for a specific project, and may not be
approriate for any other project)

In this model our primary concern is reporting observation positions that have
a variant of a diffent type (mismatch, insertion, or deletion) than the
reference at the same position.  Thus if the truth and observation have
different insertions at the same position, we do *not* report that as a false
positive."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global headLimit,reportProgress,debug

	# parse the command line

	errorModel          = "1"
	filenames           = []
	truthFilename       = None
	observationFilename = None
	regionsFilename     = None
	regionsOrigin       = None
	idName              = None
	headLimit           = None
	reportProgress      = None
	debug               = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--model=")):
			if (argVal in ["1","2"]):
				errorModel = argVal
			else:
				usage("unrecognized model in %s" % arg)
		elif (arg.startswith("--truth=")) \
		  or (arg.startswith("--true=")):
			truthFilename = argVal
		elif (arg.startswith("--observed=")) \
		  or (arg.startswith("--observation=")) \
		  or (arg.startswith("--observations=")) \
		  or (arg.startswith("--obs=")):
			observationFilename = argVal
		elif (arg.startswith("--regions=origin0:")) \
		  or (arg.startswith("--intervals=origin0:")) \
		  or (arg.startswith("--universe=origin0:")):
			argVal = arg.split(":",1)[1]
			regionsFilename = argVal
			regionsOrigin = "zero"
		elif (arg.startswith("--regions=origin1:")) \
		  or (arg.startswith("--intervals=origin1:")) \
		  or (arg.startswith("--universe=origin1:")):
			argVal = arg.split(":",1)[1]
			regionsFilename = argVal
			regionsOrigin = "one"
		elif (arg.startswith("--regions=")) \
		  or (arg.startswith("--intervals=")) \
		  or (arg.startswith("--universe=")):
			regionsFilename = argVal
			regionsOrigin = "zero"
		elif (arg.startswith("--identifier=")) or (arg.startswith("--id=")):
			idName = argVal
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
		else:
			filenames += [arg]

	if (truthFilename == None):
		truthFilename = filenames.pop(0)
	if (observationFilename == None):
		observationFilename = filenames.pop(0)
	if (filenames != []):
		usage("unrecognized option: %s (extra file name?)" % filenames[0])

	if (truthFilename == None):
		usage("you must provide a 'truth' variants file")

	if ("fp" in debug): debug += ["FP"]
	if ("tp" in debug): debug += ["TP"]
	if ("fn" in debug): debug += ["FN"]

	if (idName == None): idNameTag = ""
	else:                idNameTag = idName + " "

    # read regions of interest

	regions = None
	if (regionsFilename != None):
		if (regionsFilename.endswith(".gz")) or (regionsFilename.endswith(".gzip")):
			regionsF = gzip_open(regionsFilename,"rt")
		else:
			regionsF = file(regionsFilename,"rt")

		regions = {}
		for (lineNumber,chrom,start,end) in read_intervals(regionsF,origin=regionsOrigin):
			if (chrom not in regions): regions[chrom] = []
			regions[chrom] += [(start,end)]

		regionsF.close()

		for chrom in regions:
			regions[chrom] = non_overlapping_intervals(regions[chrom])

	# collect the ground truth

	if (truthFilename.endswith(".gz")) or (truthFilename.endswith(".gzip")):
		truthF = gzip_open(truthFilename,"rt")
	else:
		truthF = file(truthFilename,"rt")

	truth = {}
	for (lineNumber,chrom,pos,ref,alts) in read_variants(truthF):
		if (regions != None):
			if   (chrom not in regions) \
			  or (not in_intervals(chrom,regions[chrom],pos)):
				continue

		if (chrom not in truth): truth[chrom] = {}
		assert (pos not in truth[chrom]), \
		       "truth contains more than one variant record for %s %d (second is at line %d)" \
		     % (chrom,pos,lineNumber)
		truth[chrom][pos] = (ref,alts)

	truthF.close()

	universeLen = None
	if (regions != None):
		universeLen = 0
		for chrom in regions:
			universeLen += sum([e-s for (s,e) in regions[chrom]])

	# process the observations

	if (observationFilename == None):
		obsF = stdout
	elif (observationFilename.endswith(".gz")) or (observationFilename.endswith(".gzip")):
		obsF = gzip_open(observationFilename,"rt")
	else:
		obsF = file(observationFilename,"rt")

	if (errorModel == "1"):
		(truePositives,falsePositives) = error_model_1(obsF,truth,regions)
	elif (errorModel == "2"):
		(truePositives,falsePositives) = error_model_2(obsF,truth,regions)
	else:
		assert (False), "error model \"%s\" is not implemented" % errorModel

	if (obsF != stdout): obsF.close()

	# report results

	if (universeLen == None):
		print "%sfalseMismatches: %d" % (idNameTag,falsePositives["MM"])
		print "%sfalseDeletions:  %d" % (idNameTag,falsePositives["DEL"])
		print "%sfalseInsertions: %d" % (idNameTag,falsePositives["INS"])
	else:
		print "%sfalseMismatches: %d/%d %0.4f%%" \
		    % (idNameTag,falsePositives["MM"] ,universeLen,
		       100.0*falsePositives["MM"]/universeLen)
		print "%sfalseDeletions:  %d/%d %0.4f%%" \
		    % (idNameTag,falsePositives["DEL"],universeLen,
		       100.0*falsePositives["DEL"]/universeLen)
		print "%sfalseInsertions: %d/%d %0.4f%%" \
		    % (idNameTag,falsePositives["INS"],universeLen,
		       100.0*falsePositives["INS"]/universeLen)

	# report false negatives

	if ("FN" in debug):
		assert (model == "2"), "FN not implemented yet for model %s" % model
		for (chrom,pos,ref,alts) in false_negatives(truth,truePositives):
			print >>stderr, "FN %s %d %s %s" % (chrom,pos+1,ref,",".join(alts))


# error_model_1--

def error_model_1(obsF,truth,regions=None):
	falsePositives = {"MM":0 , "INS":0 , "DEL": 0}
	truePositives  = {}

	obsNumber = 0
	for (lineNumber,chrom,pos,ref,alts) in read_variants(obsF):
		obsNumber += 1
		if (reportProgress != None) and (obsNumber % reportProgress == 0):
			progressCount = commatize(obsNumber)
			print >>stderr, "progress: %s observations processed" % commatize(obsNumber)

		if (headLimit != None) and (obsNumber > headLimit):
			print >>stderr, "limit of %d observations reached" % headLimit
			break

		if (regions != None):
			if   (chrom not in regions) \
			  or (not in_intervals(chrom,regions[chrom],pos)):
				continue

		# if the ref and alt are both the same, this is not a variant;  note
		# that it might be a false negative, but we don't collect those here

		if (len(alts) == 1) and (alts[0] == ref):
			continue

		alts = tuple(alt for alt in alts if (alt != ref))

		# if the chromosome or position are not in the truth, it's a FP

		if (chrom not in truth) or (pos not in truth[chrom]):
			for (fpType,ref,alt) in falsities(ref,alts):
				falsePositives[fpType] += 1
				if ("FP" in debug):
					print >>stderr, "FP %s %s %d %s %s" \
					              % (fpType,chrom,pos+1,ref,alt)
			continue

		# if the references are different ...
		# ... if they are both single nts, it's an error
		# ... if truth is a single nt (and obs isn't), it's a FP (false deletion)
		# ... otherwise it will be a false negative, so we ignore it here

		(trueRef,trueAlts) = truth[chrom][pos]
		trueAlts = tuple(alt for alt in trueAlts if (alt != trueRef))

		if (ref != trueRef):
			if (len(trueRef) == 1) and (len(ref) == 1):
				assert (False), \
				       "ref difference at %s %d (true=%s vs observed=%s)" \
				     % (chrom,pos+1,trueRef,ref)
			if (len(trueRef) == 1):
				for (fpType,ref,alt) in falsities(ref,alts):
					falsePositives[fpType] += 1
					if ("FP" in debug):
						print >>stderr, "FP %s %s %d %s %s" \
						              % (fpType,chrom,pos+1,ref,alt)
			continue

		# if the references are the same but not a single nt, we consider this
		# as a TP (true deletion); note that we are ignoring the possibility
		# that the alts might be different
		
		if (len(trueRef) > 1):
			if (chrom not in truePositives): truePositives[chrom] = Set()
			truePositives[chrom].add(pos)
			if ("TP" in debug):
				print >>stderr, "TP %s %d %s %s" % (chrom,pos+1,ref,",".join(alts))
			continue

		# the references are the same single nt;  if they have have any alt in
		# common, we consider it a TP

		if (has_common_alt(trueAlts,alts)):
			if (chrom not in truePositives): truePositives[chrom] = Set()
			truePositives[chrom].add(pos)
			if ("TP" in debug):
				print >>stderr, "TP %s %d %s %s" % (chrom,pos+1,ref,",".join(alts))
			continue

		# otherwise, the references are the same single nt but they have have
		# no alt in common;  but if the both have a mismatch, or both have an
		# insertion, we consider it a TP

		if (has_common_alt_type(trueAlts,alts)):
			if (chrom not in truePositives): truePositives[chrom] = Set()
			truePositives[chrom].add(pos)
			if ("TP" in debug):
				print >>stderr, "TP %s %d %s %s" % (chrom,pos+1,ref,",".join(alts))
			continue

		# otherwise, it's a FP

		for (fpType,ref,alt) in falsities(ref,alts):
			falsePositives[fpType] += 1
			if ("FP" in debug):
				print >>stderr, "FP %s %s %d %s %s" \
				              % (fpType,chrom,pos+1,ref,alt)

	return (truePositives,falsePositives)

# error_model_2--

def error_model_2(obsF,truth,regions=None):
	fpMismatches = fpInsertions = fpDeletions = 0
	truePositives = {}

	obsNumber = 0
	for (lineNumber,chrom,pos,ref,alts) in read_variants(obsF):
		obsNumber += 1
		if (reportProgress != None) and (obsNumber % reportProgress == 0):
			progressCount = commatize(obsNumber)
			print >>stderr, "progress: %s observations processed" % commatize(obsNumber)

		if (headLimit != None) and (obsNumber > headLimit):
			print >>stderr, "limit of %d observations reached" % headLimit
			break

		if (regions != None):
			if   (chrom not in regions) \
			  or (not in_intervals(chrom,regions[chrom],pos)):
				continue

		# if the ref and alt are both the same, this is not a variant;  note
		# that it might be a false negative, but we don't collect those here

		if (len(alts) == 1) and (alts[0] == ref):
			continue

		alts = tuple(alt for alt in alts if (alt != ref))

		# if the chromosome or position are not in the truth, it's a FP

		if (chrom not in truth) or (pos not in truth[chrom]):
			if (len(ref) > 1):
				fpDeletions += 1
				fpType = "DEL"
			elif (has_insertion(alts)):
				fpInsertions += 1
				fpType = "INS"
			else:
				fpMismatches += 1
				fpType = "MM"
			if ("FP" in debug):
				print >>stderr, "FP %s %s %d %s %s" \
				              % (fpType,chrom,pos+1,ref,",".join(alts))
			continue

		# if the references are different ...
		# ... if they are both single nts, it's an error
		# ... if truth is a single nt (and obs isn't), it's a FP (false deletion)
		# ... otherwise it will be a false negative, so we ignore it here

		(trueRef,trueAlts) = truth[chrom][pos]
		trueAlts = tuple(alt for alt in trueAlts if (alt != trueRef))

		if (ref != trueRef):
			if (len(trueRef) == 1) and (len(ref) == 1):
				assert (False), \
				       "ref difference at %s %d (true=%s vs observed=%s)" \
				     % (chrom,pos+1,trueRef,ref)
			if (len(trueRef) == 1):
				if ("FP" in debug):
					print >>stderr, "FP %s %s %d %s %s" \
					              % ("DEL",chrom,pos+1,ref,",".join(alts))
				fpDeletions += 1
			continue

		# if the references are the same but not a single nt, we consider this
		# as a TP (true deletion); note that we are ignoring the possibility
		# that the alts might be different
		
		if (len(trueRef) > 1):
			if (chrom not in truePositives): truePositives[chrom] = Set()
			truePositives[chrom].add(pos)
			if ("TP" in debug):
				print >>stderr, "TP %s %d %s %s" % (chrom,pos+1,ref,",".join(alts))
			continue

		# the references are the same single nt;  if they have have any alt in
		# common, we consider it a TP

		if (has_common_alt(trueAlts,alts)):
			if (chrom not in truePositives): truePositives[chrom] = Set()
			truePositives[chrom].add(pos)
			if ("TP" in debug):
				print >>stderr, "TP %s %d %s %s" % (chrom,pos+1,ref,",".join(alts))
			continue

		# otherwise, the references are the same single nt but they have have
		# no alt in common;  but if the both have a mismatch, or both have an
		# insertion, we consider it a TP

		if (has_common_alt_type(trueAlts,alts)):
			if (chrom not in truePositives): truePositives[chrom] = Set()
			truePositives[chrom].add(pos)
			if ("TP" in debug):
				print >>stderr, "TP %s %d %s %s" % (chrom,pos+1,ref,",".join(alts))
			continue

		# otherwise, it's a FP

		if (has_insertion(alts)):
			fpInsertions += 1
			fpType = "INS"
		else:
			fpMismatches += 1
			fpType = "MM"
		if ("FP" in debug):
			print >>stderr, "FP %s %s %d %s %s" \
			              % (fpType,chrom,pos+1,ref,",".join(alts))

	falsePositives = {"MM":fpMismatches , "INS":fpInsertions , "DEL": fpDeletions}
	return (truePositives,falsePositives)


# falsities--

def falsities(ref,alts):
	if (type(alts) not in (tuple,list)):
		alts = (alts,)

	refLen = len(ref)
	if (refLen == 1):
		for alt in alts:
			if (len(alt) == 1): yield ("MM",ref,alt)
			else:               yield ("INS",ref,alt)
	else: # if (refLen > 1):
		for alt in alts:
			if (len(alt) < refLen):
				yield ("DEL",ref,alt)
			elif (len(alt) > refLen):
				yield ("INS",ref,alt)
			else: # if (len(alt) == refLen):
				yield ("MM",ref,alt)


# false_negatives--

def false_negatives(truth,truePositives):
	for chrom in truth:
		if (chrom not in truePositives):
			for pos in truth[chrom]:
				(ref,alts) = truth[chrom][pos]
				yield (chrom,pos,ref,alts)
		else:
			for pos in truth[chrom]:
				if (pos not in truePositives[chrom]):
					(ref,alts) = truth[chrom][pos]
					yield (chrom,pos,ref,alts)

# has_common_alt--

def has_common_alt(trueAlts,alts):
	commonAlt = None
	for trueAlt in trueAlts:
		for alt in alts:
			if (alt == trueAlt): return True
	return False


# has_common_alt_type--

def has_common_alt_type(trueAlts,alts):
	trueHasMismatch  = False
	trueHasInsertion = False
	for alt in trueAlts:
		if (len(alt) == 1): trueHasMismatch  = True
		else:               trueHasInsertion = True
	for alt in alts:
		if (len(alt) == 1) and (trueHasMismatch):  return True
		elif                   (trueHasInsertion): return True
	return False


# has_insertion--

def has_insertion(alts):
	for alt in alts:
		if (len(alt) > 1): return True
	return False


# read_variants--

def read_variants(f):
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line.startswith("#")):
			continue

		fields = line.split()
		assert (len(fields) == 4), \
		      "wrong number of columns at line %d (%d, expected %d):\n%s" \
		    % (lineNumber,len(fields),4,line)

		try:
			chrom =     fields[0]
			pos   = int(fields[1]) - 1   # convert to origin-zero
			ref   =     fields[2]
			alts  =     fields[3].split(",")
		except ValueError:
			assert (False), "bad variant (at line %d):\n%s" \
			              % (lineNumber,line)

		yield (lineNumber,chrom,pos,ref,alts)


# read_intervals--

def read_intervals(f,origin="zero"):
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line.startswith("#")):
			continue

		fields = line.split()
		assert (len(fields) >= 3), \
		      "not enough columns at line %d (%d, expected at least %d):\n%s" \
		    % (lineNumber,len(fields),3,line)

		try:
			chrom =     fields[0]
			start = int(fields[1])
			end   = int(fields[2])
			if (start >= end): raise ValueError
			if (origin == "one"): start -= 1
		except ValueError:
			assert (False), "bad interval (at line %d):\n%s" \
			              % (lineNumber,line)

		yield (lineNumber,chrom,start,end)


# non_overlapping_intervals--
#	merge into non-overlapping intervals

def non_overlapping_intervals(intervals):
	if (intervals == []): return []

	intervals.sort()
	overlaps = []
	start = end = None
	for (s,e) in intervals:
		if (start == None):
			(start,end) = (s,e)
		elif (s < end):
			end = max(end,e)
		else:
			overlaps += [(start,end)]
			(start,end) = (s,e)

	if (start != None):
		overlaps += [(start,end)]

	return overlaps


# in_intervals--
#	Determine if a position is present in a set of intervals;  the intervals
#	are assumed to be sorted and non-overlapping
#
# Since we expect positions to arrive in order, we remember the list index
# where the most recent position was found, and begin our search at that spot
# in the list if possible.
#
# $$$ Binary search would probably be faster.

iiPrevChrom = iiPrevIx = None

def in_intervals(chrom,intervals,pos):
	global iiPrevChrom,iiPrevIx

	if (intervals == []): return False

	if (chrom != iiPrevChrom):
		startIx     = 0
		iiPrevChrom = chrom
	else:
		(s,e) = intervals[iiPrevIx]
		if (s <= pos < e): return True
		if (pos < s): startIx = 0
		else:         startIx = iiPrevIx+1

	for ix in xrange(startIx,len(intervals)):
		(s,e) = intervals[ix]
		if (pos >= e): continue
		iiPrevIx = ix
		if (pos < s): return False
		return True # if (pos < e)
	return False


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
