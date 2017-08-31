#!/usr/bin/env python
"""
Use one set of genomic intervals to 'clip' another.
"""

from sys  import argv,stdin,stderr,exit
from gzip import open as gzip_open


def usage(s=None):
	message = """
usage: cat intervals_file | clip_intervals <clip_file> [options]
  <clip_file>          file of intervals to clip to
  --dilation=<length>  number of bases to add to each side of a clip interval
  --sorted             input intervals have been sorted along each chromosome;
                       this saves a lot of memory if the intervals file is
                       large
  --chrom[s]=<list>    limit intervals to specific chromosomes (comma-separated)
  --ignore=<prefix>    (cumulative) ignore any interval lines that begin with
                       this prefix
  --origin=one         intervals are origin-one, closed
  --origin=zero        intervals are origin-zero, half-open
                       (this is the default)
  --progress=<count>   report status every <count> bases

Use one set of genomic intervals to 'clip' another.

This is essentially the same as intersection, performed independently between
each interval and the set of clipping intervals.  Thus, the input may contain
overlapping intervals.  The output can contain overlapping intervals (unless
the --sorted option is used).  Further, all additional information on the lines
is kept in the clipped output."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global ignorePrefixes,chromsOfInterest,origin
	global reportProgress
	global debug

	# parse args

	clipFilename     = None
	dilation         = None
	inputIsSorted    = False
	chromsOfInterest = None
	ignorePrefixes   = None
	origin           = "zero"
	reportProgress   = None
	debug            = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--dilation=")) or (arg.startswith("--dilate=")):
			dilation = int_with_unit(argVal)
		elif (arg == "--sorted"):
			inputIsSorted = True
		elif (arg.startswith("--chrom=")) or (arg.startswith("--chroms=")):
			if (chromsOfInterest == None): chromsOfInterest = []
			chromsOfInterest += argVal.split(",")
		elif (arg.startswith("--origin=")):
			origin = argVal
			if (origin == "0"): origin = "zero"
			if (origin == "1"): origin = "one"
			assert (origin in ["zero","one"]), "can't understand %s" % arg
		elif (arg.startswith("--ignore=")):
			if (ignorePrefixes == None): ignorePrefixes = []
			ignorePrefixes += [argVal]
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (clipFilename == None):
			clipFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (clipFilename == None):
		usage("you must tell me where to find the clipping intervals")

	# read the clipping intervals

	if (reportProgress != None):
		print >>stderr, "progress: reading clipping intervals"

	if (clipFilename.endswith(".gz")) or (clipFilename.endswith(".gzip")):
		f = gzip_open(clipFilename,"rt")
	else:
		f = file(clipFilename,"rt")

	chromToClippers = {}
	for (_,chrom,start,end) in read_intervals(f,keepExtra=False):
		if (chromsOfInterest != None) and (chrom not in chromsOfInterest):
			continue
		if (chrom not in chromToClippers):
			chromToClippers[chrom] = []
		chromToClippers[chrom] += [(start,end)]

	f.close()

	if (dilation != None):
		for chrom in chromToClippers:
			clippers = []
			for (s,e) in chromToClippers[chrom]:
				s -= dilation
				e += dilation
				if (s < 0): s = 0
				clippers += [(s,e)]
				# (note that we don't need to worry here about exceeding the
				#  .. upper bound on a chromosome length;  it's OK if we do)
			chromToClippers[chrom] = clippers

	for chrom in chromToClippers:
		chromToClippers[chrom] = non_overlapping_intervals(chromToClippers[chrom])

	# clip the intervals

	if (inputIsSorted): clip_one_by_one (chromToClippers)
	else:               collect_and_clip(chromToClippers)


# read the intervals into memory, then clip them

def collect_and_clip(chromToClippers):

	# collect the intervals to be clipped

	if (reportProgress != None):
		print >>stderr, "progress: reading unsorted intervals"

	chromToIntervals = {}
	chroms = []

	for (_,chrom,start,end,extra) in read_intervals(stdin):
		if (chromsOfInterest != None) and (chrom not in chromsOfInterest):
			continue
		if (chrom not in chromToIntervals):
			chromToIntervals[chrom] = []
			chroms += [chrom]
		chromToIntervals[chrom] += [(start,end,extra)]

	for chrom in chromToIntervals:
		chromToIntervals[chrom].sort()

	# clip the intervals

	for chrom in chroms:
		if ("chroms" in debug): print >>stderr, chrom
		prevBlock = None

		if (chrom not in chromToClippers): continue
		clippers = chromToClippers[chrom]
		numClippers = len(clippers)
		if (numClippers == 0): continue

		intervals = chromToIntervals[chrom]

		clipIx = 0
		(clipStart,clipEnd) = clippers[clipIx]
		for (start,end,extra) in intervals:
			# advance to the first clipper that does not precede this interval;
			# if there is no such clipper, then all the remaining intervals are
			# also past any clipper, and we can discard them all

			if (reportProgress != None):
				if (start/reportProgress != prevBlock):
					print >>stderr, "progress: %s %s" % (chrom,commatize(start))
					prevBlock = start / reportProgress

			if ("intervals" in debug):
				print >>stderr, chrom,start,end,extra

			while (clipEnd <= start):
				clipIx += 1
				if (clipIx >= numClippers): break
				(clipStart,clipEnd) = clippers[clipIx]

			if (clipIx >= numClippers):
				break

			# if this interval is before the clipper, discard this interval

			if (end <= clipStart): continue

			# scan through subsequent clippers and clip this interval to them

			(ix,cs,ce) = (clipIx,clipStart,clipEnd)

			while (cs < end):
				s = max(cs,start)
				e = min(ce,end)
				if (extra == []): print "%s %d %d"    % (chrom,s,e)
				else:             print "%s %d %d %s" % (chrom,s,e," ".join(extra))

				ix += 1
				if (ix >= numClippers): break
				(cs,ce) = clippers[ix]


# read each interval into memory, clip it, then read the next

def clip_one_by_one(chromToClippers):
	chromSeen = {}
	prevChrom = prevBlock = None
	for (lineNumber,chrom,start,end,extra) in read_intervals(stdin):

		if (reportProgress != None):
			if (chrom != prevChrom) or (start/reportProgress != prevBlock):
				print >>stderr, "progress: %s %s" % (chrom,commatize(start))
				prevBlock = start / reportProgress

		# if this is a different chromosome, make sure we haven't see it
		# before, and set up the clipping list

		if (chrom != prevChrom):
			assert (chrom not in chromSeen), \
			       "intervals for %s are not together (lines %d and %d)" \
			     % (chrom,chromSeen[chrom],lineNumber)
			chromSeen[chrom] = lineNumber
			prevChrom = chrom

			clippers = None
			if (chrom in chromToClippers):
				clippers = chromToClippers[chrom]
				numClippers = len(clippers)
				if (numClippers == 0):
					clippers = None
				else:
					clipIx = 0
					(clipStart,clipEnd) = clippers[clipIx]
	
			if (clippers != None):
				if ("chroms" in debug): print >>stderr, chrom
			else:
				if ("chroms" in debug): print >>stderr, "%s (discarded)" % chrom
				continue

		# if we're on the same chromosome but there were no clipping intervals,
		# skip this interval

		elif (clippers == None):
			continue

		# otherwise, we're still on the same chromosome, so make sure this
		# interval is in the proper order

		else:
			assert (start >= prevEnd), \
			       "intervals out of order at line %d: %s %d %d then %s %d %d" \
			     % (lineNumber,chrom,prevStart,prevEnd,chrom,start,end)

		# remember this as the most recent interval processed

		(prevStart,prevEnd) = (start,end)

		# advance to the first clipper that does not precede this interval;
		# if there is no such clipper, then all the remaining intervals are
		# also past any clipper, and we can discard them all

		if ("intervals" in debug):
			print >>stderr, chrom,start,end,extra

		while (clipEnd <= start):
			clipIx += 1
			if (clipIx >= numClippers): break
			(clipStart,clipEnd) = clippers[clipIx]

		if (clipIx >= numClippers):
			clippers = None
			continue

		# if this interval is before the clipper, discard this interval

		if (end <= clipStart): continue

		# scan through subsequent clippers and clip this interval to them

		(ix,cs,ce) = (clipIx,clipStart,clipEnd)

		while (cs < end):
			s = max(cs,start)
			e = min(ce,end)
			if (extra == []): print "%s %d %d"    % (chrom,s,e)
			else:             print "%s %d %d %s" % (chrom,s,e," ".join(extra))

			ix += 1
			if (ix >= numClippers): break
			(cs,ce) = clippers[ix]


def read_intervals(f,keepExtra=True):
	numFields = None

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		if (ignorePrefixes != None):
			ignoreMe = False
			for prefix in ignorePrefixes:
				if (not line.startswith(prefix)): continue
				ignoreMe = True
				break
			if (ignoreMe): continue

		fields = line.split()
		if (numFields == None):
			assert (len(fields) >= 3), \
			       "not enough fields at line %d (%d, expected at least %d)" \
			     % (lineNumber,len(fields),3)
			numFields = len(fields)
		else:
			assert (len(fields) == numFields), \
			       "inconsistent number of fields at line %d (%d, expected %d)" \
			     % (lineNumber,len(fields),numFields)

		try:
			chrom = fields[0]
			start = int(fields[1])
			end   = int(fields[2])
			if (keepExtra): extra = fields[3:]
			if (end < start): raise ValueError
			if (origin == "one"): start -= 1
		except ValueError:
			assert (False), "bad line (%d): %s" % (lineNumber,line)

		if (keepExtra): yield (lineNumber,chrom,start,end,extra)
		else:           yield (lineNumber,chrom,start,end)


# merge into non-overlapping intervals

def non_overlapping_intervals(intervals):
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
