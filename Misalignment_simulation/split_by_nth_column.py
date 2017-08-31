#!/usr/bin/env python
"""
Split a file into multiple files, based on the Nth column in each line
----------------------------------------------------------------------

The input alignment (on stdin) is a white-space-delimited file.  The name_format
must contains the substring "{name}", which will be replaced by the contents of
the specified column.

usage:  cat farf | split_by_nth_column --column=<number> name_format

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys,gzip


def main():

	# parse args

	keyColumn  = 0
	removeKey  = False
	keepHeader = False
	separator  = "\t"
	fileMode   = "wt"
	nameFmt    = None

	for arg in sys.argv[1:]:
		if (arg.startswith("--key=")) or (arg.startswith("--column=")):
			keyColumn = int(arg.split("=",1)[1]) - 1
		elif (arg == "--removekey"):
			removeKey = True
		elif (arg == "--keepheader"):
			keepHeader = True
		elif (arg == "--spaces"):
			separator = " "
		elif (arg == "--newlines"):
			separator = "\n"
		elif (arg == "--append"):
			fileMode = "at"
		elif (arg.startswith("--")):
			assert (False), "unknown argument: %s" % arg
		elif (nameFmt == None):
			nameFmt = arg
		else:
			assert (False), "can't understand %s" % arg

	assert ("{name}" in nameFmt), "bad format string: %s" % nameFmt

	# process the lines from the file

	nameToFile = {}
	header = None

	lineNum = 0
	for line in sys.stdin:
		lineNum += 1
		line = line.rstrip()
		if (line.startswith("#")):
			if (keepHeader):
				if (header == None): header = line
				else:                header = header + "\n" + line
			continue
		if (line == ""): continue

		fields = line.split()
		if (len(fields) <= keyColumn): continue

		name = fields[keyColumn]
		if (name not in nameToFile):
			filename = nameFmt.replace("{name}",name)
			if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
				f = gzip.open(filename,fileMode)
			else:
				f = file(filename,fileMode)
			nameToFile[name] = f
			if (header != None): print >>f, header
		else:
			f = nameToFile[name]

		if (removeKey):
			fields = fields[:keyColumn] + fields[keyColumn+1:]
			line   = separator.join(fields)

		print >>f, line

	for name in nameToFile:
		f = nameToFile[name]
		f.close()


if __name__ == "__main__": main()
