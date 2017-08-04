#!/usr/bin/env python
"""
Given a file in which each line starts with a name, select those with given names.
"""

from sys import argv,stdin,stderr,exit

def usage(s=None):
	message = """
usage: cat <data> | select_lines_by_name <names_file> [options] > <filtered>
  <names_file>       The names of items to keep from the data.  Usually this
                     will have one name per line.  This is exclusive of names
                     defined by --name= or --names=.
  --name=<string>    (cumulative) The name of one item to keep.  Cannot be used
                     with a <names_file>.
  --names=<list>     (cumulative) Names of several items to keep (a comma-
                     separated list).  Cannot be used with a <names_file>.
  --name+=<string>   (cumulative) Same as --name=, but can be used to augment
                     a <names_file>.
  --names+=<list>    (cumulative) Same as --names=, but can be used to augment
                     a <names_file>.
  --key=<column>     The column in the data, to be compared to item names.  By
                     default, column 1 is used.
  --prefix=<string>  Prefix to add to every named item.
  --suffixes=<list>  Suffixes to add to every named item (a comma-separated
                     list).
  --keepnameorder    Output data in the order of their associated names.  By
                     default, data is output in the same order as it was input.
  --nostrip          Do not strip leading and trailing blanks from data lines.
  --keepstart=<string>  Any line that starts with <string> bypasses the normal
                     process, and it copied to the output.
  --keepnonstart=<string>  Any line that does not start with <string> bypasses
                     the normal process, and it copied to the output.
  --max=<number>     Limit the number of data lines that will be output with
                     the same name.
  --caseinsensitive  consider upper/lower case differences to be matches

The basic premise of this program is that we have a list of named items that
are used to select lines from a file (on stdin)."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	namesFile     = None
	namesGiven    = []
	extraNames    = []
	keyColumn     = 0
	namePrefix    = None
	nameSuffixes  = []
	maxPerName    = None
	keepNameOrder = False
	stripLines    = True
	keepStarts    = None
	keepNonStarts = None
	caseSensitive = True

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--key=")) or (arg.startswith("--column=")):
			if ("," not in argVal):
				keyColumn = int(argVal) - 1
			else:
				keyColumn = tuple([int(col)-1 for col in argVal.split(",")])
		elif (arg == "--atstart"):   # (very slow)
			keyColumn = "atstart"
		elif (arg == "--anywhere"):  # (even slower)
			keyColumn = "anywhere"
		elif (arg.startswith("--prefix=")):
			namePrefix = argVal
		elif (arg.startswith("--suffix=")) or (arg.startswith("--suffixes=")):
			nameSuffixes += argVal.split(",")
		elif (arg.startswith("--max=")):
			maxPerName = int(argVal)
			assert (maxPerName > 0)
		elif (arg == "--keepnameorder"):
			keepNameOrder = True
		elif (arg == "--nostrip"):
			stripLines = False
		elif (arg.startswith("--name=")):
			namesGiven += [argVal]
		elif (arg.startswith("--names=")):
			arg = argVal
			newNames = arg.split(",")
			for name in newNames:
				assert (name != ""), "empty name not allowed"
			namesGiven += newNames
		elif (arg.startswith("--name+=")):
			extraNames += [argVal]
		elif (arg.startswith("--names+=")):
			arg = argVal
			newNames = arg.split(",")
			for name in newNames:
				assert (name != ""), "empty name not allowed"
			extraNames += newNames
		elif (arg.startswith("--keepstart=")):
			if (keepStarts == None): keepStarts = []
			keepStarts += [argVal]
		elif (arg.startswith("--keepnonstart=")):
			if (keepNonStarts == None): keepNonStarts = []
			keepNonStarts += [argVal]
		elif (arg == "--caseinsensitive"):
			caseSensitive = False
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (namesFile == None):
			namesFile = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (namesGiven != []) and (type(keyColumn) == tuple):
		assert (False), "--name or --names can only be used with a single key"
	if (extraNames != []) and (type(keyColumn) == tuple):
		assert (False), "--name+ or --names+ can only be used with a single key"

	if (namesGiven == []):
		assert (namesFile != None), \
		       "you have to tell me the names you're interested in"
	else:
		assert (namesFile == None), \
		       "you can't give me names in two different ways (e.g. %s and %s)" \
		     % (namesFile,namesGiven[0])

	if (keepStarts != None):
		assert (len(keepStarts) == 1), \
		       "current implementation allows no more than one keep start"
	if (keepNonStarts != None):
		assert (len(keepNonStarts) == 1), \
		       "current implementation allows no more than one keep-non start"
	if (keepStarts != None) and (keepNonStarts != None):
		assert (False), \
		       "current implementation doesn't allow both keep and keep-non starts"

	if (not caseSensitive) and (namePrefix != None):
		namePrefix = namePrefix.upper()

	if (not caseSensitive) and (nameSuffixes != []):
		nameSuffixes = [suffix.upper() for suffix in nameSuffixes]

	# read the names

	names = {}
	if (keepNameOrder): namesInOrder = []

	if (namesFile == None):
		for name in namesGiven:
			if (not caseSensitive): name = name.upper()
			if (namePrefix != None): name = namePrefix + name
			if (nameSuffixes == []):
				names[name] = 0
				if (keepNameOrder): namesInOrder += [name]
			else:
				for suffix in nameSuffixes:
					names[name+suffix] = 0
					if (keepNameOrder): namesInOrder += [name+suffix]
	else:
		namesFileKeyColumn = None
		if (":" in namesFile):
			(namesFileKeyColumn,namesFile) = namesFile.split(":",1)
			if (type(keyColumn) != tuple):
				namesFileKeyColumn = int(namesFileKeyColumn) - 1
			else:
				namesFileKeyColumns = namesFileKeyColumn.split(",")
				assert (len(namesFileKeyColumns) == len(keyColumn)), \
				       "you specify %d key columns for the data, but only %d for \"%s\"" \
				     % (len(keyColumn),len(namesFileKeyColumns),namesFile)
				namesFileKeyColumn = tuple([int(col)-1 for col in namesFileKeyColumns])

		if (namesFileKeyColumn == None) and (type(keyColumn) == tuple):
			namesFileKeyColumn = keyColumn

		f = file(namesFile,"rt")

		lineNumber = 0
		for line in f:
			lineNumber += 1
			line = line.strip()
			if (line.startswith("#")): continue

			if (namesFileKeyColumn == None):
				pass
			elif (type(namesFileKeyColumn) != tuple):
				line = line.split()[namesFileKeyColumn]
			else:
				fields = line.split()
				assert (len(fields) > max(namesFileKeyColumn)), \
				       "you specify %d key columns for \"%s\", but line %d has only %d" \
				     % (len(namesFileKeyColumn),namesFile,lineNumber,len(fields))
				line = tuple([fields[col] for col in namesFileKeyColumn])

			if (not caseSensitive):
				if (type(line) != tuple):
					line = line.upper()
				else:
					line = tuple([field.upper() for field in line])

			if (namePrefix != None):
				if (type(line) != tuple):
					line = namePrefix + line
				else:
					line = tuple([namePrefix+field for field in line])

			if (nameSuffixes == []):
				names[line] = 0
				if (keepNameOrder): namesInOrder += [line]
			else:
				for suffix in nameSuffixes:
					if (type(line) != tuple):
						linewithSuffix = line+suffix
					else:
						linewithSuffix = tuple([field+suffix for field in line])
					names[linewithSuffix] = 0
					if (keepNameOrder): namesInOrder += [linewithSuffix]

		f.close()

	for name in extraNames:
		if (namePrefix != None): name = namePrefix + name
		if (nameSuffixes == []):
			names[name] = 0
			if (keepNameOrder): namesInOrder += [name]
		else:
			for suffix in nameSuffixes:
				names[name+suffix] = 0
				if (keepNameOrder): namesInOrder += [name+suffix]

	# process the lines

	nameToLines = {}

	lineNumber = 0
	for line in stdin:
		lineNumber += 1
		if (stripLines): line = line.strip()
		else:            line = line.rstrip("\n")

		if (keepStarts != None) and (line.startswith(keepStarts[0])):
			print line
			continue

		if (keepNonStarts != None) and (not line.startswith(keepNonStarts[0])):
			print line
			continue

		if (line.startswith("#")): continue

		if (type(keyColumn) == int):
			fields = line.split()
			if (len(fields) <= keyColumn): continue
			key = fields[keyColumn]
			if (not caseSensitive): key = key.upper()
			if (key not in names): continue
			name = fields[keyColumn]
		elif (type(keyColumn) == tuple):
			fields = line.split()
			if (len(fields) <= max(keyColumn)): continue
			key = []
			for col in keyColumn:
				k = fields[col]
				if (not caseSensitive): k = k.upper()
				key += [k]
			key = tuple(key)
			if (key not in names): continue
			name = tuple([fields[col] for col in keyColumn])
		elif (keyColumn == "atstart"):
			found = False
			for name in names:
				key = line
				if (not caseSensitive): key = key.upper()
				if (key.startswith(name)):
					found = True
					break
			if (not found): continue
		elif (keyColumn == "anywhere"):
			found = False
			for name in names:
				key = line
				if (not caseSensitive): key = key.upper()
				if (name in key):
					found = True
					break
			if (not found): continue

		if (keepNameOrder):
			if (line not in nameToLines): nameToLines[name] =  [line]
			else:                         nameToLines[name] += [line]
		else:
			print line

		if (maxPerName != None):
			if (names[name] < maxPerName): names[name] += 1
			else:                          del names[name]

	if (keepNameOrder):
		for name in namesInOrder:
			if (name not in nameToLines): continue
			print "\n".join(nameToLines[name])


if __name__ == "__main__": main()
