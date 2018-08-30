import sys
import xml.etree.ElementTree as ET

######################################
#python python extract_from_xml.py metadata.xml
#Extracts the chemistry type from pacbio metadata.xml metadata_file
######################################
metadata_file=sys.argv[1]

tree = ET.parse(metadata_file)
root = tree.getroot()

for child in root:
	if ('BindingKit' in child.tag):
		for prop in child:
			if ('Name' in prop.tag):
				print(metadata_file + " " + prop.text)

