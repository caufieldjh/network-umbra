#!/usr/bin/python
#Network_umbra_protein_core.py
'''
A tool for determining a core proteome set from 
a given set of proteomes with orthology annotations.
INPUT: 
	'TAXID-target.txt' - one or more tab-delimeted files containing, 
	on each line, an OG and a unique protein-coding locus.
	The TAXID is a unique NCBI taxonomy ID.
OUTPUT: 
	A list of unique OGs and their occurrence count across all input files.
	'core-proteome.txt' - a single tab-delimited file containing, 
	on each line, an OG and a list of taxids containing the OG.
	This file represents the core proteome as determined from the input.

The threshold for inclusion in the core proteome is determined by
presence of an OG in more than one input proteome file.
'''

import glob, re, sys

#Options
threshold = 5	#The number of proteomes an OG must be in to be included

#Functions
def setup_core():
	all_ogs = []
	core_ogs = []
	core_ogs_unique = []
	target_file_list = glob.glob('*target.txt')
	if len(target_file_list) == 0:
		print("No target proteome files found, or may not be named properly.")
		sys.exit(0)
	for filename in target_file_list:
		this_ogs = []
		this_ogs_unique = []
		taxid = (re.split('-', filename))[0]
		targetfile = open(filename)
		for line in targetfile:
			target_og = (re.split('\t', line))[0]
			this_ogs.append(target_og)
		targetfile.close()
		for og in this_ogs:
			if og not in this_ogs_unique:
				this_ogs_unique.append(og)
		for og in this_ogs_unique:
			all_ogs.append(og)		#If an og appears in this list more than once, it's in more than one proteome
	for og in all_ogs:
		this_count = all_ogs.count(og)
		if this_count >= threshold:
			core_ogs.append([og, this_count])
	for og_and_count in core_ogs:
		if og_and_count not in core_ogs_unique:
			core_ogs_unique.append(og_and_count)
	return core_ogs_unique
	
#Main
if __name__ == "__main__":
	print("Establishing core proteome.")
	print("Threshold is " + str(threshold) + " proteomes.")
	core = setup_core()
	try:
		core_file = open("core-proteome.txt", 'w')
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
		sys.exit(0)
	for og_and_count in core:
		core_file.write(("\t".join(str(x) for x in og_and_count) + "\n"))
	sys.exit(0)
