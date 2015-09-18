#!/usr/bin/python
#Network_umbra_stats.py
'''
A tool for counting occurrences of interactions in a
consensus network.
INPUT: '*consensus*.txt' - a tab-delimited file for the consensus network (interactions between OGs)
		See details in Network_umbra.py
OUTPUT: 'Network_statistics_out.txt', containing the following:
	Counts of total times an OG pair appears, under the following conditions:
		Full redundancy (actual PPI count across all species).
		No redundancy across species (yields count of species OG v. OG interaction is seen in)
		No redundancy across interaction type (yields count of different methods OG v. OG interaction is seen with)
		(Last option in progress)
		
Treats all interactions as undirected for consistency, so A vs. B and B vs. A count as a single unique interaction.
'''

import glob, re, sys

#Options

#Functions
def print_header(datatype): 
	#Set up the output format for the general stats
	if datatype == "network":
		stats_header = ("All interactions\tUnique OG interactions\t" + \
						"Unique Taxids\tUnique OGs\n")
	elif datatype == "interaction":
		stats_header = ("Interactors\tTotal count\tNumber of species\t" + \
						"Number of interaction types\n")
	activefile.write(stats_header)
	
def get_network_stats(network = []):
	OG_i = []
	taxid_list = []
	OG_list = []
	for ppi in network:
		ppi_taxid_pair = re.split(r' vs. ', ppi[1])
		for taxid in ppi_taxid_pair:
			if taxid not in taxid_list:
				taxid_list.append(taxid)
		this_OG_i = [ppi[0], ppi[3]]
		this_reverse_OG_i = [ppi[3], ppi[0]]
		for OG in this_OG_i:
			if OG not in OG_list:
				OG_list.append(OG)
		if this_OG_i not in OG_i:
			if this_reverse_OG_i not in OG_i:
				OG_i.append(this_OG_i)
	list_of_lists = [OG_i, taxid_list, OG_list]
	total = len(network)
	OG_i_total = len(OG_i)
	taxid_total = len(taxid_list) - 1 #Account for the input header	
	OG_list_total = len(OG_list)
	totals_list = [total, OG_i_total, taxid_total, OG_list_total]
	for item in totals_list:
		activefile.write(str(item) + "\t"),
	return list_of_lists
	
def get_interaction_stats(network = [], input_lists = []): #Relies on previous method to save time
	for interacting_pair in input_lists[0]:	#For each OG-OG pair in the list of unique OG interactors
		activefile.write(" ".join(interacting_pair) + "\t"),
		interacting_pair_count = 0
		interacting_pair_count_species = 0
		last_taxid_pair = ""
		for ppi in network:
			if interacting_pair[0] == ppi[0] and interacting_pair[1] == ppi[3]:
				interacting_pair_count = interacting_pair_count +1
				this_taxid_pair = ppi[1]
				if this_taxid_pair != last_taxid_pair:
					interacting_pair_count_species = interacting_pair_count_species +1
				last_taxid_pair = this_taxid_pair
			elif interacting_pair[0] == ppi[3] and interacting_pair[1] == ppi[0]:
				interacting_pair_count = interacting_pair_count +1
				this_taxid_pair = ppi[1]
				if this_taxid_pair != last_taxid_pair:
					interacting_pair_count_species = interacting_pair_count_species +1
				last_taxid_pair = this_taxid_pair
		activefile.write(str(interacting_pair_count) + "\t" + str(interacting_pair_count_species) + "\n")
	
#Main
if __name__ == "__main__":
	consensus_file_list = glob.glob('*consensus*.txt')
	if len(consensus_file_list) >1:
		print("One consensus network at a time, please!")
		sys.exit(0)
	if len(consensus_file_list) == 0:
		print("No consensus network file found or may not be named properly.")
		sys.exit(0)
	try:
		consensusfile = open(consensus_file_list[0])
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
		sys.exit(0)
		
	activefile_name = "Network_statistics_out.txt"
	try:
		activefile = open(activefile_name, 'w')
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
	
	print("Counting interactions in " + consensusfile.name)
	activefile.write("Statistics for " + consensusfile.name)
	consensusPPI = []
	for line in consensusfile:
		one_consensusPPI = re.split(r'\t+', line.rstrip('\t\n'))
		consensusPPI.append(one_consensusPPI) 
		
	print_header("network")
	network_lists = get_network_stats(consensusPPI)
	activefile.write("\n----- Individual interaction counts begin here. -----\n")
	print_header("interaction") 
	get_interaction_stats(consensusPPI, network_lists)
	print("Complete.")
	sys.exit(0)
