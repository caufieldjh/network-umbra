#!/usr/bin/python
#Network_umbra.py
'''
Predicts interactions in a protein interaction network based off a consensus network.
REQUIRES: Biopython 1.65 or more recent
INPUT: '*consensus*.txt' - a tab-delimited file for the consensus network (interactions between OGs)
			Each line includes five strings, in this order:
			an OG, 
			a species ID or a pair of IDs (format, "197 vs. 192222")
			a species name or a pair of species names, ideally with a strain identifier, 
			an interacting OG,
			and an interaction type.
			Species IDs are NCBI taxon IDs.
			Interaction types are those defined by the PSI-MI vocabulary.  
			An interaction may appear more than once in the consensus network.
		'TAXID-target.txt' - a list of proteins and corresponding OGs present in the target genome/proteome.
			This list should not contain duplicate proteins but can contain duplicate OGs.
			Each line contains an OG first, a tab, and a unique protein-coding locus.
			The filename must have the format "TAXID-target.txt"
			An alternate filename may be specified as the first ARGV parameter.
			Or, use batch mode and it will open usable files with the above format.
OUTPUT: 'predicted_interactions_[taxid].txt' - a file of the predicted interactions for the target species
			Each line resembles those in the consensus network in content.
			An additional identifier is included to describe the type of prediction:
				Experimental results - this PPI has been found in results using this specific species and strain
				Experimental results, spoke expansion - as above but originally from a spoke expansion model
				Experimental results, related strain - this PPI has been found in results using the same species but different taxid
				Predicted interolog - this PPI has been found in results using a different species.
		'noninteracting_OGs_[taxid].txt' - file containing OGs and loci from the target list not found 
			in any predicted interacions. They may still be in the consensus network somewhere.
		'umbra_taxid_db' - storage for similarities between genomes corresponding to taxon IDs.
			Created if it does not yet exist and appended if new (across all local sessions) taxon IDs are used.

At the moment, this only makes predictions based off presence of the same OGs as in the consensus network.
It needs to verify that both OGs in the predicted PPI are present in the target species.
Redundant predictions (the same interaction from the same taxon ID) are merged.

PLANS:
This pipeline is really two main elements: 
	generation of a consensus network and
	given a target proteome, generation of a new network from that set.
*Get counts of matched proteins, OGs, and pairs, and the same for all unmatched.
	To be biologically meaningful, this should include protein-coding loci w/o COG annotations.
	PPI matches should be counted as merged by species, too.
*Incorporate pair odds as per Rodgers-Melnick et al. 2013 BMC Genomics (ENTS PPI network prediction)
	These values need to be calculated per-species...or do they?
	Should we instead use a log-likelihood function and just find the "unique" interactions?
*Parse the IntAct output directly or get it using PSICQUIC (not ideal for grabbing many PPI but good for details).
*Couple predicted interactions back to specific protein-coding genes in the target genome.
*Weight interactions by:
	Number of species seen in
	Distance between observed species (requires distance calculations and need to use unique species IDs)
		Will work on using Biopython modules for this
		Could just use existing NCBI taxonomy as a proxy but will have to pre-download as NCBI 
			doesn't like large queries. Do so at runtime and save locally for future runs.
		Can also handle with MUMer v.3 as per Deloger et al. 2008 J Bac
			(should still only compare genomes once and save to local database)
			But how to get it to play nice with Python?
		There is also the GGDC method - http://ggdc.dsmz.de/ - but it only has an HTML interface
	Distance between proteins in OGs (requres OG comparisons)
	Essentiality in the target species (requires essentiality data)
'''

import sys, re, glob, os
from Bio import Entrez
Entrez.email = 'caufieldjh@vcu.edu'

#Options
taxonomy_compare = 0
batch_mode = 1

#Functions

def print_header(): 
	#Set up the output format for the general stats
	stats_header = ("Name\tUnique proteins\tProteins without OG\tProteins Not in PPI\tUnique OGs\tTotal predicted PPI\tExperimental PPI\tUnique PPI in Predicted Network\tUnique OGs in Predicted Network\n")
	print(stats_header)

def get_lineage(taxid):
	#Get information about target from Entrez, but mostly name and parent taxid
	target_handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
	target_records = Entrez.read(target_handle)
	name = target_records[0]["ScientificName"]
	lineage = target_records[0]["Lineage"] #Don't use this yet but might need it
	parent = target_records[0]["ParentTaxId"]
	return [name, lineage, parent]
	
def network_store(target, target_taxid):	
	#Store target proteins and OGs
	target_loci = [] #Both OG and corresponding unique locus
	target_ogs = [] #Each OG will be unique in this list
	target_proteins = [] #Each locus, really
	proteins_are_na = 0
	for line in target:
		og_and_prot = re.split(r'\t+', line.rstrip('\t\n'))					
		target_loci.append(og_and_prot)
		if og_and_prot[0] == "NA":
			proteins_are_na = proteins_are_na + 1
		if og_and_prot[0] not in target_ogs:
			target_ogs.append(og_and_prot[0])
		target_proteins.append(og_and_prot[0])
	[target_name, target_lineage, parent_taxid] = get_lineage(target_taxid)
	activefile_name = "Predicted_interactions_" + str(target_taxid) + ".txt"
	try:
		activefile = open(activefile_name, 'w')
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
	noninteracting_file_name = "Noninteracting_OGs_" + str(target_taxid) + ".txt"
	try: 
		noninteracting_file = open(noninteracting_file_name, 'w')
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
	#Set up predicted interaction network and search for presence of interactors
	predicted_net = []
	match_count = 0
	for og in target_ogs:
		for ppi in consensusPPI:
			if og == ppi[0] and ppi[3] in target_ogs:
				predicted_net.append(ppi)
				match_count = match_count +1
	#Remove duplicate predictions per species
	predicted_net_unique = []
	predicted_net_unique_alltaxid = []
	predicted_OG_coverage = []
	predicted_OG_coverage_unique = []
	experimental_count = 0 #The number of PPI already found for this taxid
	for ppi in predicted_net:
	    if ppi not in predicted_net_unique:
	        predicted_net_unique.append(ppi)
		this_prediction = [ppi[0], ppi[2]]
		if this_prediction not in predicted_net_unique_alltaxid:
			predicted_net_unique_alltaxid.append(this_prediction)
	predicted_net = predicted_net_unique
	predicted_OG_coverage = [x for y in predicted_net_unique_alltaxid for x in y]
	for og in predicted_OG_coverage:
		if og not in predicted_OG_coverage_unique:
			predicted_OG_coverage_unique.append(og)
	#Send to output and classify each interaction as experimental or predicted
	for ppi in predicted_net:
		ppi_taxid_pair = re.split(r' vs. ', ppi[1])			
		if (target_taxid in ppi_taxid_pair or parent_taxid in ppi_taxid_pair):
			if target_taxid == ppi_taxid_pair[0] and target_taxid == ppi_taxid_pair[1]:
				if ppi[4] == "association":
					method = "Experimental results, spoke expansion"
				else:	
					method = "Experimental results"
			else:
				method = "Experimental results, related strain"
			experimental_count = experimental_count +1	
		else:
			if taxonomy_compare == 1:	#Determine taxonomic lineage-based similarity if asked
				taxid_db.seek(0,0)
				level_match_score = 0
				target_lineage_line = re.split(r';+', target_lineage)
				for line in taxid_db:
					split_line = re.split(r'\t+', line.rstrip('\t\n'))
					if int(ppi[1]) == int(split_line[0]):
						lineage_line = re.split(r';+', split_line[1])
						for level in target_lineage_line:
							if level in lineage_line:
								level_match_score = level_match_score +1
				method = "Predicted interolog, level " + str(level_match_score)
			else:
				method = "Predicted interolog"
		out_string = ("\t".join(ppi) + "\t" + method + "\n")
		activefile.write(out_string)
	activefile.close()
	non_interactor_count = 0
	for og_and_prot in target_loci:
		if og_and_prot[0] not in predicted_OG_coverage_unique:
			out_string = ("\t".join(og_and_prot) + "\n")
			noninteracting_file.write(out_string)
			non_interactor_count = non_interactor_count + 1
	noninteracting_file.close()
	stats_output = [target_name, str(len(target_proteins)), str(proteins_are_na), str(non_interactor_count), str(len(target_ogs)), str(match_count), str(experimental_count), str(len(predicted_net_unique_alltaxid)), str(len(predicted_OG_coverage_unique))]
	print("\t".join(stats_output) + "\n")
	
#Load consensus network file
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
print("Using " + consensusfile.name + " as the consensus network.")
consensusPPI = []
for line in consensusfile:
	one_consensusPPI = re.split(r'\t+', line.rstrip('\t\n'))
	consensusPPI.append(one_consensusPPI)
	
print_header()

#Load target file as default or as stated in argv
#Run the network_store method to do the actual work
if batch_mode == 0:	
	if (len(sys.argv)>1):
		try:
			targetfile = open(str(sys.argv[1]))
		except IOError as e:
			print("I/O error({0}): {1}".format(e.errno, e.strerror))
		target_taxid = (re.split('-', sys.argv[1]))[0]
		network_store(targetfile, target_taxid)
	else:
		print("No target specified.")
		sys.exit(0)	
else:
	target_file_list = glob.glob('*target.txt')
	if len(target_file_list) == 0:
		print("No target proteome files found, or may not be named properly.")
		sys.exit(0)
	for filename in target_file_list:
		taxid = (re.split('-', filename))[0]
		targetfile = open(filename)
		network_store(targetfile, taxid)
		targetfile.close()
		
#Option - use taxids to build set of similarities
#Just gets NCBI Taxonomy lineage for now
#Currently producing errors
#if taxonomy_compare == 1:
	#taxid_unique = []
	#taxid_new = []
	#try:
		#taxid_db = open('umbra_taxid_db', 'r+')
		#for line in taxid_db:
			#split_line = re.split(r'\t+', line.rstrip('\t\n'))
			#if split_line[0] not in taxid_unique:
				#taxid_unique.append(split_line[0])
	#except IOError as e:
		#print("Could not open similarity database or file not found.\nCreating new file.")
		#taxid_db = open('umbra_taxid_db', 'w')
		#print("Please wait...")
	#if target_taxid not in taxid_unique:
		#taxid_new.append(target_taxid)
	#for ppi in consensusPPI:
		#if ppi[1] not in taxid_unique:
			#taxid_new.append(ppi[1])
		#for taxid in taxid_new:
			#target_handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
			#target_records = Entrez.read(target_handle)
			#target_name = target_records[0]["Lineage"]
			#print(target_name)
			#taxid_db.write(str(taxid) + "\t" + target_name + "\n")

sys.exit(0)
