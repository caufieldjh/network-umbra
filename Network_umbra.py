#!/usr/bin/python
#Network_umbra.py
#Predicts interactions in a protein interaction network based off a consensus network.
#REQUIRES: Biopython 1.65 or more recent
#INPUT: 'consensus.sif' - a SIF file for the consensus network (interactions between OGs)
#			Each line includes an OG, a species ID, and an interacting OG.
#			Species IDs are names now but will be NCBI taxon IDs.
#			An interaction may appear more than once in the consensus network.
#		'target.txt' - a list of proteins and corresponding OGs present in the target genome/proteome.
#			This list should not contain duplicate proteins but can contain duplicate OGs.
#			Each line contains an OG first, a tab, and a unique protein-coding locus.
#			An alternate filename may be specified as the first ARGV parameter.
#		The second ARGV parameter is the NCBI taxon ID corresponding to the target.
#OUTPUT: 'predicted_interactions_[taxid].txt' - a file of the predicted interactions for the target species
#			Each line includes an OG, a type of prediction,
#			an interacting OG, and the taxonomy ID of the original interaction source (usually a species).
#			Interologs are evaluated by similarity of taxonomic lineage (not great but informative);
#			Higher values indicate greater similarity to the target species and potentially
#			better predictions or just broadly-conserved interactions.
#		'noninteracting_OGs_[taxid].txt' - file containing OGs and loci from the target list not found 
#			in any predicted interacions. They may still be in the consensus network somewhere.
#		'umbra_taxid_db' - storage for similarities between genomes corresponding to taxon IDs.
#			Created if it does not yet exist and appended if new (across all local sessions) taxon IDs are used.
#
#At the moment, this only makes predictions based off presence of the same OGs as in the consensus network.
#It needs to verify that both OGs in the predicted PPI are present in the target species.
#Redundant predictions (the same interaction from the same taxon ID) are merged.
#
#PLANS:
#This pipeline is really two main elements: 
#	generation of a consensus network and
#	given a target proteome, generation of a new network from that set.
#*Get counts of matched proteins, OGs, and pairs, and the same for all unmatched.
#*Set up to be batch-run for many proteomes.
#*Incorporate pair odds as per Rodgers-Melnick et al. 2013 BMC Genomics (ENTS PPI network prediction)
#	These values need to be calculated per-species...or do they?
#	Should we instead use a log-likelihood function and just find the "unique" interactions?
#*Parse the IntAct output directly or get it using PSICQUIC (not ideal for grabbing many PPI but good for details).
#*Couple predicted interactions back to specific protein-coding genes in the target genome.
#*Weight interactions by:
#	Number of species seen in
#	Distance between observed species (requires distance calculations and need to use unique species IDs)
#		Will work on using Biopython modules for this
#		Could just use existing NCBI taxonomy as a proxy but will have to pre-download as NCBI 
#			doesn't like large queries. Do so at runtime and save locally for future runs.
#		Can also handle with MUMer v.3 as per Deloger et al. 2008 J Bac
#			(should still only compare genomes once and save to local database)
#			But how to get it to play nice with Python?
#		There is also the GGDC method - http://ggdc.dsmz.de/ - but it only has an HTML interface
#	Distance between proteins in OGs (requres OG comparisons)
#	Essentiality in the target species (requires essentiality data)
#*Include partial matches so "missing" interactors could be identified

import sys, re
from Bio import Entrez
Entrez.email = 'caufieldjh@vcu.edu'

#Options
taxonomy_compare = 0
       
#Load consensus network file
try:
	consensusfile = open("consensus.sif")
except IOError as e:
	print("I/O error({0}): {1}".format(e.errno, e.strerror))	
consensusPPI = []
for line in consensusfile:
	#print line
	one_consensusPPI = re.split(r'\t+', line.rstrip('\t\n'))
	consensusPPI.append(one_consensusPPI)

#Load target file as default or as stated in argv	
if (len(sys.argv)>1):
	#print str(sys.argv[1])
	targetfile = open(str(sys.argv[1]))
	target_taxid = int(sys.argv[2])
else:
	try:
		targetfile = open("target.txt")	#The test file is from H. pylori 26695, taxid 85962
		target_taxid = 85962
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
		
#Store target proteins and OGs
#Get information about target from Entrez
target_loci = [] #Both OG and corresponding unique locus
target_ogs = [] #Each OG will be unique in this list
target_proteins = [] #Each locus, really
for line in targetfile:
	og_and_prot = re.split(r'\t+', line.rstrip('\t\n'))					
	target_loci.append(og_and_prot)
	if og_and_prot[0] not in target_ogs:
		target_ogs.append(og_and_prot[0])
	target_proteins.append(og_and_prot[0])
target_handle = Entrez.efetch(db="Taxonomy", id=str(target_taxid), retmode="xml")
target_records = Entrez.read(target_handle)
target_name = target_records[0]["ScientificName"]
target_lineage = target_records[0]["Lineage"]
print("Predicting interactions for " + target_name +"\n")
print("Unique OGs\tUnique proteins")
print("%s\t\t%s") % (len(target_ogs), len(target_proteins))
activefile_name = "Predicted_interactions_" + str(target_taxid) + ".txt"
try:
	activefile = open(activefile_name, 'w')
except IOError as e:
	print("I/O error({0}): {1}".format(e.errno, e.strerror))
noninteracting_file_name = "noninteracting_OGs_" + str(target_taxid) + ".txt"
try: 
	noninteracting_file = open(noninteracting_file_name, 'w')
except IOError as e:
	print("I/O error({0}): {1}".format(e.errno, e.strerror))
	
#Option - use taxids to build set of similarities
#Just gets NCBI Taxonomy lineage for now
#Currently producing errors
if taxonomy_compare == 1:
	taxid_unique = []
	taxid_new = []
	try:
		taxid_db = open('umbra_taxid_db', 'r+')
		for line in taxid_db:
			split_line = re.split(r'\t+', line.rstrip('\t\n'))
			if split_line[0] not in taxid_unique:
				taxid_unique.append(split_line[0])
	except IOError as e:
		print("Could not open similarity database or file not found.\nCreating new file.")
		taxid_db = open('umbra_taxid_db', 'w')
		print("Please wait...")
	if target_taxid not in taxid_unique:
		taxid_new.append(target_taxid)
	for ppi in consensusPPI:
		if ppi[1] not in taxid_unique:
			taxid_new.append(ppi[1])
		for taxid in taxid_new:
			target_handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
			target_records = Entrez.read(target_handle)
			target_name = target_records[0]["Lineage"]
			print(target_name)
			taxid_db.write(str(taxid) + "\t" + target_name + "\n")

#Set up predicted interaction network and search for presence of interactors
predicted_net = []
match_count = 0
for og in target_ogs:
	for ppi in consensusPPI:
		if og == ppi[0] and ppi[2] in target_ogs:
			predicted_net.append(ppi)
			match_count = match_count +1

#Remove duplicates predictions per species
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

#Send to output
for ppi in predicted_net:			
	if int(ppi[1]) == target_taxid:
		method = "Experimental results"
		experimental_count = experimental_count +1
	else:
		if taxonomy_compare == 1:				#Determine taxonomic lineage-based similarity if asked
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
	out_string = (str(ppi[0]) + "\t" + method + "\t" + str(ppi[2]) + "\t" + str(ppi[1]) + "\n")
	activefile.write(out_string)
for og_and_prot in target_loci:
	if og_and_prot[0] not in predicted_OG_coverage_unique:
		out_string = (str(og_and_prot[0]) + "\t" + str(og_and_prot[1] + "\n"))
		noninteracting_file.write(out_string)
	
activefile.close()
noninteracting_file.close()

print("Found " + str(match_count) + " interaction predictions (redundant by species).")
print("Of these, at least " + str(experimental_count) + " are from experimental data for this species.")
#it's "at least" because different studies use different taxids for the same species (i.e. different strains)
print(str(len(predicted_net_unique_alltaxid)) + " interactions are in the predicted network.")
print("These interactions include " + str(len(predicted_OG_coverage_unique)) + " unique OGs.")
print("See the results in " + activefile.name)
sys.exit(0)
