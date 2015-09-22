#!/usr/bin/python
#Network_umbra.py
'''
Predicts interactions in a protein interaction network based off a meta-interactome network.
Uses eggNOG v.4.1.

REQUIRES: Biopython 1.65 or more recent
			Also needs ~600 MB of available disk space to accomodate data files and output
			More space may be necessary for proteome files.

INPUT: Downloads all available protein-protein interactions for bacteria from IntAct.
		Alternatively, uses a provided PPI data file in PSI-MI TAB 2.7 format.
		REMOVE THE HEADER ROW if it's present!
		Downloads highest-level (LUCA) and bacteria-specific Uniprot ID to NOG mappings from eggNOG v.4.1.
		Downloads highest-level (LUCA) bacteria-specific NOG annotations from eggNOG v.4.1.

OUTPUT: 
'metainteractome[date].txt'
			A meta-interactome composed of all available bacterial protein-protein interactions.
			Follows PSI-MI Tab27 format, with the addition of two ortholog identifiers per row.
			See format description at https://code.google.com/p/psimi/wiki/PsimiTab27Format
			
'meta_statistcs[date].txt'
			Contains statistics relevant to the produced meta-interactome.
			
'consensus[date].txt'
			A consensus meta-interactome composed of all available bacterial protein-protein interactions.
			This set of interactions compresses all unique proteins into their corresponding orthologous groups.
			Data in each column is the following, from left to right:
InteractorA		The first interactor. Usually an OG.
InteractorB		The second interactor. Usually an OG.
InteractionCount		Count of individual PROTEIN interactions contributing to this consensus interaction, as per the meta-interactome.
TaxonCount		Count of different taxons (here, a proxy for species) corresponding to the interaction.
			Similar taxons have been grouped together where possible, e.g. two different E. coli K-12 strains are just considered E. coli K-12.
Taxons		The taxons corresponding to this interaction.
FuncCatA		Functional category of the first interactor
DescA		Description of the first interactor
FuncCatB		Functional category of the second interactor
DescB		Description of the second interactor

'cons_statistics[date].txt'
			Contains statistics relevant to the produced consensus meta-interactome.

At the moment, this only makes predictions based off presence of the same OGs as in the consensus network.
It needs to verify that both OGs in the predicted PPI are present in the target species.
Redundant predictions (the same interaction from the same taxon ID) are merged.

Uses PSIQUIC service to retrieve IntAct data - see https://github.com/micommunity/psicquic

CHANGES COMPLETE:
Downloads eggNOG map file (LUCA-level and bacteria specific) and IntAct interactions (just bacteria specific)
Generates meta-interactome and rudimentary consensus meta-interactome.
IntAct data cleaned before using (removes "intact" and "chebi" interactors)
A few basic counts (interactors and interactions) are made for meta-interactome and consensus sets
Counts for all consensus interactions are also made across the whole meta-interactome and provided in consensus network
Downloads the eggNOG annotation file for all NOGs but doesn't do anything with it yet
Gets taxon IDs, names, and parent taxon IDs. Adds them to interactions in consensus network but doesn't compare to eliminate redundant taxids
	Checks for parent and child relationships between taxon IDs to limit redundancy.
Gets and maps FuncCat and description annotations (for both LUCA-level and bacteria) to OGs. Use them in the consensus network. 
Provides the option to use a local file in lieu of a downloaded one (especially as PSIQUIC filtering may not provide what we want).
Filters out any input interactions involving non-bacterial taxons (at the meta-interactome building step).
Removes true cross-species interactions at the consensus network building step

IN PROGRESS:
*Are priorities

*Offer the option to append PPI data sets together, as long as they are all in the proper format.
	{{Add the primary Synechocystis and M. loti sets.}}
	{{Generate consensus after this step - will use in subsequent steps if OK.}}
*Filter by FuncCat and produce subsets.
	For subsets, would like to know similar interactions at the protein level.
	E.g., if an OG UF interacts with the same kinds of proteins, what proteins are they AND what else interacts with them in different species?
Get counts and statistics for input data and various interactomes.
*Do interactome prediction for a given proteome. Download proteome on request.
Use protein and species count from eggNOG (it's in the annotation file).
Output interaction sets, filtered by FuncCat (and especially OG UFs).
Perform ANOVA between different FuncCats to see consensus interaction patterns.
	Or at least get interaction frequencies by FuncCat (as in filtering goal above)
*Assign methods to interactions (more general than original data, so we can detect spoke expansion)
	Spoke expansion could actually be filtered out in the input set (use complex:"-") but would rather keep it.
Download a proteome with a search query and set up OG mapping for it.
'''

import glob, gzip, os, re, sys, urllib2, zipfile
from Bio import Entrez
from datetime import date

Entrez.email = 'caufieldjh@vcu.edu'

#Options


#Functions

def get_eggnog_maps(): 
	#Download and unzip the eggNOG Uniprot ID maps
	baseURL = "http://eggnogdb.embl.de/download/eggnog_4.1/id_mappings/uniprot/"
	bactmapfilename = "uniprot-15-May-2015.Bacteria.tsv.gz"	#The Bacteria-specific mapping file
	lucamapfilename = "uniprot-15-May-2015.LUCA.tsv.gz"	#The LUCA mapping file - more generic NOGs
	
	for mapfilename in [bactmapfilename, lucamapfilename]:
		mapfilepath = baseURL + mapfilename
		outfilepath = mapfilename[0:-3]
		if os.path.isfile(mapfilename): 
			print("Found compressed map file on disk: " + mapfilename)
		else:
			response = urllib2.urlopen(mapfilepath)
			print("Downloading from " + mapfilepath)
			compressed_file = open(os.path.basename(mapfilename), "w+b") #Start local compressed file
			chunk = 1048576
			while 1:
				data = (response.read(chunk)) #Read one Mb at a time
				compressed_file.write(data)
				if not data:
					print("\n" + mapfilename + " file download complete.")
					compressed_file.close()
					break
				sys.stdout.write(".")
			
		print("Decompressing map file.")
		with gzip.open(mapfilename) as infile: #Open that compressed file, read and write to uncompressed file
			file_content = infile.read()
			outfile = open(outfilepath, "w+b")
			outfile.write(file_content)
			infile.close()
		outfile.close()
	
def get_interactions():
	#Download and unzip the most recent IntAct version, filtered for bacteria, using REST
	#Just uses IntAct for consistency, but could theoretically include other PSIQUIC compatible DB's
	#May need to add more interactions to the file if not present in IntAct
	#The PSIQUIC interface may also not retrieve all available interactions or may not filter as desired,
	#so script prompts for option.
	#See format description here: https://code.google.com/p/psimi/wiki/PsimiTab27Format
	
	baseURL = "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:%22taxid:2%22?format=tab27"
	intfilename = "protein-interactions.tab"
	
	if os.path.isfile(intfilename): 
		print("Found interaction file on disk: " + intfilename)
	else:
		response = urllib2.urlopen(baseURL)
		print("Downloading from IntAct.")
		intfile = open(os.path.basename(intfilename), "w+b") #Start local file
		chunk = 1048576
		while 1:
			data = (response.read(chunk)) #Read one Mb at a time
			intfile.write(data)
			if not data:
				print("\nInteraction file download complete.")
				intfile.close()
				break
			sys.stdout.write(".")

def get_eggnog_annotations():
	#Downloads and extracts the eggNOG NOG annotations. 
	baseURLs = ["http://eggnogdb.embl.de/download/latest/data/bactNOG/", "http://eggnogdb.embl.de/download/latest/data/NOG/"]
	bactannfilename = "bactNOG.annotations.tsv.gz"	#The annotations for bacteria-specific NOGs
	lucaannfilename = "NOG.annotations.tsv.gz"	#The annotations for other NOGs, but not bacteria-specific NOGs
	
	this_url = 0
	for annfilename in [bactannfilename, lucaannfilename]:
		annfilepath = baseURLs[this_url] + annfilename
		this_url = this_url +1
		outfilepath = annfilename[0:-3]
		if os.path.isfile(annfilename): 
			print("Found compressed annotation file on disk: " + annfilename)
		else:
			response = urllib2.urlopen(annfilepath)
			print("Downloading from " + annfilepath)
			compressed_file = open(os.path.basename(annfilename), "w+b") #Start local compressed file
			chunk = 1048576
			while 1:
				data = (response.read(chunk)) #Read one Mb at a time
				compressed_file.write(data)
				if not data:
					print("\n" + annfilename + " file download complete.")
					compressed_file.close()
					break
				sys.stdout.write(".")
		
		print("Decompressing annotation file.")
		with gzip.open(annfilename) as infile: #Open that compressed file, read and write to uncompressed file
			file_content = infile.read()
			outfile = open(outfilepath, "w+b")
			outfile.write(file_content)
			infile.close()
		outfile.close()
	
def build_meta(mapping_file_list, ppi_data): 
	#Sets up the meta-interactome network.
	#Also creates statistics file about the meta-interactome.
	#This means unique proteins become referred to by their OGs.
	#Interactions are still unique, so two OGs may interact multiple times.
	
	nowstring = (date.today()).isoformat()
	meta_network_filename = "metainteractome" + nowstring + ".txt"
	meta_network_file = open(meta_network_filename, "w")
	
	all_taxids = []
	all_filtered_taxids = []	#Will remove non-bacterial taxids
	map_array = []
	interaction_file = open(ppi_data)
	interaction_array = []
	interaction_filtered_array = []
	protein_array = []
	taxid_species = {} #Dictionary to store taxids and their name and PARENT taxid.
	
	print("Building meta-interactome...")
	print("Arraying map files.")
	for input_map_file in mapping_file_list:
		try:
			map_file = open(input_map_file)
		except IOError as e:
			print("I/O error({0}): {1}".format(e.errno, e.strerror))
		for line in map_file:
			map_array.append((line.rstrip()).split("\t"))
		map_file.close()
		
	print("Arraying interaction file and creating lists of proteins and taxids.")
	
	for line in interaction_file:
		one_interaction = (line.rstrip()).split("\t")
		for taxid in [one_interaction[9], one_interaction[10]]:
			taxid = (((((taxid.split("|"))[0]).lstrip("taxid:")).split("("))[0])
			if taxid not in all_taxids and taxid != "-": #Need to start filtering taxids here so we don't pass bad values to Entrez
				#Why would we get this value anyway? Could be malformed entry as all interactions should have taxids
				all_taxids.append(taxid)	#Just the raw taxid list
		interaction_array.append(one_interaction)	#This is just the raw interaction list at this point
	
	interaction_file.close()
	
	
	print("Finding details for interactor taxids.")
	
	for taxid in all_taxids:
		target_handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
		target_records = Entrez.read(target_handle)
		taxid_name = target_records[0]["ScientificName"]
		taxid_parent = target_records[0]["ParentTaxId"]
		taxid_division = target_records[0]["Division"]
		if taxid_division == "Bacteria":	#Restrict the set to bacteria!
			taxid_species[taxid] = [taxid_name, taxid_parent, taxid_division]
			if taxid not in all_filtered_taxids:
				all_filtered_taxids.append(taxid)
				sys.stdout.write(".")
			#print(taxid_species[taxid])	
	
	print("\nCleaning up data by removing non-protein and non-bacterial interactors.")
	interactions_removed = 0
	for interaction in interaction_array:
		interaction_ok = 1
		for taxid in [interaction[9], interaction[10]]:
			taxid = (((((taxid.split("|"))[0]).lstrip("taxid:")).split("("))[0])
			if taxid not in all_filtered_taxids:
				interaction_ok = 0	#Ensure the interaction won't be kept later
				break				#Ignore this interactor and the other in the pair
		if interaction_ok == 1:
			for protein in interaction[0:2]:			#both proteins in the interacting pair
				if protein[0:6] == "intact" or protein[0:5] == "chebi":	#If the interactor isn't a single protein...
					interaction_ok = 0	#Ensure the interaction won't be kept later
					break				#Ignore this interactor and the other in the pair
				this_protein = protein.lstrip("uniprotkb:")
				if this_protein not in protein_array:
					protein_array.append(this_protein)
			if interaction not in interaction_filtered_array:
				interaction_filtered_array.append(interaction)
		else:
			interactions_removed = interactions_removed +1
				
		
	print("Total taxids: " + str(len(all_filtered_taxids)))
	print("Total raw interactions: " + str(len(interaction_array)))		
	print("Interactions removed: " + str(interactions_removed))	
	
	print("Mapping OGs to " + str(len(protein_array)) + " proteins in " + str(len(interaction_filtered_array)) + " interactions. Proteins mapped:")
	
	protein_OG_maps = {} #Dictionary to save protein-OG mapping specific for this interaction set
	mapped_count = 0
	proteins_without_OG = 0
	for protein in protein_array:
		mapped_count = mapped_count +1
		if mapped_count % 10 == 0:
			sys.stdout.write(".")
			if mapped_count % 100 == 0:
				sys.stdout.write(str(mapped_count))
		matching_OG = protein	#If the protein doesn't map to an OG it retains its original ID
		for mapping in map_array:
			if protein == mapping[0]:
				matching_OG = mapping[1]
				break
		if matching_OG == protein:	#If there wasn't an OG match, count it as a protein without an OG
			proteins_without_OG = proteins_without_OG +1
		protein_OG_maps[protein] = matching_OG
	
	print("\nWriting meta-interactome file. Interactions complete: ")
	interaction_count = 0
	for interaction in interaction_array:		#Write OGs for all interactions.
		interaction_count = interaction_count +1
		if interaction_count % 10 == 0:
			sys.stdout.write(".")
			if interaction_count % 100 == 0:
				sys.stdout.write(str(interaction_count))
				
		for protein in interaction[0:2]:	#Get matching OGs for both proteins in the pair.
			matching_OG = protein_OG_maps[protein.lstrip("uniprotkb:")]
			interaction.append(matching_OG)
			
		interaction_out = "\t".join(interaction) + "\n"
		meta_network_file.write(interaction_out)
		
		#print("mapped " + str(interaction_count) + " - " + matching_OG_A + " vs. " + matching_OG_B)
	
	meta_network_file.close()
	
	meta_stats_filename = "meta_statistics_" + nowstring + ".txt"
	meta_stats_file = open(meta_stats_filename, "w")
	stats_header = ("Unique proteins\tInteractions\tProteins without OG\n")
	meta_stats_file.write(stats_header)
	meta_statistics = []
	for meta_stat in [len(protein_array), interaction_count, proteins_without_OG]:
		meta_statistics.append(str(meta_stat))
	meta_stats_file.write("\t".join(meta_statistics))
	print("\nWrote meta-interactome statistics to " + meta_stats_filename)
	meta_stats_file.close()

	return [meta_network_filename, taxid_species]
	
def build_consensus(metafile, annotation_file_list, taxid_species): 
	#Sets up the consensus meta-interactome network.
	#This is identical to the meta-interactome but compresses interactions into their respective OGs.
	#Interactors without OG assignment are retained and considered single-member OGs.
	
	#This should also count the number of taxons corresponding to an interaction, though it does not ATM
	
	nowstring = (date.today()).isoformat()
	consensus_network_filename = "consensus" + nowstring + ".txt"
	consensus_network_file = open(consensus_network_filename, "w")
	
	consensus_interactors = []	#Consensus interactome interactors - an OG where mapping is possible
	all_interactions = []	#Meta-interactome interactions - unique by protein ID but contain OGs too.
	consensus_interactions = []	#Consensus interactome interactions and associated data. Get written to output file.
	all_annotations = [] #Annotations in file - a bit inefficient to load the whole thing but more searchable this way
	consensus_annotations = {} #Dictionary to store functional category and descriptions of OG interactors.
	
	#First pass: create a list of unique interactions only, using OG IDs
	print("Finding unique interactions.")
	cons_interaction_count = 0
	for line in metafile:
		one_interaction = ((line.rstrip()).split("\t"))
		new_interaction = 0
		if one_interaction[42] not in consensus_interactors:	#Interactor A
			consensus_interactors.append(one_interaction[42])
			new_interaction = 1
		if one_interaction[43] not in consensus_interactors:	#Interactor B	
			consensus_interactors.append(one_interaction[43])
			new_interaction = 1
		if new_interaction == 1:
			cons_interaction_count = cons_interaction_count +1
			consensus_interactions.append([one_interaction[42], one_interaction[43]])
		taxid_A = (((((one_interaction[9].split("|"))[0]).lstrip("taxid:")).split("("))[0])
		taxid_B = (((((one_interaction[10].split("|"))[0]).lstrip("taxid:")).split("("))[0])
		if taxid_A != taxid_B:	#If the two taxids aren't identical, they may still be related or may truly be cross-species.
			#Cross-species PPI get removed.
			taxid_mismatch = 1
			if (taxid_species[taxid_A])[1] == (taxid_species[taxid_B])[1]: #Check if taxids share a parent
				taxid_mismatch = 0
			if (taxid_species[taxid_A])[1] == taxid_B or (taxid_species[taxid_B])[1] == taxid_A: #Check for parent-child relationship
				taxid_mismatch = 0	
		if taxid_mismatch != 1:
			all_interactions.append([one_interaction[42], one_interaction[43], taxid_A, taxid_B])	#Get interactors AND taxid of source of both interactors
			
	#Second pass: count the number of interactions contributing to each consensus
	#Compare taxids across interactions to see how many different sources interaction is seen for (i.e., X different species) 
	#Add counts to each item in consensus_interactions
	#The first count is the total occurence of the given interaction across the full meta-interactome
	print("Counting interaction contributions.")
	all_consensus_taxids = []
	for interaction in consensus_interactions:
		interaction_sources = []	#The list of taxids found to correspond to this interaction.
		original_count = 0
		#This gets a bit complicated.
		for original_interaction in all_interactions:	#For each interaction in the set of all (not OG-compressed consensus) meta-interactome interactions...
			if original_interaction[0] in interaction and original_interaction[1] in interaction: #If both interactors from the meta-interactome interaction are in the consensus interaction...
				original_count = original_count +1	#Add to the count of this interaction across meta-interactome.
				
				for taxid in original_interaction[2:4]:	#For both taxids corresponding to the meta-interactome interaction...
					if taxid not in interaction_sources and taxid != "-":	#If the taxid isn't in the source taxids for this interaction yet...and isn't empty...
						if (taxid_species[taxid])[1] not in interaction_sources:	#Check to see if the sources contain the taxid's parent taxid (if so, it's redundant)
							for source in interaction_sources:
								if (taxid_species[source])[1] == taxid:	#Check to see if taxid is a parent of existing sources (if so, remove children and just use parent)
									interaction_sources.remove(source)
								if (taxid_species[source])[1] == (taxid_species[taxid])[1]: #Check if taxids share a parent (if so, use parent taxid and drop children)
									taxid = (taxid_species[source])[1]	#The problem here is that the parent taxid may not be in taxid_species since it's new to us
									
									#So we look it up and add it!
									target_handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
									target_records = Entrez.read(target_handle)
									taxid_name = target_records[0]["ScientificName"]
									taxid_parent = target_records[0]["ParentTaxId"]
									taxid_species[taxid] = [taxid_name, taxid_parent]
									#This really should be its own function to limit redundancy
									
									interaction_sources.remove(source)
							interaction_sources.append(taxid)
		for source in interaction_sources:	#List all the taxids used across the consensus - does NOT care about parent/child relationships
			if source not in all_consensus_taxids:
				all_consensus_taxids.append(source)
		interaction.append(str(original_count))
		interaction.append(str(len(interaction_sources)))
		interaction.append(" ".join(interaction_sources))
	
	#Third pass: get the functional categories and descriptions of all interactors
	print("Adding interactor annotations.")
	for input_ann_file in annotation_file_list:
		try:
			ann_file = open(input_ann_file)
		except IOError as e:
			print("I/O error({0}): {1}".format(e.errno, e.strerror))
		for line in ann_file:
			all_annotations.append((line.rstrip()).split("\t"))
		ann_file.close()
	
	for interactor in consensus_interactors:
		consensus_annotations[interactor] = ["NA", "NA"] #Should only happen if OG not in description file (e.g. if it's unmapped to an OG)
		for annotation in all_annotations:
			if interactor == annotation[1]:
				consensus_annotations[interactor] = [annotation[4], annotation[5]] #FuncCat and description
				break
	for interaction in consensus_interactions:
		for interactor in interaction[0:2]:
			interaction.append("\t".join(consensus_annotations[interactor]))
				
	print("Writing consensus meta-interactome file.")
	for interaction in consensus_interactions:
		consensus_network_file.write("\t".join(interaction) + "\n")
	consensus_network_file.close()
	
	cons_stats_filename = "cons_statistics_" + nowstring + ".txt"
	cons_stats_file = open(cons_stats_filename, "w")
	stats_header = ("Unique interactors\tInteractions\tTaxids\n")
	cons_stats_file.write(stats_header)
	cons_statistics = []
	for cons_stat in [len(consensus_interactors), cons_interaction_count, len(all_consensus_taxids)]:
		cons_statistics.append(str(cons_stat))
	cons_stats_file.write("\t".join(cons_statistics))
	print("Wrote consensus meta-interactome statistics to " + cons_stats_filename)
	cons_stats_file.close()
	
	return consensus_network_filename
	
'''
def network_analyze(target, target_taxid):	#Not updated yet, just cannibalizing code
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
		target_proteins.append(og_and_prot)
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
	OG_i = []				#All unique OG-OG interactions in the predicted network
	experimental_OG_i = [] #unique OG-OG interactions specific to this and parental taxids
	experimental_OG_i_unique = []
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
		this_OG_i = [ppi[0], ppi[3]]
		if this_OG_i not in OG_i:
			OG_i.append(this_OG_i)			
		if (target_taxid in ppi_taxid_pair or parent_taxid in ppi_taxid_pair):
			if target_taxid == ppi_taxid_pair[0] and target_taxid == ppi_taxid_pair[1]:
				if ppi[4] == "association":
					method = "Experimental results, spoke expansion"
				else:	
					method = "Experimental results"
			else:
				method = "Experimental results, related strain"
			this_reverse_OG_i = [ppi[3], ppi[0]]
			experimental_OG_i.append(this_OG_i)
			if this_OG_i not in experimental_OG_i_unique:
				if this_reverse_OG_i not in experimental_OG_i_unique:
					experimental_OG_i_unique.append(this_OG_i)
		else:
			if taxonomy_compare == 1:	#Determine taxonomic lineage-based similarity if asked
				#Doesn't work quite right yet
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
		if og_and_prot[0] not in predicted_OG_coverage_unique or og_and_prot[0] == "NA":
			out_string = ("\t".join(og_and_prot) + "\n")
			noninteracting_file.write(out_string)
			non_interactor_count = non_interactor_count + 1
	noninteracting_file.close()
	if output_mode == 1:
		stats_output = [target_name, str(len(target_proteins)), str(proteins_are_na), str(non_interactor_count), str(len(target_ogs)), str(match_count), str(len(experimental_OG_i)), str(len(experimental_OG_i_unique)), str(len(predicted_net_unique_alltaxid)), str(len(predicted_OG_coverage_unique))]
	elif output_mode == 2:
		i_exp_proteins = find_proteins_from_OGs(experimental_OG_i_unique, target_proteins)
		i_pred_proteins = find_proteins_from_OGs(OG_i,target_proteins)
		i_exp_OGs = get_OGs_from_network(experimental_OG_i_unique)
		stats_output = [target_name, str(len(target_proteins)), str(len(i_exp_proteins)), str(len(i_pred_proteins)), str(non_interactor_count), str(len(target_ogs)), str(len(i_exp_OGs)), str(len(predicted_OG_coverage_unique)), str(len(target_ogs) - len(predicted_OG_coverage_unique)), str(len(experimental_OG_i_unique)), str(len(OG_i))]
	print("\t".join(stats_output) + "\n")
	
def find_proteins_from_OGs(these_interactions = [], target = []): 
	#this retrieves ALL POSSIBLE protein interactors based on OGs
	unique_OGs = []
	proteins = []
	proteins_out = []
	for interaction in these_interactions:
		for OG in interaction:
			if OG not in unique_OGs:
				unique_OGs.append(OG)
	for OG in unique_OGs:
		for interactor in target:
			if OG in interactor:
				proteins.append(interactor[1])
	for protein in proteins:
		if protein not in proteins_out:
			proteins_out.append(protein)
	return proteins_out
	
def get_OGs_from_network(these_interactions = []):
	unique_OGs = []
	for interaction in these_interactions:
		for OG in interaction:
			if OG not in unique_OGs:
				unique_OGs.append(OG)
	return unique_OGs
'''

def merge_data(list_of_filenames):
	
	nowstring = (date.today()).isoformat()
	merged_file_name = "interactions" + nowstring + ".txt"
	merged_file = open(merged_file_name, "w")
	
	for item in list_of_filenames:
		this_file = open(item)
		line_count = 0
		for line in this_file:
			write_ok = 1
			line_count = line_count +1
			line_contents = ((line.rstrip()).split("\t"))
			for interactor in line_contents[0:2]:
				if interactor == "-":
					print("Empty interactor in " + item + " in line " + str(line_count))
					write_ok = 0
					#Unmapped interactors might be denoted with a -. Don't add them.
			if len(line_contents) != 42:
				print("Format problem in " + item + " line " + str(line_count))
				write_ok = 0
				#Just checking to see if the right number of columns are there
				#Won't write problem lines to the merged file
			if write_ok == 1:
				merged_file.write(line)
		this_file.close()
	return merged_file_name
	
#Main

#Check for eggNOG mapping file and get if needed
mapping_file_list = glob.glob('uniprot*.tsv')
if len(mapping_file_list) >2:
	sys.exit("Only expected two Uniprot to NOG mapping files. Check for duplicates.")
if len(mapping_file_list) <2:
	print("No eggNOG mapping files found or they're incomplete. Retrieving them.")
	get_eggnog_maps()
	mapping_file_list = glob.glob('uniprot*.tsv')
	
#Check for eggNOG annotation file and get if needed
annotation_file_list = glob.glob('*annotations.tsv')
if len(annotation_file_list) >2:
	sys.exit("Only expected two eggNOG annotation files. Check for duplicates.")
if len(annotation_file_list) <2:
	print("No eggNOG annotation files found or they're incomplete. Retrieving them.")
	get_eggnog_annotations()
	annotation_file_list = glob.glob('*annotations.tsv')
	
#Prompt for choice of protein interactions.
#May provide manually or may download, but downloaded set may not be filtered properly.
#Don't need to get interactions if we already have a meta-interactome.
meta_file_list = glob.glob('*metainteractome*.txt')
if len(meta_file_list) >1:
	sys.exit("More than one meta-interactome found. Please use just one at a time.")
if len(meta_file_list) == 0:
	print("No meta-interactome found.")
	ppi_data_option = raw_input("Retreive IntAct bacterial PPI or use local file(s) to build meta-interactome?\n"
	"Enter:\n R for retrieval\n L for local file, or\n M for multiple inputs.\n")
	if ppi_data_option in ["R", "r"]:	#Downloads PPI data from IntAct server. 
		#May not include all PPI available through HTTP IntAct interface.
		ppi_data_filename = "protein-interactions.tab"
		interaction_file_list = glob.glob(ppi_data_filename)
		if len(interaction_file_list) >1:
			sys.exit("One protein interaction file at a time, please! Check for duplicates.")
		if len(interaction_file_list) == 0:
			print("No protein interaction file found. Retrieving it.")
			get_interactions()
			interaction_file_list = glob.glob(ppi_data_filename)
		try:
			interactionfile = open(interaction_file_list[0])
		except IOError as e:
			print("I/O error({0}): {1}".format(e.errno, e.strerror))
	if ppi_data_option in ["L", "l"]:	#Uses a local file, usually a downloaded IntAct PPI set, in PSI-MI Tab27 format
		ppi_data_filename = raw_input("Please provide local filename.\n")
		interaction_file_list = glob.glob(ppi_data_filename)
		if len(interaction_file_list) == 0:
			sys.exit("Can't find a file with that filename.")	
	if ppi_data_option in ["M", "m"]:	#Uses multiple local files in PSI-MI Tab27 format
		adding_files = 1
		interaction_file_list = []
		while adding_files:
			ppi_data_filename = raw_input("Please provide local filename.\n")
			files_present = glob.glob(ppi_data_filename)
			if len(files_present) >0:
				interaction_file_list.append(ppi_data_filename)	#Can be expanded easily later to do batch processing
				print("Added " + ppi_data_filename + " to input list.")
			else:
				print("Can't find a file with that filename. Didn't add.")
			ask_again = raw_input("Add another? Y/N\n")
			if ask_again in ["N", "n"]:
				adding_files = 0
		print("Using the following inputs for the meta-interactome:\n")
		if len(interaction_file_list) == 0:
			sys.exit("Input list is empty. Exiting...")
		for item in interaction_file_list:
			print(item)
		print("Merging into a single file.")
		ppi_data_filename = merge_data(interaction_file_list)
	else:
		sys.exit("Not an option.")

#Load meta-interactome network file
#Needs to be built first.
new_meta = 0
if len(meta_file_list) == 0:
	build_meta_network = raw_input("Build a new meta-interactome? Y/N ")
	if build_meta_network in ["Y", "y"]:
		new_meta_result = build_meta(mapping_file_list, ppi_data_filename)
		new_meta_filename = new_meta_result[0]
		taxids_and_context = new_meta_result[1]
		new_meta = 1
	else:
		sys.exit("Meta-network needed. Exiting.")
try:
	if new_meta == 1:
		metafile = open(new_meta_filename)
	else:
		metafile = open(meta_file_list[0])
except IOError as e:
	print("I/O error({0}): {1}".format(e.errno, e.strerror))	
print("\nUsing " + metafile.name + " as the meta-interactome network.")

#Load consensus network file
#Needs to be built first.
consensus_file_list = glob.glob('*consensus*.txt')
new_consensus = 0
if len(consensus_file_list) >1:
	sys.exit("One consensus network at a time, please!")
if len(consensus_file_list) == 0:
	print("No consensus network file found. Building one.")
	description_file = open("bactNOG.annotations.tsv")
	new_consensus_filename = build_consensus(metafile, annotation_file_list, taxids_and_context)
	new_consensus = 1
try:
	if new_consensus == 1:
		consensusfile = open(new_consensus_filename)
	else:
		consensusfile = open(consensus_file_list[0])
except IOError as e:
	print("I/O error({0}): {1}".format(e.errno, e.strerror))	
print("\nUsing " + consensusfile.name + " as the consensus network.")

#Quit now or ask for next step.
requested = 0
while requested == 0:
	print("\n------------------------------------------------------------")
	request_next = raw_input("Choose from the following options.\n" 
		"A: Generate a subset of the consensus network.\n"
		"B: Generate a predicted interactome for one or more proteomes.\n"
		"X: Exit.\n") 
	if request_next in ["x", "X"]:
		sys.exit("Exiting...")
	if request_next in ["a", "A"]:
		print("Option under construction.")
	if request_next in ["b", "B"]:
		print("Option under construction.")
	print("Choose from the list, please.")

'''	
target_file_list = glob.glob('*target.txt')
if len(target_file_list) == 0:
	print("No target proteome files found, or may not be named properly.")
	one_target_file = raw_input("Provide the filename of one target proteome file")
	else:
		target_file_list.append(one_target_file)
		single_proteome = 1
for filename in target_file_list:
	taxid = (re.split('-', filename))[0]
	try:
		targetfile = open(filename)
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
	network_analyze(targetfile, taxid)
	targetfile.close()
'''
			
sys.exit(0)
