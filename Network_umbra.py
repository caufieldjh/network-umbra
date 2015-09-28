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
			
'meta_statistics[date].txt'
			Contains statistics relevant to the produced meta-interactome.

'taxid_context[date].txt'
			Contains NCBI taxonomy IDs, names, parent IDs, and domains for all input interactions.
			All domains should be Bacteria.
			Used as a reference if meta-interactome and consensus meta-interactome not
			built during the same session.
			
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
			
'subgraph_expansion_[FuncCat]_[date].txt'
			A set of subgraphs of the consensus meta-interactome, filtered by conservation of interactions and function of interactors.
			Each line is one interaction between a consensus interactor and a unique protein, accompanied by the source taxid of the unique protein.
			
'subgraph_expansion_[FuncCat]_nodes_[date].txt'
			Annotation file for the nodes in the expanded subgraphs.


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
Can append input data sets into a single set of interactions. Checks for proper format.
Multiple-OG interactors are handled as single unique OGs in the consensus set and receive annotations.
Verified that self-interactions aren't counted incorrectly.
Subgraph expansion and filtering module is complete.

IN PROGRESS:
*Are priorities

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
from collections import Counter
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
	taxid_context_filename = "taxid_context" + nowstring + ".txt"
	meta_network_file = open(meta_network_filename, "w")
	taxid_context_file = open(taxid_context_filename, "w")
	
	all_taxids = []
	all_filtered_taxids = []	#Will remove non-bacterial taxids
	map_dict = {}	#Uniprot ID to OG dictionary for ALL IDs
	interaction_file = open(ppi_data)
	interaction_array = []
	interaction_filtered_array = []
	protein_array = []
	taxid_species = {} #Dictionary to store taxids and their name and PARENT taxid.
	
	print("Building meta-interactome...")
	print("Setting up protein to OG maps.")
	for input_map_file in mapping_file_list:
		try:
			map_file = open(input_map_file)
		except IOError as e:
			print("I/O error({0}): {1}".format(e.errno, e.strerror))
		for line in map_file:
			one_map = ((line.rstrip()).split("\t"))
			map_dict[one_map[0]] = one_map[1]
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
	
	
	print("Finding details for interactor taxids. This will take some time.")
	
	for taxid in all_taxids:
		unique_taxid_count = 0
		target_handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
		target_records = Entrez.read(target_handle)
		taxid_name = target_records[0]["ScientificName"]
		taxid_parent = target_records[0]["ParentTaxId"]
		taxid_division = target_records[0]["Division"]
		if taxid_division == "Bacteria":	#Restrict the set to bacteria!
			taxid_species[taxid] = [taxid_name, taxid_parent, taxid_division]
			taxid_context_file.write(str(taxid) + "\t" + "\t".join(taxid_species[taxid])+ "\n")
			if taxid not in all_filtered_taxids:
				all_filtered_taxids.append(taxid)
				sys.stdout.write(".")
				unique_taxid_count = unique_taxid_count +1
				if unique_taxid_count % 100 == 0:
					sys.stdout.write(str(unique_taxid_count))
			#print(taxid_species[taxid])
	taxid_context_file.close()	
	
	print("\nCleaning up data by removing non-protein and non-bacterial interactors.")
	interactions_removed = 0
	for interaction in interaction_array:
		interaction_ok = 1
		for taxid in [interaction[9], interaction[10]]:
			taxid = (((((taxid.split("|"))[0]).lstrip("taxid:")).split("("))[0])
			if taxid not in all_filtered_taxids:	#This is where the non-bacterial interactions get removed
				interaction_ok = 0	#Ensure the interaction won't be kept later
				break				#Ignore this interactor and the other in the pair
		if interaction_ok == 1:	#Don't bother to filter proteins if this didn't pass the first filter
			for protein in interaction[0:2]:			#both proteins in the interacting pair
				if protein[0:9] != "uniprotkb":	#Only keep proteins with uniprot IDs
					#Might be nice to keep other protein IDs too but they're rare
					interaction_ok = 0	#Ensure the interaction won't be kept later
					break				#Ignore this interactor and the other in the pair
				this_protein = protein.lstrip("uniprotkb:")
				if this_protein not in protein_array:
					protein_array.append(this_protein)
			if interaction_ok == 1:	# The last filter check for cleaning
				if interaction not in interaction_filtered_array:
					interaction_filtered_array.append(interaction)
			else:
				interactions_removed = interactions_removed +1
		else:
			interactions_removed = interactions_removed +1
				
		
	print("Total taxids: " + str(len(all_filtered_taxids)))
	print("Total raw interactions: " + str(len(interaction_array)))		
	print("Interactions removed: " + str(interactions_removed))	
	
	print("Mapping OGs to " + str(len(protein_array)) + " proteins in " + str(len(interaction_filtered_array)) + " interactions.")
	
	protein_OG_maps = {} #Dictionary to save protein-OG mapping specific for this interaction set
	mapped_count = 0
	proteins_without_OG = 0
	for protein in protein_array:
		mapped_count = mapped_count +1
		
		if protein in map_dict:
			matching_OG = map_dict[protein]
		else:
			matching_OG = protein	#If the protein doesn't map to an OG it retains its original ID
			proteins_without_OG = proteins_without_OG +1
			
		protein_OG_maps[protein] = matching_OG
	
	print("\nWriting meta-interactome file.")
	interaction_count = 0
	for interaction in interaction_filtered_array:		#Write OGs for all (filtered) interactions.
		interaction_count = interaction_count +1
				
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
	
	#First pass: create a list of unique interactors and interactions, using OG IDs
	print("Finding unique interactors and interactions.")
	cons_interaction_count = 0
	for line in metafile:
		one_interaction = ((line.rstrip()).split("\t"))
		new_interaction = 0
		if one_interaction[42] not in consensus_interactors:	#Interactor A's OG or ID if no OG mapped
			consensus_interactors.append(one_interaction[42])
			new_interaction = 1
		if one_interaction[43] not in consensus_interactors:	#Interactor B's OG or ID if no OG mapped	
			consensus_interactors.append(one_interaction[43])
			new_interaction = 1
		if new_interaction == 1:
			cons_interaction_count = cons_interaction_count +1
			consensus_interactions.append([one_interaction[42], one_interaction[43]])
		taxid_A = (((((one_interaction[9].split("|"))[0]).lstrip("taxid:")).split("("))[0])
		taxid_B = (((((one_interaction[10].split("|"))[0]).lstrip("taxid:")).split("("))[0])
		taxid_mismatch = 0	#Assume that the two taxids are the same by default
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
	#The second count is the number of different, unique taxids (species or at least distant strains) corresponding to the interaction
	
	print("Counting interaction contributions. This may take a while.")
	print("Consensus interactions checked, out of " + str(len(consensus_interactions)) +":")
	
	all_consensus_taxids = []
	con_interactions_checked = 0
	for interaction in consensus_interactions:
		con_interactions_checked = con_interactions_checked +1
		if con_interactions_checked % 10 == 0:
			sys.stdout.write(".")
		if con_interactions_checked % 100 == 0:
			sys.stdout.write(str(con_interactions_checked))
		interaction_sources = []	#The list of taxids found to correspond to this interaction.
		original_count = 0
		#This gets a bit complicated.
		for original_interaction in all_interactions:	#For each interaction in the set of all (not OG-compressed consensus) meta-interactome interactions...
			original_interaction_slim = original_interaction[0:2]
			original_interaction_slim_rev = [original_interaction[1], original_interaction[0]]
			if interaction == original_interaction_slim or interaction == original_interaction_slim_rev:	#If the original interaction OR its reverse matches the consensus interaction...
				original_count = original_count +1	#Add to the count of this interaction across the meta-interactome.
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
	print("\nAdding interactor annotations.")
	for input_ann_file in annotation_file_list:
		try:
			ann_file = open(input_ann_file)
		except IOError as e:
			print("I/O error({0}): {1}".format(e.errno, e.strerror))
		for line in ann_file:
			all_annotations.append((line.rstrip()).split("\t"))
		ann_file.close()
	
	multiple_og_count = 0	#The count of interactors mapping to >1 OG. Are treated as single OGs as this may be biologically meaningful
	for interactor in consensus_interactors:
		consensus_annotations[interactor] = ["NA", "NA"] #Should only happen if OG not in description file (e.g. if it's unmapped to an OG)
		if "," in interactor:	#Meaning it maps to multiple OGs, so we need to annotate all of them
			#print(interactor)
			this_mult_og_count = 0	#Keeps track of multiple OG sets. Usually just two or three different OGs at most.
			consensus_annotations[interactor] = ["", ""]
			multiple_og_count = multiple_og_count +1
			multiple_ogs = interactor.split(",")
			for og in multiple_ogs:
				this_mult_og_count = this_mult_og_count +1
				for annotation in all_annotations:
					if og == annotation[1]:
						#Concatenate each FuncCat, separated by |
						(consensus_annotations[interactor])[0] = (consensus_annotations[interactor])[0] + annotation[4]
						if this_mult_og_count != len(multiple_ogs):
							(consensus_annotations[interactor])[0] = (consensus_annotations[interactor])[0] + "|"
						#Concatenate each description, separated by |
						(consensus_annotations[interactor])[1] = (consensus_annotations[interactor])[1] + annotation[5]
						if this_mult_og_count != len(multiple_ogs):
							(consensus_annotations[interactor])[1] = (consensus_annotations[interactor])[1] + "|"
						break
		else:
			for annotation in all_annotations:
				if interactor == annotation[1]:
					consensus_annotations[interactor] = [annotation[4], annotation[5]] #FuncCat and description
					break
	for interaction in consensus_interactions:
		for interactor in interaction[0:2]:
			interaction.append("\t".join(consensus_annotations[interactor]))
	
	print("Consensus meta-interactome involves " + str(len(consensus_interactors)) +
			" interactors and " + str(cons_interaction_count) + " interactions.")
	print("It involves " + str(len(all_consensus_taxids)) + " unique taxids, " +
			"though some may be closely related.")
	print(str(multiple_og_count) + " interactors map to more than one OG.")
	
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
	
def subset_expansion(metafile, consensusfile):
	print("\nSubset expansion will filter consensus interactions by functional category" +
			 " and by conservation across taxonomic groups.\n" +
			 "It will produce a set of subgraphs, where each graph involves a consensus" +
			 " interactor and ALL of its interactions in the meta-interactome.\n" +
			 "These graphs will contain taxonomy annotations for each interaction" +
			 " and can be split in network analysis software, e.g. Cytoscape.\n")
	print("Functional categories:\n" 
			"INFORMATION STORAGE AND PROCESSING\n"
			"[J] Translation, ribosomal structure and biogenesis\n" 
			"[A] RNA processing and modification\n" 
			"[K] Transcription\n"
			"[L] Replication, recombination and repair\n" 
			"[B] Chromatin structure and dynamics\n"
			"CELLULAR PROCESSES AND SIGNALING\n"
			"[D] Cell cycle control, cell division, chromosome partitioning\n" 
			"[Y] Nuclear structure\n"
			"[V] Defense mechanisms\n"
			"[T] Signal transduction mechanisms\n"
			"[M] Cell wall/membrane/envelope biogenesis\n" 
			"[N] Cell motility\n" 
			"[Z] Cytoskeleton\n"
			"[W] Extracellular structures\n"
			"[U] Intracellular trafficking, secretion, and vesicular transport\n" 
			"[O] Posttranslational modification, protein turnover, chaperones\n" 
			"METABOLISM\n"
			"[C] Energy production and conversion\n"
			"[G] Carbohydrate transport and metabolism\n" 
			"[E] Amino acid transport and metabolism\n"
			"[F] Nucleotide transport and metabolism\n" 
			"[H] Coenzyme transport and metabolism\n" 
			"[I] Lipid transport and metabolism\n"
			"[P] Inorganic ion transport and metabolism\n"
			"[Q] Secondary metabolites biosynthesis, transport and catabolism\n" 
			"POORLY CHARACTERIZED\n"
			"[R] General function prediction only\n"
			"[S] Function unknown\n")
			
	func_filter = raw_input("Filter interactors for which functional category? (Type X for interactors of unknown function.)\n")
	search_unknowns = 0
	if func_filter in ["x", "X"]:
		search_unknowns = 1
		print("Filtering for interactors of unknown function. Interactors marked NA will not be included.")
	
	consensus_interactions = []	#Contains whole line (one interaction) from consensus
	consensus_interactions_filtered = [] #Interactions filtered by FuncCat
	consensus_interactions_taxfilt = []	#Interactions filtered by FuncCat and taxids
	consensus_interactors_filtered = {}	#Contains interactor, FuncCat, and description (filtered by FuncCat)
	consensus_interactors_taxfilt = {}	#Contains interactor, FuncCat, and description (filtered by FuncCat and 
										#by participation in an interaction passing the taxid filter)
	expanded_interactions = {}			#Keys are consensus interactors. Values are all unique proteins (and sources) they interact with.
										#Actually a dict of lists of lists. Fun.
	max_taxon_range = 1	#The greatest count of different taxids per interaction, across the whole consensus
	all_interactions = []	#All interactions in the meta-interactome
	protein_annotations = {}	#Annotations (from IntAct) for unique proteins. No FuncCats here.
	
	print("Filtering consensus interactors by function.")
	for line in consensusfile:	#Filter interactors by FuncCat
		one_consensus_interaction = ((line.rstrip()).split("\t"))
		consensus_interactions.append(one_consensus_interaction)
		
		if one_consensus_interaction[5] != "NA":	#For interactor A
			if search_unknowns == 1:
				if "R" in one_consensus_interaction[5] or "S" in one_consensus_interaction[5]:
					consensus_interactors_filtered[one_consensus_interaction[0]] = [one_consensus_interaction[5], one_consensus_interaction[6]]
			else:
				if func_filter in one_consensus_interaction[5]:
					consensus_interactors_filtered[one_consensus_interaction[0]] = [one_consensus_interaction[5], one_consensus_interaction[6]]
		if one_consensus_interaction[7] != "NA":	#For interactor B
			if search_unknowns == 1:
				if "R" in one_consensus_interaction[7] or "S" in one_consensus_interaction[7]:
					consensus_interactors_filtered[one_consensus_interaction[1]] = [one_consensus_interaction[7], one_consensus_interaction[8]]
			else:
				if func_filter in one_consensus_interaction[7]:
					consensus_interactors_filtered[one_consensus_interaction[1]] = [one_consensus_interaction[7], one_consensus_interaction[8]]
					
	consensusfile.close()	#May not want to close file if we plan on doing multiple filters during same session.
	
	for interaction in consensus_interactions:
		if interaction[0] in consensus_interactors_filtered or interaction[1] in consensus_interactors_filtered:
			consensus_interactions_filtered.append(interaction)
			if interaction[3] > max_taxon_range:
				max_taxon_range = interaction[3]
		
	print("The maximum for this filter will be " + str(max_taxon_range) + " different taxids.")
	tax_filter = raw_input("Select for at least how many different taxonomic groups?\n")
	
	for interaction in consensus_interactions_filtered:
		if interaction[3] >= tax_filter:
			consensus_interactions_taxfilt.append(interaction)
			
			if interaction[5] != "NA":	#For interactor A
				if search_unknowns == 1:
					if "R" in interaction[5] or "S" in interaction[5]:
						consensus_interactors_taxfilt[interaction[0]] = [interaction[5], interaction[6]]
				else:
					if func_filter in interaction[5]:
						consensus_interactors_taxfilt[interaction[0]] = [interaction[5], interaction[6]]
			if interaction[7] != "NA":	#For interactor B
				if search_unknowns == 1:
					if "R" in interaction[7] or "S" in interaction[7]:
						consensus_interactors_taxfilt[interaction[1]] = [interaction[7], interaction[8]]
				else:
					if func_filter in interaction[7]:
						consensus_interactors_taxfilt[interaction[1]] = [interaction[7], interaction[8]]
	
	print("Generated list of filtered consensus interactors. Searching meta-interactome.")
	
	for line in metafile:	#Set up the meta-interactome file first
		one_interaction = (line.rstrip()).split("\t")
		all_interactions.append(one_interaction)	#This is just the raw interaction list at this point
	
	metafile.close()
		
	for interactor in consensus_interactors_taxfilt:
		expanded_interactions[interactor] = []
		for interaction in all_interactions:	#Search meta-interactome for matching interactions; return unique all proteins and corresponding organisms
			if interaction[42] == interactor:
				taxid = (((((interaction[10].split("|"))[0]).lstrip("taxid:")).split("("))[0])
				protein = interaction[1].lstrip("uniprotkb:")
				if protein not in protein_annotations:
					protein_annotations[protein] = [interaction[23], interaction[43]]
				protein_and_source = [protein, taxid]
				(expanded_interactions[interactor]).append(protein_and_source)
			if interaction[43] == interactor:
				if interaction[0] != interaction[1]:	#Avoid adding self-interactions twice.
					taxid = (((((interaction[9].split("|"))[0]).lstrip("taxid:")).split("("))[0])
					protein = interaction[0].lstrip("uniprotkb:")
					if protein not in protein_annotations:
						protein_annotations[protein] = [interaction[22], interaction[42]]
					protein_and_source = [protein, taxid]
					(expanded_interactions[interactor]).append(protein_and_source)
	
	nowstring = (date.today()).isoformat()
	subgraph_file_name = "subgraph_expansion_" + func_filter + "_" + nowstring + ".txt"
	subgraph_node_file_name = "subgraph_expansion_" + func_filter + "_nodes_" + nowstring + ".txt"
	subgraph_file = open(subgraph_file_name, "w")
	subgraph_node_file = open(subgraph_node_file_name, "w")
	
	#print(consensus_interactors_taxfilt)	
	#print(expanded_interactions)
	
	print("Writing subgraph expansion file and node annotation file.")
	for consensus_interactor in expanded_interactions:
		for interaction in expanded_interactions[consensus_interactor]:
			#print(consensus_interactor + "\t" + "\t".join(interaction))
			subgraph_file.write(consensus_interactor + "\t" + "\t".join(interaction) + "\n")
	
	#Protein annotations are kind of a mess but that's because the interaction data table combines interactor annotations into single columns.
	#It's also difficult to know what kind of annotation to expect. 
	#All are included here, for now.
	
	for interactor in consensus_interactors_taxfilt:
		subgraph_node_file.write(interactor + "\t" + "\t".join(consensus_interactors_taxfilt[interactor]) + "\t-\n")
	for protein in protein_annotations:
		subgraph_node_file.write(protein + "\t-\t" + "\t".join(protein_annotations[protein]) + "\n")
		
	print("Done.")
	
	
def predict_interactome():
	print("Interactome prediction under construction.")
	
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
ppi_data_filename = ""
if len(meta_file_list) >1:
	sys.exit("More than one meta-interactome found. Please use just one at a time.")
if len(meta_file_list) == 0:
	print("No meta-interactome found.")
	while ppi_data_filename == "":
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
	if new_meta == 1:	#If we just build a meta-interactome we have taxid details already
		new_consensus_filename = build_consensus(metafile, annotation_file_list, taxids_and_context)
	else:	#Otherwise we need to read taxid details from file - just rebuild dict from it
		taxid_ref_list = glob.glob('taxid_context*.txt')
		taxids_and_context = {}
		if len(taxid_ref_list) >1:
			sys.exit("Something went wrong - more than one taxid context file found.")
		if len(taxid_ref_list) == 0:
			sys.exit("Something went wrong - no taxid context file found.")
		taxid_ref_file = open(taxid_ref_list[0])
		for line in taxid_ref_file:
			content = ((line.rstrip()).split("\t"))
			taxids_and_context[content[0]] = [content[1], content[2], content[3]]
		taxid_ref_file.close()
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
		"A: Generate expanded subgraphs of the consensus network, filtering by function.\n"
		"B: Generate a predicted interactome for one or more proteomes.\n"
		"X: Exit.\n") 
	if request_next in ["x", "X"]:
		sys.exit("Exiting...")
	if request_next in ["a", "A"]:
		subset_expansion(metafile, consensusfile)
	if request_next in ["b", "B"]:
		predict_interactome()
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
