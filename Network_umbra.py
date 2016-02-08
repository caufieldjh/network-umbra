#!/usr/bin/python
#Network_umbra.py
'''
Predicts interactions in a protein interaction network based off a meta-interactome network.
Uses eggNOG v.4.1.
Written for Python 2.7. Not tested with Python 3.

REQUIRES: Biopython 1.65 or more recent
			Also needs at least 5 GB of available disk space to accomodate data files and output
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
			
'interactome_statistics_[date].txt'
			Counts of interactors - proteins and OGs - participating in predicted interactomes.
			Contains the following counts per input proteome:
			Name, taxid, Proteins, ProteinsNotInPPI, ProteinsWithExpPPI, ProteinsWithPredPPI, 
			UniqueOGs, OGsWithoutInteractions, OGsWithExpInt, OGsWithPredInt, ExpOGIntNet, OGIntInPredNet

Uses PSIQUIC service to retrieve IntAct data - see https://github.com/micommunity/psicquic

CHANGES COMPLETE:
Downloads eggNOG map file (LUCA-level and bacteria specific) and IntAct interactions (just bacteria specific)
Generates meta-interactome and rudimentary consensus meta-interactome.
IntAct data cleaned before using (removes "intact" and "chebi" interactors)
A few basic counts (interactors and interactions) are made for meta-interactome and consensus sets
Counts for all consensus interactions are also made across the whole meta-interactome and provided in consensus network
Downloads the eggNOG annotation file for all NOGs
Gets taxon IDs, names, and parent taxon IDs. Adds them to interactions in consensus network but doesn't compare to eliminate redundant taxids
	Checks for parent and child relationships between taxon IDs to limit redundancy.
Gets and maps FuncCat and description annotations (for both LUCA-level and bacteria) to OGs. Uses them in the consensus network. 
Provides the option to use a local file in lieu of a downloaded one (especially as PSIQUIC filtering may not provide what we want).
Filters out any input interactions involving non-bacterial taxons (at the meta-interactome building step).
Removes true cross-species interactions at the consensus network building step
Can append input data sets into a single set of interactions. Checks for proper format.
Multiple-OG interactors are handled as single unique OGs in the consensus set and receive annotations.
Verifies that self-interactions aren't counted incorrectly.
Subgraph expansion and filtering module is complete.
Can provide contribution counts for each taxid in the consensus network.
Basic interactome prediction (based on consensus) is complete. Will download proteomes from Uniprot on request.
Added unique proteins to predicted interactome output.
Output statistics about predicted networks after building them.
The main ID conversion file is now used in lieu of the Uniprot-specific flat file.
Viral proteins can now be used as long as the useViruses option is on (though host vs. virus interactions are still filtered out)
	Phage proteins won't work very well right now as the eggNOG viral NOGs don't include very many phage proteins 
	(e.g. just 37 out of 66 lambda proteins and that's one of the phages with better coverage)
	But there are *some* phage proteins in NOGs - they may just not be in the eggNOG protein_id_conversion file

IN PROGRESS:
* <- Are priorities

Check on inputs for useViruses option - right now, most viral proteins don't get mapped to NOGs
	It is likely that viral proteins are not included in the eggNOG protein ID conversion file.
	
Evaluate predictions on the basis of OG members (large OGs have less predictive power, at least without sequence alignments)
Use protein and species count from eggNOG (it's in the annotation file).
Perform ANOVA between different FuncCats to see consensus interaction patterns.
	Or at least get interaction frequencies by FuncCat (as in filtering goal above)
*Assign methods to interactions (more general than original data, so we can detect spoke expansion)
	Spoke expansion could actually be filtered out in the input set (use complex:"-") but would rather keep it.

'''

import glob, gzip, operator, os, re, requests, sys, urllib2, zipfile
from Bio import Entrez
from bs4 import BeautifulSoup
from collections import Counter
from datetime import date

Entrez.email = 'caufieldjh@vcu.edu'

#Options
useViruses = False	#Option for using eggNOG's viral OGs. Requires the filters permitting only Bacteria to be modified
					#Also requires the viral OGs to be downloaded and added.
					#This option needs to be set True BEFORE the Uniprot to OG map is built or it won't include proteins from viruses
					
useNonRefProteomes = True	#Option to search non-reference Uniprot proteomes in the interactome prediction module
#Retrieving non-reference proteomes sometimes returns an empty response.
#This happens with proteomes only in UniParc (e.g., if they are redundant)
#In those cases, we reject the search result.

#Functions

def get_eggnog_maps(): 
	#Download and unzip the eggNOG ID conversion file 
	#Filters file to just Uniprot IDs; the resulting file is the map file.
	#One Uniprot ID may correspond to multiple OGs - e.g. COG1234,COG3810,COG9313. 
	#these cases are considered OGs in their own right as this may indicate a pattern of conserved sequences on its own 
	baseURL = "http://eggnogdb.embl.de/download/eggnog_4.1/"
	convfilename = "eggnog4.protein_id_conversion.tsv.gz"	#File contains ALL database identifiers and corresponding proteins
	
	convfilepath = baseURL + convfilename
	outfilepath = convfilename[0:-3]
	dl_convfile = 1	#If 1, we need to download
	if os.path.isfile(convfilename): #Already have the compressed file, don't download
		print("Found compressed ID conversion file on disk: %s" % convfilename)
		decompress_convfile = 1
		dl_convfile = 0
	if os.path.isfile(outfilepath): #Already have the decompressed file don't download
		print("Found ID conversion file on disk: %s" % outfilepath)
		decompress_convfile = 0
		dl_convfile = 0
	
	if dl_convfile == 1:
		print("Downloading ID mapping file - this file is ~400 Mb compressed so this may take some time.")
		print("Downloading from %s" % convfilepath)
		response = urllib2.urlopen(convfilepath)
		compressed_file = open(os.path.basename(convfilename), "w+b") #Start local compressed file
		chunk = 1048576
		while 1:
			data = (response.read(chunk)) #Read one Mb at a time
			compressed_file.write(data)
			if not data:
				print("\n%s file download complete." % convfilename)
				compressed_file.close()
				break
			sys.stdout.flush()
			sys.stdout.write(".")
		decompress_convfile = 1
		
	if decompress_convfile == 1:
		print("Decompressing map file. Lines written, in millions:")
		#Done in chunks since it's a large file
		with gzip.open(convfilename) as infile: #Open that compressed file, read and write to uncompressed file
			outfile = open(outfilepath, "w+b")
			linecount = 0
			for line in infile:
				outfile.write(line)
				linecount = linecount +1
				if linecount % 100000 == 0:
						sys.stdout.write(".")
				if linecount % 1000000 == 0:
						sys.stdout.flush()
						sys.stdout.write(str(linecount/1000000))
			infile.close()
		newconvfilename = outfilepath
		outfile.close()
	
	#Download and decompress member NOG files (2 of them)
	nogURL = baseURL + "data/NOG/"
	nogfilename = "NOG.members.tsv.gz"
	bactnogURL = baseURL + "data/bactNOG/"
	bactnogfilename = "bactNOG.members.tsv.gz" 
	all_nog_locations = [[nogURL, nogfilename], [bactnogURL, bactnogfilename]]
	
	if useViruses == True:
		virnogURL = baseURL + "data/viruses/Viruses/"
		virnogfilename = "Viruses.members.tsv.gz"
		all_nog_locations.append([virnogURL, virnogfilename])
	
	for location in all_nog_locations:
		baseURL = location[0]
		memberfilename = location[1]
		memberfilepath = baseURL + memberfilename
		outfilepath = memberfilename[0:-3]
		if os.path.isfile(memberfilename): 
			print("\nFound compressed NOG membership file on disk: %s" % memberfilename)
			decompress_memberfile = 1
		if os.path.isfile(outfilepath): 
			print("\nFound NOG membership file on disk: %s" % outfilepath)
			decompress_memberfile = 0
		else:
			print("\nDownloading NOG membership file - this may take some time.")
			print("Downloading from %s" % memberfilepath)
			response = urllib2.urlopen(memberfilepath)
			compressed_file = open(os.path.basename(memberfilename), "w+b") #Start local compressed file
			chunk = 1048576
			while 1:
				data = (response.read(chunk)) #Read one Mb at a time
				compressed_file.write(data)
				if not data:
					print("\n%s file download complete." % memberfilename)
					compressed_file.close()
					break
				sys.stdout.flush()
				sys.stdout.write(".")
			decompress_memberfile = 1
			
		if decompress_memberfile == 1:
			print("Decompressing NOG membership file %s" % memberfilename)
			#Done in chunks since it's a large file
			with gzip.open(memberfilename) as infile: #Open that compressed file, read and write to uncompressed file
				outfile = open(outfilepath, "w+b")
				linecount = 0
				for line in infile:
					outfile.write(line)
					linecount = linecount +1
					if linecount % 100000 == 0:
						sys.stdout.write(".")
					if linecount % 1000000 == 0:
						sys.stdout.flush()
						sys.stdout.write(str(linecount/1000000))
				infile.close()
			outfile.close()
			
	#Clean up by removing compressed files
	print("\nRemoving compressed files.")
	all_compressed_files = [convfilename, nogfilename, bactnogfilename]
	if useViruses == True:
		all_compressed_files.append(virnogfilename)
	for filename in all_compressed_files:
		if os.path.isfile(filename):
			os.remove(filename)
	
	#Load and filter the ID conversion file as dictionary
	print("Parsing ID conversion file. Lines read, in millions:")
	with open(convfilename[0:-3]) as infile:
		id_dict = {}	#Dictionary of eggNOG protein IDs with database IDs as keys
		#Gets filtered down to relevant database IDs (i.e., Uniprot IDs)
		linecount = 0
		for line in infile:
			linecount = linecount +1
			line_raw = ((line.rstrip()).split("\t"))	#Protein IDs are split for some reason; merge them
			one_id_set = [line_raw[0] + "." + line_raw[1], line_raw[2], line_raw[3]]
			if "UniProt_AC" in one_id_set[2]:
				id_dict[one_id_set[1]] = one_id_set[0]
			if linecount % 100000 == 0:
				sys.stdout.write(".")
			if linecount % 1000000 == 0:
				sys.stdout.flush()
				sys.stdout.write(str(linecount/1000000))
		infile.close()

	#Use filtered ID conversion input to map to NOG members
	print("\nReading NOG membership files.")
	all_nog_filenames = [nogfilename[0:-3], bactnogfilename[0:-3]]
	nog_members = {}	#Dictionary of NOG ids with protein IDs as keys (need to split entries for each)
	nog_count = 0
	for filename in all_nog_filenames:
		temp_nog_members = {}	#We will have duplicates within each set but don't want to lose the information.
		print("Reading from %s" % filename)
		with open(filename) as infile:
			for line in infile:
				nog_count = nog_count +1
				line_raw = ((line.rstrip()).split("\t"))
				nog_id = line_raw[1]
				line_members = line_raw[5].split(",")
				for protein_id in line_members:			#The same protein could be in more than one OG at the same level
					if protein_id in temp_nog_members:
						temp_nog_members[protein_id] = temp_nog_members[protein_id] + "," + nog_id
					else:
						temp_nog_members[protein_id] = nog_id
			infile.close()
		nog_members.update(temp_nog_members)
	
	upids_length = str(len(id_dict))
	nogs_length = str(nog_count)
	proteins_length = str(len(nog_members))
	
	print("Mapping %s Uniprot IDs to %s NOGs through %s eggNOG protein IDs:" % (upids_length, nogs_length, proteins_length))
	upid_to_NOG = {}	#Conversion dictionary. Values are OGs, keys are UPIDs.
	mapped_count = 0	#upids mapped to nogs.
	for upid in id_dict:
		if id_dict[upid] in nog_members:
			upid_to_NOG[upid] = nog_members[id_dict[upid]]
			mapped_count = mapped_count +1
			if mapped_count % 100000 == 0:
				sys.stdout.write(".")
			if mapped_count % 1000000 == 0:
				sys.stdout.flush()
				sys.stdout.write(str(mapped_count/1000000))
		
	#Use this mapping to build map file, named "uniprot_og_maps_*.txt"
	print("Writing map file.")
	nowstring = (date.today()).isoformat()
	mapfilename = "uniprot_og_maps_" + nowstring + ".txt"
	mapfile = open(mapfilename, "w+b")
	for mapping in upid_to_NOG:
		mapfile.write(mapping + "\t" + upid_to_NOG[mapping] + "\n")	#Each line is a uniprot ID and an OG id
	mapfile.close() 
	
	
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
		print("Found interaction file on disk: %s" % intfilename)
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
	annfilenames = [bactannfilename, lucaannfilename]
	
	if useViruses == True:
		baseURLs.append("http://eggnogdb.embl.de/download/latest/data/viruses/Viruses/")
		annfilenames.append("Viruses.annotations.tsv.gz")
	
	this_url = 0
	for annfilename in annfilenames:
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
		
	print("\nRemoving compressed files.")
	all_compressed_files = [bactannfilename, lucaannfilename]
	for filename in all_compressed_files:
		os.remove(filename)
	
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
	all_filtered_taxids = []	#Will remove non-bacterial taxids, unless useViruses is on
	map_dict = {}	#Uniprot ID to OG dictionary for ALL IDs
	interaction_file = open(ppi_data)
	interaction_array = []
	interaction_filtered_array = []
	protein_array = []
	taxid_species = {} #Dictionary to store taxids and their name and PARENT taxid.
	
	print("Building meta-interactome...")
	print("Setting up protein to OG maps.")
	#We preferentially use bacteria mapping first
	for input_map_file in mapping_file_list:
		try:
			map_file = open(input_map_file)
		except IOError as e:
			print("I/O error({0}): {1}".format(e.errno, e.strerror))
		for line in map_file:
			one_map = ((line.rstrip()).split("\t"))
			map_dict[one_map[0]] = one_map[1]
		map_file.close()
	
	#for item in map_dict:
	#	print(item + "\t" + map_dict[item])
	
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
		
		if taxid in ["Taxid interactor A", "Taxid interactor B"]:	#This means the header wasn't removed.
			continue
			
		unique_taxid_count = 0
		#print(str(taxid))
		target_handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
		target_records = Entrez.read(target_handle)
		taxid_name = target_records[0]["ScientificName"]
		taxid_parent = target_records[0]["ParentTaxId"]
		taxid_division = target_records[0]["Division"]
		#print(taxid_division)
		taxid_filter = ["Bacteria"]
		if useViruses == True:
			virus_types = ["Phages", "Viruses"]
			for virus_type in virus_types:
				taxid_filter.append(virus_type)	
		if taxid_division in taxid_filter:	#Restrict the set to bacteria, unless useViruses is on
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
	
	if useViruses == False:
		print("\nCleaning up data by removing non-protein and non-bacterial interactors.")
	else:
		print("\nCleaning up data by removing non-protein and non-bacterial and non-viral interactors.")
	interactions_removed = 0
	for interaction in interaction_array:
		interaction_ok = 1
		for taxid in [interaction[9], interaction[10]]:
			taxid = (((((taxid.split("|"))[0]).lstrip("taxid:")).split("("))[0])
			if taxid not in all_filtered_taxids:	#This is where the non-bacterial (and non-viral, if useViruses is on) interactions get removed
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
				
		
	print("Total taxids: %s" % (len(all_filtered_taxids)))
	print("Total raw interactions: %s" % (len(interaction_array)))		
	print("Interactions removed: %s" % (interactions_removed))	
	
	print("Mapping OGs to %s proteins in %s interactions." % (len(protein_array), len(interaction_filtered_array)))
	
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
	print("\nWrote meta-interactome statistics to %s" % meta_stats_filename)
	meta_stats_file.close()

	return [meta_network_filename, taxid_species]
	
def build_consensus(metafile, annotation_file_list, taxid_species): 
	#Sets up the consensus meta-interactome network.
	#This is identical to the meta-interactome but compresses interactions into their respective OGs.
	#Interactors without OG assignment are retained and considered single-member OGs.
	
	nowstring = (date.today()).isoformat()
	consensus_network_filename = "consensus" + nowstring + ".txt"
	consensus_network_file = open(consensus_network_filename, "w")
	
	consensus_interactors = []	#Consensus interactome interactors - an OG where mapping is possible
	all_interactions_taxids = []	#All interactions in the meta-interactome, but just with OGs and taxids 
	all_interactions_simple = [] #All interactions in the meta-interactome, but just with OGs
	consensus_interactions = []	#Consensus interactome interactions first, then associated data. Get written to output file.
	all_annotations = [] #Annotations in file - a bit inefficient to load the whole thing but more searchable this way
	consensus_annotations = {} #Dictionary to store functional category and descriptions of OG interactors.
	
	#Load the meta-interactome file, removing true cross-species interactions
	for line in metafile:
		one_interaction = ((line.rstrip()).split("\t"))
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
			all_interactions_taxids.append([one_interaction[42], one_interaction[43], taxid_A, taxid_B])
			all_interactions_simple.append([one_interaction[42], one_interaction[43]])
	
	#First pass: create a list of unique interactors and interactions, using OG IDs
	print("Finding unique interactors and interactions. Interactions found:")
	cons_interaction_count = 0
	
	for interaction in all_interactions_simple:
		interaction_rev = [interaction[1], interaction[0]]
		new_interaction = 0
		for interactor in interaction:	#Interactor A or B's OG or ID if no OG mapped
			if interactor not in consensus_interactors:
				consensus_interactors.append(interactor)
		if interaction not in consensus_interactions and interaction_rev not in consensus_interactions:
			new_interaction = 1
		if new_interaction == 1:
			cons_interaction_count = cons_interaction_count +1
			if cons_interaction_count % 100 == 0:
					sys.stdout.write(".")
			if cons_interaction_count % 1000 == 0:
					sys.stdout.write(str(cons_interaction_count))
			consensus_interactions.append(interaction)
			
	#Second pass: count the number of interactions contributing to each consensus
	#Compare taxids across interactions to see how many different sources interaction is seen for (i.e., X different species) 
	#Add counts to each item in consensus_interactions
	#The first count is the total occurence of the given interaction across the full meta-interactome
	#The second count is the number of different, unique taxids (species or at least distant strains) corresponding to the interaction
	
	print("\nCounting interaction contributions. This may take a while.")
	print("Consensus interactions checked, out of %s:" % (len(consensus_interactions)))
	
	all_consensus_taxids = []
	con_interactions_counted = 0
	for interaction in consensus_interactions:
		con_interactions_counted = con_interactions_counted +1
		if con_interactions_counted % 100 == 0:
			sys.stdout.write(".")
		if con_interactions_counted % 1000 == 0:
			sys.stdout.write(str(con_interactions_counted))
		interaction_sources = []	#The list of taxids found to correspond to this interaction.
		original_count = 0
		#This gets a bit complicated.
		for original_interaction in all_interactions_taxids:	#For each interaction in the set of all (not OG-compressed consensus) meta-interactome interactions...
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
		#print(interaction)
	
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
	annotation_count = 0
	for interactor in consensus_interactors:
		annotation_count = annotation_count +1
		if annotation_count % 100 == 0:
			sys.stdout.write(".")
		if annotation_count % 1000 == 0:
			sys.stdout.write(str(annotation_count))
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
	
	print("\nConsensus meta-interactome involves %s" +
			" interactors and %s interactions." % (len(consensus_interactors), cons_interaction_count))
	print("It involves %s unique taxids, " +
			"though some may be closely related." % (len(all_consensus_taxids)))
	print("%s interactors map to more than one OG." % multiple_og_count)
	
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
					print("Empty interactor in %s in line %s" + str() % (item, line_count))
					write_ok = 0
					#Unmapped interactors might be denoted with a -. Don't add them.
			if len(line_contents) != 42:
				print("Format problem in %s line %s" % (item, line_count))
				write_ok = 0
				#Just checking to see if the right number of columns are there
				#Won't write problem lines to the merged file
			if write_ok == 1:
				merged_file.write(line)
		this_file.close()
	return merged_file_name
	
def subgraph_expansion(metafile, consensusfile):
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
	consensusfile.seek(0)	#In case we've been using the file already
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
	
def predict_interactome(mapping_file_list, metafile, consensusfile):
	cwd = os.getcwd()
	storage_path = "proteomes"
	if not os.path.isdir(storage_path):
		try: 
			os.mkdir(storage_path)
			print("Setting up proteome directory.")
		except OSError:
			if not os.path.isdir(storage_path):
				raise
				
	pred_interactome_path = "predicted_interactomes"
	if not os.path.isdir(pred_interactome_path):
		try: 
			os.mkdir(pred_interactome_path)
			print("Setting up predicted interactome directory.")
		except OSError:
			if not os.path.isdir(pred_interactome_path):
				raise
	
	getting_proteomes = 1			#Can retrieve proteome entries from Uniprot and will map to OGs.
	while getting_proteomes == 1:
		get_new_proteomes = raw_input("Get a proteome from Uniprot? (Y/N)\n")
		if get_new_proteomes in ["Y", "y"]:
			get_a_proteome()	#run get_a_proteome() method
		else:
			print("Will now map proteomes to OGs.")
			break
	
	proteome_list = glob.glob('proteomes\proteome_raw_*.txt')	#Raw proteomes, from Uniprot, in list format, labeled with taxid
	
	map_dict = {}	#Dictionary for Uniprot to OG maps
	
	if len(proteome_list) > 0:	#Only need the OG map if we have raw proteomes to be processed
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
	
	for proteome_filename in proteome_list:					#Map all available raw proteomes to OGs.
		#Proteins without OG mappings retain their Uniprot IDs but we keep track of it in an extra column, too
		print("Mapping proteins in " + proteome_filename)
		proteome_proteins = []
		proteome_map_filename = proteome_filename.replace("raw", "map")
		try:
			proteome_file = open(proteome_filename)
			proteome_map_file = open(proteome_map_filename, "w")
		except IOError as e:
			print("I/O error({0}): {1}".format(e.errno, e.strerror))
		for line in proteome_file:
			one_protein = line.rstrip()
			proteome_proteins.append(one_protein)
		total_proteins = 0
		total_proteins_mapped = 0
		for protein in proteome_proteins:
			og_mapped = 0	#All proteins are unmapped to OGs by default
			total_proteins = total_proteins +1
			if protein in map_dict:
				og_mapped = 1
				total_proteins_mapped = total_proteins_mapped +1
				proteome_map_file.write(map_dict[protein] + "\t" + protein + "\t" + str(og_mapped) + "\n")
			else:
				proteome_map_file.write(protein + "\t" + protein + "\t" + str(og_mapped) + "\n")
		
		proteome_file.close()
		os.remove(proteome_filename)	#Remove the raw file as it's redundant now.
		
		print(proteome_filename + " contains " + str(total_proteins) + " proteins. " 
				+ str(total_proteins_mapped) + " map to OGs.")
		if total_proteins_mapped == 0:
			print("WARNING: No proteins in this proteome map to OGs.")
	
	os.chdir(storage_path)
	proteome_map_list = glob.glob('proteome_map_*.txt')	#Proteomes mapped to eggNOG OGs, labeled with taxid
	os.chdir("..")
	
	#Uses Entrez here for more info about taxid corresponding to proteome.
	print("\nAvailable proteome maps:")
	taxid_context = {}	#We'll keep the taxonomy information for later.
	for proteome_filename in proteome_map_list:
		taxid = ((proteome_filename.split("_"))[2]).rstrip(".txt")
		target_handle = Entrez.efetch(db="Taxonomy", id=str(taxid), retmode="xml")
		target_records = Entrez.read(target_handle)
		#print(target_records)
		taxid_name = target_records[0]["ScientificName"]
		taxid_parent = target_records[0]["ParentTaxId"]
		taxid_division = target_records[0]["Division"]
		if taxid_division != "Bacteria" and useViruses == False:
			print(taxid_name + "\t\t" + proteome_filename + "\tNOTE: Not Bacteria! May not work well with bacterial consensus networks.")
		if taxid_division == "Viruses" and useViruses == True:
			print(taxid_name + "\t\t" + proteome_filename + "\tNOTE: This is a viral proteome. Ensure your meta-interactome uses viral proteins.")
		else:
			print(taxid_name + "\t\t" + proteome_filename)
		taxid_context[taxid] = [taxid_name, taxid_parent, taxid_division]	#This is critical as we'll need it shortly
	
	#Also retrieve taxid details from the taxid context file.
	taxid_context_filenames = glob.glob("taxid_context*.txt")
	if len(taxid_context_filenames) > 1:
		print("More than one taxid context file found. Check for duplicates.")
		return None
	if len(taxid_context_filenames) == 0:
		print("Cannot find taxid context file. Rebuild meta-interactome.")
		return None
	taxid_context_file = open(taxid_context_filenames[0])
	for line in taxid_context_file:
		one_context = (line.rstrip()).split("\t")
		this_taxid = one_context[0]
		taxid_name = one_context[1]
		taxid_parent = one_context[2]
		taxid_division = one_context[3]
		taxid_context[this_taxid] = [taxid_name, taxid_parent, taxid_division]
	
	#Now that we have proteomes, we can use them to predict interactomes.
	
	continue_prediction = raw_input("Predict interactomes for all above? (Y/N)\n")
	if continue_prediction in ["Y", "y"]:
		print("Predicting interactomes for the above species.")
	else:
		return None
		
	print("Loading meta-interactome files.")	#Only uses the consensus right now, but full meta-interactome may be needed
	consensus_interactions = []					#if we want to filter by spoke expansion or know individual proteins
	all_interactions = []
	
	for line in consensusfile:
		one_consensus_interaction = (line.rstrip()).split("\t")
		consensus_interactions.append(one_consensus_interaction)
	consensusfile.close()
	
	for line in metafile:
		one_interaction = (line.rstrip()).split("\t")
		all_interactions.append(one_interaction)	#This is just the raw interaction list at this point
	metafile.close()
	
	#Interactome prediction starts here, iterating through each OG-mapped proteome.
	interactome_stats = {}	#Uses taxid as key
	
	for proteome_filename in proteome_map_list:	#Go through each of the available OG-mapped proteomes
		print("\nPredicting interactome for " + proteome_filename + ".")
		taxid = ((proteome_filename.split("_"))[2]).rstrip(".txt")
		this_proteome_map = {}	#A dictionary of OGs to multiple proteins, since >1 protein may map to an OG.
		this_proteome = []	#A list of just proteins
		this_og_eome = []	#A list of just OGs in the proteome
		this_pred_interactome = []	#Actually the interactome at any one time - the whole prediction is written to file
		this_pred_interactome_detailed = [] #The same interactome, but with contextual details
		#It will also include a prediction category.
		os.chdir(storage_path)
		proteome_map_file = open(proteome_filename)
		for line in proteome_map_file:
			contents = (line.rstrip()).split("\t")
			one_protein = contents[1]
			one_og = contents[0]
			if one_og not in this_proteome_map:
				this_proteome_map[one_og] = [one_protein]
			else:
				this_proteome_map[one_og].append(one_protein)
			this_proteome.append(one_protein)	#Each line in the input should already contain a unique protein ID
			if one_og not in this_og_eome:
				this_og_eome.append(one_og)	#Should be the same as the keys in this_proteome_map
		proteome_map_file.close()
			
		os.chdir("..")
		os.chdir(pred_interactome_path)
		pred_interactome_filename = proteome_filename.replace("proteome_map", "pred_interactome")
		pred_interactome_file = open(pred_interactome_filename, "w")
		
		pred_ppi_count = 0	#The count of PPI from predictions
		exp_ppi_count = 0	#The count of PPI from experimental results, counting spoke expansion
		
		#First pass: check the meta-interactome for exact match PPI
		#This is mostly to account for proteins without OG matches, but we count them all 
		#as we want to distingish between interactions seen already and new predictions.
		#Like with building the consensus set, we need to check for related taxids.
		
		print("Checking for experimental interactions.")
		for interaction in all_interactions:
			same_species = 0	#Well, not the same, but same as the target species OR related
			parent_taxid = taxid_context[taxid][1]
			taxid_A = (((((interaction[9].split("|"))[0]).lstrip("taxid:")).split("("))[0])
			taxid_B = (((((interaction[10].split("|"))[0]).lstrip("taxid:")).split("("))[0])
			if taxid == taxid_A or parent_taxid == taxid_A:	#If taxids are the same as target or its parent
				if taxid == taxid_B or parent_taxid == taxid_B:	
					same_species = 1
			elif taxid == taxid_A or taxid == taxid_context[taxid_A][1]:	#If taxids are child of target
				if taxid == taxid_B or taxid == taxid_context[taxid_B][1]:
					same_species = 1
			elif taxid == taxid_A or parent_taxid == taxid_context[taxid_A][1]:	#If taxids share parent
				if taxid == taxid_B or parent_taxid == taxid_context[taxid_B][1]:
					same_species = 1
			#May throw KeyError here, indicating taxid not in taxid_context - look up if needed?
			if same_species == 1:
				proteinA = interaction[0].lstrip("uniprotkb:")
				proteinB = interaction[1].lstrip("uniprotkb:")
				ogA = interaction[42]
				ogB = interaction[43]
				unique_interaction = [proteinA, proteinB, ogA, ogB]
				unique_interaction_detailed = [proteinA, proteinB, ogA, ogB, "Experimental"]
				if unique_interaction not in this_pred_interactome:
					exp_ppi_count = exp_ppi_count +1
					if exp_ppi_count % 10 == 0:
						sys.stdout.write(".")
					if exp_ppi_count % 100 == 0:
						sys.stdout.write(str(exp_ppi_count))
					this_pred_interactome.append(unique_interaction)
					this_pred_interactome_detailed.append(unique_interaction_detailed)
		sys.stdout.write(str(exp_ppi_count))
		
		print("\nMaking interaction predictions.")
		#Second pass: make predictions based on OGs and the consensus interactome.
		#That is, if two proteins interact, predict all proteins in their two OGs interact.
		#All experimental interactions should be covered in the consensus, so don't care about species here
		#Don't need to handle protein vs. protein as we should have seen it in the meta-interactome already
		
		for interaction in consensus_interactions:
			if interaction[0] in this_og_eome and interaction[1] in this_og_eome:	#Check for OG vs. OG first
				for proteinA in this_proteome_map[interaction[0]]:	#Expand interaction to all possible proteins with OG matches
					for proteinB in this_proteome_map[interaction[1]]:
						unique_interaction = [proteinA, proteinB, interaction[0], interaction[1]]
						unique_interaction_detailed = [proteinA, proteinB, interaction[0], interaction[1], "Predicted"]
						if unique_interaction not in this_pred_interactome:
							pred_ppi_count = pred_ppi_count +1
							this_pred_interactome.append(unique_interaction)
							this_pred_interactome_detailed.append(unique_interaction_detailed)
							if pred_ppi_count % 10 == 0:
								sys.stdout.write(".")
							if pred_ppi_count % 100 == 0:
								sys.stdout.write(str(pred_ppi_count))
			elif interaction[0] in this_og_eome and interaction[1] in this_proteome:	#Check if it's an OG and a protein
				for proteinA in this_proteome_map[interaction[0]]:	#Expand interaction to all possible proteins with OG matches
					proteinB = interaction[1]
					unique_interaction = [proteinA, proteinB, interaction[0], interaction[1]]
					unique_interaction_detailed = [proteinA, proteinB, interaction[0], interaction[1], "Predicted"]
					if unique_interaction not in this_pred_interactome:
						pred_ppi_count = pred_ppi_count +1
						this_pred_interactome.append(unique_interaction)
						this_pred_interactome_detailed.append(unique_interaction_detailed)
						if pred_ppi_count % 10 == 0:
							sys.stdout.write(".")
						if pred_ppi_count % 100 == 0:
							sys.stdout.write(str(pred_ppi_count))
			elif interaction[0] in this_proteome and interaction[1] in interaction[1] in this_og_eome:	#Check if it's a protein and an OG
				proteinA = interaction[0]
				for proteinB in this_proteome_map[interaction[1]]:	#Expand interaction to all possible proteins with OG matches
					unique_interaction = [proteinA, proteinB, interaction[0], interaction[1]]
					unique_interaction_detailed = [proteinA, proteinB, interaction[0], interaction[1], "Predicted"]
					if unique_interaction not in this_pred_interactome:
						pred_ppi_count = pred_ppi_count +1
						this_pred_interactome.append(unique_interaction)
						this_pred_interactome_detailed.append(unique_interaction_detailed)
						if pred_ppi_count % 10 == 0:
							sys.stdout.write(".")
						if pred_ppi_count % 100 == 0:
							sys.stdout.write(str(pred_ppi_count))
		sys.stdout.write(str(pred_ppi_count))
		
		#Finally - get a few more protein and OG counts.
		#These counts won't be right if we just use the proteome, as PPI may include related species
		#Isn't a problem for predictions as they're all based off one proteome
		#But for experimental results we just get every Uniprot ID
		proteins_in_interactions = []
		proteins_w_exp_ppi = []
		proteins_w_pred_ppi = []
		ogs_in_interactions = []
		ogs_w_exp_int = []
		ogs_w_pred_int = []
		exp_og_int = []
		pred_og_int = []
		for interaction in this_pred_interactome_detailed:
			og_pair = interaction[2:4] 
			rev_og_pair = [interaction[3], interaction[2]]	#We don't care about interaction direction.
			if interaction[4] == "Experimental":
				for protein in interaction[0:2]:
					if protein not in proteins_w_exp_ppi:
						proteins_w_exp_ppi.append(protein)
				for og in og_pair:
					if og not in ogs_w_exp_int:
						ogs_w_exp_int.append(og)
				if og_pair not in exp_og_int and rev_og_pair not in exp_og_int:
					exp_og_int.append(og_pair)
			for protein in this_proteome:
				both_interactors = 0
				if protein in interaction[0:2] and protein not in proteins_in_interactions:
					both_interactors = both_interactors +1
					proteins_in_interactions.append(protein)
					if interaction[4] == "Predicted" and protein not in proteins_w_pred_ppi:
						proteins_w_pred_ppi.append(protein)
				if both_interactors == 2:
					break
			for og in this_og_eome:
				both_interactors = 0
				if og in interaction[2:4] and og not in ogs_in_interactions:
					both_interactors = both_interactors +1
					ogs_in_interactions.append(og)
					if interaction[4] == "Predicted" and og not in ogs_w_pred_int:
						ogs_w_pred_int.append(og)
				if both_interactors == 2:
					break
			if og_pair not in exp_og_int and rev_og_pair not in pred_og_int:
				pred_og_int.append(og_pair)
					
		proteins_not_in_interactions = len(this_proteome) - len(proteins_in_interactions)	#Just a count, here
		ogs_not_in_interactions = len(this_og_eome) - len(ogs_in_interactions)
		#interactome_stats contains statistics used in batch output. Contains:
		#Name, taxid, Proteins, ProteinsNotInPPI, ProteinsWithExpPPI, ProteinsWithPredPPI, 
		#UniqueOGs, OGsWithoutInteractions, OGsWithExpInt, OGsWithPredInt, ExpOGIntNet, OGIntInPredNet
		interactome_stats[taxid] = [taxid_context[taxid][0], taxid, str(len(this_proteome)), str(proteins_not_in_interactions), 
									str(len(proteins_w_exp_ppi)), str(len(proteins_w_pred_ppi)), str(len(this_og_eome)),
									str(ogs_not_in_interactions), str(len(ogs_w_exp_int)), str(len(ogs_w_pred_int)),
									str(len(exp_og_int)), str(len(pred_og_int))]
		
		#This is just for testing.
		'''
		for interaction in all_interactions:
			taxid_A = (((((interaction[9].split("|"))[0]).lstrip("taxid:")).split("("))[0])
			taxid_B = (((((interaction[10].split("|"))[0]).lstrip("taxid:")).split("("))[0])
			if taxid == taxid_A or taxid == taxid_B:
				proteinA = interaction[0].lstrip("uniprotkb:")
				proteinB = interaction[1].lstrip("uniprotkb:")
				ogA = interaction[42]
				ogB = interaction[43]
				this_meta_interaction = [proteinA, proteinB, ogA, ogB]
				if this_meta_interaction not in this_pred_interactome:
					print(this_meta_interaction)
		'''
			
		#Write interactome file.
		for interaction in this_pred_interactome_detailed:
			pred_interactome_file.write("\t".join(interaction) + "\n")
		
		print("\nFound " + str(exp_ppi_count)  + " experimental interactions (including spoke" + 
				" expansion) and made " + str(pred_ppi_count) + " interaction predictions" +
				" for " + taxid_context[taxid][0] + ".")
		
		pred_interactome_file.close()
		os.chdir("..")
		
	#Once all the interactome predictions for all proteomes are done, do summary statistics.
	nowstring = (date.today()).isoformat()
	multi_inter_stats_file_name = "interactome_statistics_" + nowstring + ".txt"
	multi_inter_stats_file = open(multi_inter_stats_file_name, "w")
	os.chdir("predicted_interactomes")
	interactome_filenames = glob.glob("pred_interactome_*.txt")
	
	stats_outlines = []	#This is where the output for each species will go, one interactome per line
	for filename in interactome_filenames:
		taxid = ((filename.rstrip(".txt")).split("_"))[2]
		outline = "\t".join(interactome_stats[taxid]) + "\n"
		stats_outlines.append(outline)
		
	os.chdir("..")
	#Text header
	multi_inter_stats_file.write("Species\tTaxid\tProteins\tProteinsNotInPPI\tProteinsWithExpPPI\tProteinsWithPredPPI\t" +
								"UniqueOGs\tOGsWithoutInteractions\tOGsWithExpInt\tOGsWithPredInt\tOGIntInExpNet\tOGIntInPredNet\n")
	for outline in stats_outlines:
		multi_inter_stats_file.write(outline)
	
	print("\nWrote summary statistics for these interactomes to " + multi_inter_stats_file_name)
	print("\nComplete.\n")
	
def get_a_proteome():	#Does what it says.	Much more organized than the rest of this since I wrote it a while ago.
	
	def get_search_url(query, fil):
		search_url = "http://www.uniprot.org/proteomes/?query=" + query + \
					"&fil=" + fil + "&sort=score"
		return search_url
	
	def parse_search(up_input):
		search_results = []
		soup = BeautifulSoup(up_input)
		for child in (soup.find_all('tr')):
			single_result = child.get_text("\t")
			search_results.append(single_result)
		del search_results[0:2]
		if len(search_results) == 0:
			print("No results found.")
			return None
		return search_results
	
	def get_proteome_url(entry, format_choice):
		proteome_url = "http://www.uniprot.org/uniprot/?sort=&desc=&query=proteome:" + entry + "&force=no&format=" + format_choice
		return proteome_url
		
	def parse_proteome_entry(up_input):
		if not up_input:
			entry_text = "EMPTY"
		else:
			soup = BeautifulSoup(up_input)
			entry_text = (soup.p.get_text())
		return entry_text
		
	def save_proteome(text,taxid):
		os.chdir("proteomes")
		filename = "proteome_raw_" + str(taxid) + ".txt"
		try:
			outfile = open(filename, 'wb')
		except IOError as e:
			print("I/O error({0}): {1}".format(e.errno, e.strerror))
			sys.exit()
		for line in text:
			outfile.write(line)
		print("File written to " + filename)
		outfile.close()
		os.chdir("..")

	#Retrieve proteomes on a query
	query = (raw_input("Please specify a full or partial species name.\n")).rstrip()
	ref_filter = "reference%3Ayes"
	if useNonRefProteomes == True:
		ref_filter = ""
	search_results_url = get_search_url(query, ref_filter) #Leave filter as "" to get non-reference proteomes too
								#Other option: taxonomy%3A"Bacteria+%5B2%5D" for just bacteria
	
	search_response = requests.get(search_results_url)
	
	#Output the query results
	print(search_response)
	proteome_entries = parse_search(search_response.text)
	if proteome_entries == None:
		return None
	i = 0
	print("Result\tAccession\tName")
	for entry in proteome_entries:
		print(str(i) + "\t" + entry)
		i = i +1
	
	#Choose a single proteome and output to file
	choice = raw_input('Please choose a search result.\n')
	if not re.match("^[0-9]*$", choice):
		print("Numbers only, please.")
		sys.exit()
	chosen_entry = (proteome_entries[int(choice)]).split("\t")
	print("Retrieving proteome for " + chosen_entry[1])
	proteome_url = get_proteome_url(chosen_entry[0], "list") #Options include list, txt, tab
	proteome_response = requests.get(proteome_url)
	proteome_text = parse_proteome_entry(proteome_response.text)
	if proteome_text == "EMPTY":
		print("Could not retrieve this proteome. See the Uniprot entry for %s." % chosen_entry[0])
	else:
		save_proteome(proteome_text, chosen_entry[2])

def describe_consensus(consensusfile):
	cons_stats_filenames = glob.glob("cons_statistics_*.txt")
	if len(cons_stats_filenames) > 1:
		print("More than one consensus statistics file found. Check for duplicates.")
		return None
	if len(cons_stats_filenames) == 0: 
		print("No consensus statistics file found. Will skip basic counts.")
	else:
		con_stats = open(cons_stats_filenames[0])
		for line in con_stats:
			print(line)
		print("\n")
	
	consensus_interactions = []
	for line in consensusfile:
		one_interaction = (line.rstrip()).split("\t")
		consensus_interactions.append(one_interaction)
		
	taxids_and_context = {}
	taxid_ref_list = glob.glob('taxid_context*.txt')
	if len(taxid_ref_list) >1:
		sys.exit("Something went wrong - more than one taxid context file found.")
	if len(taxid_ref_list) == 0:
		sys.exit("Something went wrong - no taxid context file found.")
	taxid_ref_file = open(taxid_ref_list[0])
	for line in taxid_ref_file:
		content = ((line.rstrip()).split("\t"))
		taxids_and_context[content[0]] = [content[1], content[2], content[3]]
	taxid_ref_file.close()
	
	print("Top taxid contributions, in number of consensus interactions corresponding to the taxid.")
	print("Name\tTaxid\tNumber of interactions")
	all_taxids = {}	#All taxids AND their counts.
	
	for interaction in consensus_interactions:	#Check each interaction for contributing taxids
		these_sources = interaction[4].split()
		for taxid in these_sources:
			if taxid not in all_taxids:
				all_taxids[taxid] = 1
			all_taxids[taxid] = all_taxids[taxid] + 1

	sorted_taxids = sorted(all_taxids.items(), key=operator.itemgetter(1), reverse=True)
	top_taxids = sorted_taxids[0:15]
	
	for taxid in top_taxids:
		taxid_only = taxid[0]
		taxid_name = taxids_and_context[taxid_only][0]
		print(taxid_name + "\t" + taxid_only + "\t" + str(taxid[1]))
			
#Main

#Check for eggNOG mapping file and get if needed
#Requires downloading several files and building new mapping file from them
mapping_file_list = glob.glob('uniprot_og_maps*.txt')
if len(mapping_file_list) >2:
	sys.exit("Found more than one mapping file. Check for duplicates.")
if len(mapping_file_list) == 0:
	print("No eggNOG mapping files found or they're incomplete. Rebuilding them.")
	get_eggnog_maps()
	mapping_file_list = glob.glob('uniprot_og_maps*.txt')
	
#Check for eggNOG annotation file and get if needed
annotation_file_list = glob.glob('*annotations.tsv')
expected_filecount = 2
if useViruses == True:
	expected_filecount = 3
if len(annotation_file_list) > expected_filecount:
	sys.exit("Found more eggNOG annotation files than expected. Check for duplicates.")
if len(annotation_file_list) < expected_filecount:
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
	print("\nNo meta-interactome found.")
	while ppi_data_filename == "":
		ppi_data_option = raw_input("Retreive IntAct bacterial PPI or use local file(s) to build meta-interactome?\n"
		"Enter:\n R for retrieval\n L for local file, \n M for multiple inputs, \n or X to quit.\n")
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
			print("Will use single local file. Note that it should be in PSI-MI TAB 2.7 format and have no header row.")
			ppi_data_filename = raw_input("Please provide local filename.\n")
			interaction_file_list = glob.glob(ppi_data_filename)
			if len(interaction_file_list) == 0:
				sys.exit("Can't find a file with that filename.")	
		if ppi_data_option in ["M", "m"]:	#Uses multiple local files in PSI-MI Tab27 format
			print("Will append multiple local files. Note that each should be in PSI-MI TAB 2.7 format and have no header row.")
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
		if ppi_data_option in ["X", "x"]:
			sys.exit("Exiting...")

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
	request_next = raw_input("\nChoose from the following options.\n" 
		"A: Generate expanded subgraphs of the consensus network, filtering by function.\n"
		"B: Generate a predicted interactome for one or more proteomes.\n"
		"C: Get statistics for the consensus meta-interactome.\n"
		"X: Exit.\n") 
	if request_next in ["x", "X"]:
		sys.exit("Exiting...")
	if request_next in ["a", "A"]:
		subgraph_expansion(metafile, consensusfile)
	if request_next in ["b", "B"]:
		predict_interactome(mapping_file_list, metafile, consensusfile)
	if request_next in ["c", "C"]:
		describe_consensus(consensusfile)
	print("\nChoose from the list, please.")
			
sys.exit(0)
