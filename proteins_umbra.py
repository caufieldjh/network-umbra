#!/usr/bin/python
#proteins_umbra.py
'''
Downloads a reference proteome and assigns an orthologous group (OG) to each.
Uses eggNOG v.4.1 - or whatever the most recent version is.

REQUIRES: Biopython 1.65 or more recent
			Needs ~5 GB of disk space for ID conversion files.
			Needs an additional ~29 GB for mapping virus protein IDs.
			 (This is because it uses the full Uniprot ID mapping database,
			 which is excessive but more reliable than their mapping server
			 for large mapping requests)

INPUT: Downloads a reference proteome from Uniprot.
		Produces OG maps if needed or uses those produced by Network_umbra.py.

OUTPUT: 
'proteome_map_[taxid].txt' - on each line:
  a single OG membership
  the UniProtAC of the protein
  a binary value indicating whether the protein maps to an OG (0 if no, 1 if yes)
  
'''

import glob, gzip, operator, os, re, requests, sys, urllib2, zipfile
from Bio import Entrez
from bs4 import BeautifulSoup
from collections import Counter
from datetime import date

Entrez.email = ''

#Options
useViruses = False	#Option for using eggNOG's viral OGs. Requires the filters permitting only Bacteria to be modified
					#Also requires the viral OGs to be downloaded and added.
					#This option needs to be set True BEFORE the Uniprot to OG map is built or it won't include proteins from viruses
#NOTE: Viruses are not currently in the eggNOG ID conversion file
#The eggNOG protein IDs vary from protein to protein but are often Uniprot IDs (mnemonic identifiers, i.e. A9J730_BPLUZ)
#We still check the ID conversion file for them in case it gets updated soon
					
useNonRefProteomes = True	#Option to search non-reference Uniprot proteomes
#Many of these proteomes have been made redundant in Uniprot
#and this script ignores redundant results, so they will not be seen

#Functions

def chunkit(input_seq, chunk_size):
    return (input_seq[position:position + chunk_size] for position in xrange(0, len(input_seq), chunk_size))

def get_eggnog_maps(version): 
	#Download and unzip the eggNOG ID conversion file 
	#Filters file to just Uniprot IDs; the resulting file is the map file.
	#One Uniprot ID may correspond to multiple OGs - e.g. COG1234,COG3810,COG9313. 
	#these cases are considered OGs in their own right as this may indicate a pattern of conserved sequences on its own 
	baseURL = "http://eggnogdb.embl.de/download/" + version + "/"
	convfilename = "eggnog4.protein_id_conversion.tsv.gz"	#File contains ALL database identifiers and corresponding proteins
	
	convfilepath = baseURL + convfilename
	outfilepath = convfilename[0:-3]
	dl_convfile = 1	#If 1, we need to download
	if os.path.isfile(convfilename): #Already have the compressed file, don't download
		print("Found compressed ID conversion file on disk: %s" % convfilename)
		decompress_convfile = 1
		dl_convfile = 0
	if os.path.isfile(outfilepath): #Already have the decompressed file, don't download
		print("Found ID conversion file on disk: %s" % outfilepath)
		decompress_convfile = 0
		dl_convfile = 0
	
	if dl_convfile == 1:
		print("Downloading ID mapping file - this file is large so this may take some time.")
		print("Downloading from %s" % convfilepath)
		response = urllib2.urlopen(convfilepath)
		filesize = response.info()['Content-Length']
		compressed_file = open(os.path.basename(convfilename), "w+b") #Start local compressed file
		chunk = 2097152
		totaldata = 0
		while 1:
			data = (response.read(chunk)) #Read two Mb at a time
			compressed_file.write(data)
			totaldata = totaldata + chunk
			if not data:
				print("\n%s file download complete." % convfilename)
				compressed_file.close()
				break
			sys.stdout.flush()
			sys.stdout.write("\r%s out of %s bytes" % (totaldata, filesize))
		decompress_convfile = 1
		
	if decompress_convfile == 1:
		print("Decompressing map file. Lines written, in millions:")
		with gzip.open(convfilename) as infile: #Open that compressed file, read and write to uncompressed file
			outfile = open(outfilepath, "w+b")
			linecount = 0
			for line in infile:
				outfile.write(line)
				linecount = linecount +1
				if linecount % 1000000 == 0:
						sys.stdout.flush()
						sys.stdout.write("\r%s" % (linecount/1000000))
			infile.close()
		newconvfilename = outfilepath
		outfile.close()
	
	#Download and decompress member NOG files (at least 2 of them)
	nogURL = baseURL + "data/NOG/"
	nogfilename = "NOG.members.tsv.gz"
	bactnogURL = baseURL + "data/bactNOG/"
	bactnogfilename = "bactNOG.members.tsv.gz" 
	all_nog_locations = [[nogURL, nogfilename], [bactnogURL, bactnogfilename]]
	
	if useViruses == True:	#Need some additional files to handle viral proteins
		virnogURL = baseURL + "data/viruses/Viruses/"
		virnogfilename = "Viruses.members.tsv.gz"
		all_nog_locations.append([virnogURL, virnogfilename])
		up_baseURL = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/"
		up_mapping_filename = "idmapping.dat.gz"
		up_mapping_filepath = up_baseURL + up_mapping_filename
		up_outfilepath = up_mapping_filename[0:-3]
		
		dl_up_mapping_file = 1	#If 1, we need to download
		if os.path.isfile(up_mapping_filename): #Already have the compressed file, don't download
			print("Found compressed Uniprot ID conversion file on disk: %s" % up_mapping_filename)
			decompress_up_mapping_file = 1
			dl_up_mapping_file = 0
		if os.path.isfile(up_outfilepath): #Already have the decompressed file, don't download
			print("Found ID conversion file on disk: %s" % up_outfilepath)
			decompress_up_mapping_file = 0
			dl_up_mapping_file = 0
		
		if dl_up_mapping_file == 1:
			print("\nDownloading Uniprot ID mapping file for viral protein mapping. Please wait as this file is large.")
			print("Downloading from %s" % up_mapping_filepath)
			response = urllib2.urlopen(up_mapping_filepath)
			filesize = response.info()['Content-Length']
			compressed_file = open(os.path.basename(up_mapping_filename), "w+b") #Start local compressed file
			chunk = 2097152
			totaldata = 0
			while 1:
				data = (response.read(chunk)) #Read two Mb at a time
				compressed_file.write(data)
				if not data:
					print("\n%s file download complete." % up_mapping_filename)
					compressed_file.close()
					break
				sys.stdout.flush()
				sys.stdout.write("\r%s out of %s bytes" % (totaldata, filesize))
			decompress_up_mapping_file = 1
			
		if decompress_up_mapping_file == 1:
			print("Decompressing Uniprot ID mapping file. Lines written, in millions:")
			with gzip.open(up_mapping_filename) as infile: #Open that compressed file, read and write to uncompressed file
				outfile = open(up_outfilepath, "w+b")
				linecount = 0
				for line in infile:
					outfile.write(line)
					linecount = linecount +1
					if linecount % 1000000 == 0:
							sys.stdout.flush()
							sys.stdout.write("\r%s" % (linecount/1000000))
				infile.close()
				outfile.close()
	
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
			filesize = response.info()['Content-Length']
			compressed_file = open(os.path.basename(memberfilename), "w+b") #Start local compressed file
			chunk = 2097152
			totaldata = 0
			while 1:
				data = (response.read(chunk)) #Read two Mb at a time
				compressed_file.write(data)
				if not data:
					print("\n%s file download complete." % memberfilename)
					compressed_file.close()
					break
				sys.stdout.flush()
				sys.stdout.write("\r%s out of %s bytes" % (totaldata, filesize))
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
					if linecount % 1000000 == 0:
						sys.stdout.flush()
						sys.stdout.write("\r%s" % (linecount/1000000))
				infile.close()
			outfile.close()
			
	#Clean up by removing compressed files
	print("\nRemoving compressed files.")
	all_compressed_files = [convfilename, nogfilename, bactnogfilename]
	if useViruses == True:
		for this_filename in [virnogfilename, up_mapping_filename]:
			all_compressed_files.append(this_filename)
	for filename in all_compressed_files:
		if os.path.isfile(filename):
			os.remove(filename)
	
	#Load and filter the ID conversion file as dictionary
	print("Parsing ID conversion file. Lines read, in millions:")
	with open(convfilename[0:-3]) as infile:
		id_dict = {}	#Dictionary of eggNOG protein IDs as values and database IDs (UniprotAC) as keys
		#Gets filtered down to relevant database IDs (i.e., Uniprot IDs)
		linecount = 0
		for line in infile:
			linecount = linecount +1
			line_raw = ((line.rstrip()).split("\t"))	
			one_id_set = [line_raw[0] + "." + line_raw[1], line_raw[2], line_raw[3]] #Protein IDs are split for some reason; merge them
			if "UniProt_AC" in one_id_set[2]:
				id_dict[one_id_set[1]] = one_id_set[0]
			if linecount % 1000000 == 0:
				sys.stdout.flush()
				sys.stdout.write("\r%s" % (linecount/1000000))
		infile.close()

	#Use filtered ID conversion input to map to NOG members
	print("\nReading NOG membership files.")
	all_nog_filenames = [nogfilename[0:-3], bactnogfilename[0:-3]]
	if useViruses == True:
		all_nog_filenames.append(virnogfilename[0:-3])
	nog_members = {}	#Dictionary of NOG ids with protein IDs as keys (need to split entries for each)
	nog_count = 0
	for filename in all_nog_filenames:
		temp_nog_members = {}	#We will have duplicates within each set but don't want to lose the information.
		print("Reading from %s" % filename)
		with open(filename) as infile:
			membercol = 5	#The column where the NOG members are
			if filename == virnogfilename[0:-3]:	#The virus members file has a different format as there is no FuncCat column
				infile.readline()	#Skip the first line
				membercol = 4
				viral_ids = [] #A list of viral eggNOG protein IDs, some of which are Uniprot IDs to be converted to ACs
			for line in infile:
				nog_count = nog_count +1
				line_raw = ((line.rstrip()).split("\t"))
				nog_id = line_raw[1]
				line_members = line_raw[membercol].split(",")
				for protein_id in line_members:			
					if filename == virnogfilename[0:-3]: #If Viruses, we need to convert IDs as they aren't in the eggNOG ID conversion file.
						if protein_id not in viral_ids:
							viral_ids.append(protein_id)
					if protein_id in temp_nog_members: #The same protein could be in more than one OG at the same level
						temp_nog_members[protein_id] = temp_nog_members[protein_id] + "," + nog_id
					else:
						temp_nog_members[protein_id] = nog_id
			infile.close()
		nog_members.update(temp_nog_members)
		
	
	if useViruses == True:
		
		#We use three different dictionaries here.
		#The first is Uniprot IDs to UniprotACs (just for viral proteins)
		#The second is eggNOG IDs to Uniprot IDs.
		#The third is UniprotACs to eggNOG IDs - this is id_dict{} already.
		uniprotID_to_uniprotAC = {}
		eggnog_to_uniprotID = {}
		unmapped_ids = [] #eggNOG protein IDs which may not contain Uniprot IDs
		
		#We go through the viral protein IDs twice, first to get Uniprot IDs
		#and then to add them to id_dict.
		
		for viral_id in viral_ids:
			eggnog_to_uniprotID[viral_id] = (viral_id.split("."))[1]	#Remove the taxid
			#Some of the eggNOG IDs may not include UniprodIDs, but many do
		
		#This data file is too large to efficiently much of it in a dict.
		#Luckily we just got the IDs we need here to filter it
		print("Parsing Uniprot ID mapping file. Lines read, in millions:")
		with open(up_outfilepath) as infile:
			linecount = 0
			for line in infile:
				linecount = linecount +1
				line_raw = ((line.rstrip()).split("\t"))	
				if line_raw[1] == "UniProtKB-ID" and line_raw[2] in eggnog_to_uniprotID.values():
					uniprotID_to_uniprotAC[line_raw[2]] = line_raw[0]
				if linecount % 1000000 == 0:
					sys.stdout.flush()
					sys.stdout.write("\r%s" % (linecount/1000000))
			infile.close()
		
		print("Finding identifiers for %s viral proteins." % len(viral_ids))
		
		for viral_id in viral_ids:
			upid = eggnog_to_uniprotID[viral_id]
			if upid in uniprotID_to_uniprotAC:
				upid_ac = uniprotID_to_uniprotAC[upid]
				id_dict[upid_ac] = viral_id
			else:
				unmapped_ids.append(viral_id)
				
		print("Done mapping viral proteins.")
		print("The following entries were not recognized as Uniprot IDs:")
		print(unmapped_ids)
	
	#Get counts of how many identifiers we have now
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
				sys.stdout.flush()
				sys.stdout.write(".")
			if mapped_count % 1000000 == 0:
				sys.stdout.flush()
				sys.stdout.write(str(mapped_count/1000000))
		
	#Use this mapping to build map file, named "uniprot_og_maps_*.txt"
	print("\nWriting map file.")
	nowstring = (date.today()).isoformat()
	mapfilename = "uniprot_og_maps_" + nowstring + ".txt"
	mapfile = open(mapfilename, "w+b")
	for mapping in upid_to_NOG:
		mapfile.write(mapping + "\t" + upid_to_NOG[mapping] + "\n")	#Each line is a uniprot ID and an OG id
	mapfile.close() 

def get_eggnog_annotations():
	#Downloads and extracts the eggNOG NOG annotations. 
	baseURLs = ["http://eggnogdb.embl.de/download/latest/data/bactNOG/", "http://eggnogdb.embl.de/download/latest/data/NOG/"]
	bactannfilename = "bactNOG.annotations.tsv.gz"	#The annotations for bacteria-specific NOGs
	lucaannfilename = "NOG.annotations.tsv.gz"	#The annotations for other NOGs, but not bacteria-specific NOGs
	annfilenames = [bactannfilename, lucaannfilename]
	
	if useViruses == True:
		virannfilename = "Viruses.annotations.tsv.gz"
		baseURLs.append("http://eggnogdb.embl.de/download/latest/data/viruses/Viruses/")
		annfilenames.append(virannfilename)
	
	this_url = 0
	for annfilename in annfilenames:
		annfilepath = baseURLs[this_url] + annfilename
		this_url = this_url +1
		outfilepath = annfilename[0:-3]
		if os.path.isfile(annfilename): 
			print("Found compressed annotation file on disk: " + annfilename)
		else:
			response = urllib2.urlopen(annfilepath)
			filesize = response.info()['Content-Length']
			print("Downloading from " + annfilepath)
			compressed_file = open(os.path.basename(annfilename), "w+b") #Start local compressed file
			chunk = 2097152
			totaldata = 0
			while 1:
				data = (response.read(chunk)) #Read two Mb at a time
				compressed_file.write(data)
				if not data:
					print("\n" + annfilename + " file download complete.")
					compressed_file.close()
					break
				sys.stdout.flush()
				sys.stdout.write("\r%s out of %s bytes" % (totaldata, filesize))
		
		print("Decompressing annotation file.")
		with gzip.open(annfilename) as infile: #Open that compressed file, read and write to uncompressed file
			file_content = infile.read()
			outfile = open(outfilepath, "w+b")
			outfile.write(file_content)
			infile.close()
		outfile.close()
		
	print("\nRemoving compressed files.")
	all_compressed_files = [bactannfilename, lucaannfilename]
	if useViruses == True:
		all_compressed_files.append(virannfilename)
	for filename in all_compressed_files:
		os.remove(filename)
		
def get_mapped_proteome(mapping_file_list, custom):
	cwd = os.getcwd()
	storage_path = "proteomes"
	proteome_list = []
	
	if custom:
		custom_list = raw_input("Please provide the name of a file"
								" containing one UniprotAC per line.\n")
		proteome_list.append(custom_list)
		
	if not os.path.isdir(storage_path):
		try: 
			os.mkdir(storage_path)
			print("Setting up proteome directory.")
		except OSError:
			if not os.path.isdir(storage_path):
				raise
	
	if not custom:
		getting_proteomes = 1			#Can retrieve proteome entries from Uniprot and will map to OGs.
		while getting_proteomes == 1:
			get_new_proteomes = raw_input("Get a proteome from Uniprot? (Y/N)\n")
			if get_new_proteomes in ["Y", "y"]:
				get_a_proteome()	#run get_a_proteome() method
			else:
				print("Will now map proteomes to OGs.")
				break
		
		os.chdir(storage_path)
		proteome_list = glob.glob('proteome_raw_*.txt')	#Raw proteomes, from Uniprot, in list format, labeled with taxid
		os.chdir("..")
	
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
	
	unmapped_taxids = [] #These are the taxids without any OG mapping
	for proteome_filename in proteome_list:					#Map all available raw proteomes to OGs.
		#Proteins without OG mappings retain their Uniprot IDs but we keep track of it in an extra column, too
		
		print("Mapping proteins in " + proteome_filename)
		
		if not custom:
			taxid = ((proteome_filename.split("_"))[2]).rstrip(".txt")
			proteome_map_filename = proteome_filename.replace("raw", "map")
			os.chdir(storage_path)
		else:
			proteome_map_filename = ("%s%s%s" % (proteome_filename[:-4],
										"_map", proteome_filename[-4:]))
			
		proteome_proteins = []
		
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
			unmapped_taxids.append(taxid)
		os.chdir("..")
	
	if not custom:
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
				nameline = (taxid_name + "\t\t" + proteome_filename + "\tNOTE: Not Bacteria! May not work well with bacterial consensus networks.")
			if taxid_division == "Viruses" and useViruses == True:
				nameline = (taxid_name + "\t\t" + proteome_filename + "\tNOTE: This is a viral proteome. Ensure your meta-interactome uses viral proteins.")
			else:
				nameline = (taxid_name + "\t\t" + proteome_filename)
			if taxid in unmapped_taxids:
				nameline = "** " + nameline
			print(nameline)
			taxid_context[taxid] = [taxid_name, taxid_parent, taxid_division]	#Not used at the moment
			
		if len(unmapped_taxids) > 0:
			print("Maps marked with ** have no OG mappings.")
	
	print("\nComplete.\n")
	
def get_a_proteome():	#Does what it says.	Much more organized than the rest of this since I wrote it a while ago.
	
	def get_search_url(query, fil):
		search_url = "http://www.uniprot.org/proteomes/?query=" + query + \
					"+redundant%3Ano&fil=" + fil + "&sort=score"
		return search_url
	
	def parse_search(up_input):
		search_results = []
		soup = BeautifulSoup(up_input, "lxml")
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
			soup = BeautifulSoup(up_input, "lxml")
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
	#print(search_response)
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
		print("Could not retrieve this proteome. It may be a redundant entry. See the Uniprot entry for %s." % chosen_entry[0])
	else:
		save_proteome(proteome_text, chosen_entry[2])
		
def main():
	#Check for eggNOG mapping file and get if needed
	#Requires downloading several files and building new mapping file from them
	mapping_file_list = glob.glob('uniprot_og_maps*.txt')
	if len(mapping_file_list) >2:
		sys.exit("Found more than one mapping file. Check for duplicates.")
	if len(mapping_file_list) == 0:
		print("No eggNOG mapping files found or they're incomplete. Rebuilding them.")
		version = "latest"
		version = raw_input("Which eggNOG version would you prefer? Default is latest version.\n")
		if version not in ["4.5","4.1","4.0"]:
			version = "latest"
		else:
			version = "eggnog_" + version
		get_eggnog_maps(version)
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
		
	#Quit now or ask for next step.
	requested = 0
	while requested == 0:
		request_next = raw_input("\nChoose from the following options.\n" 
			"A: Download a reference proteome and map to OGs.\n"
			"B: Provide a list of UniprotAC IDs to map to OGs.\n"
			"X: Exit.\n") 
		if request_next in ["x", "X"]:
			sys.exit("Exiting...")
		if request_next in ["a", "A"]:
			custom = False
			get_mapped_proteome(mapping_file_list, custom)
		if request_next in ["b", "B"]:
			custom = True
			get_mapped_proteome(mapping_file_list, custom)
	
		print("\nChoose from the list, please.")


if __name__ == "__main__":
	sys.exit(main())
