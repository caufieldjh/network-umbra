# network-umbra

Predicts protein interactions using gene orthology and a combined consensus network.

INPUT: 

			'*consensus*.txt' - a tab-delimited file for the consensus network (interactions between OGs)
			
			Each line includes five strings, in this order:

			an OG, 

			a species ID or a pair of IDs (format, "197 vs. 192222")

			a species name, ideally with a strain identifier, 

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

			

OUTPUT: 

			'predicted_interactions_[taxid].txt' - a file of the predicted interactions for the target species

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

