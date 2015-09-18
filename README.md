# network-umbra.py

Predicts interactions in a protein interaction network based off a meta-interactome network.
Uses eggNOG v.4.1.

REQUIRES: 
Biopython 1.65 or more recent
Also needs ~500 MB of available disk space to accomodate data files amd output
More space may be necessary for proteome files.

INPUT: 

Downloads all available protein-protein interactions for bacteria from IntAct.
Downloads highest-level (LUCA) and bacteria-specific Uniprot ID to NOG mappings from eggNOG v.4.1.
Downloads highest-level (LUCA) bacteria-specific NOG annotations from eggNOG v.4.1.

OUTPUT: 

'metainteractome[date].txt'
			A meta-interactome composed of all available bacterial protein-protein interactions.
			Follows PSI-MI Tab27 format, with the addition of two ortholog identifiers per row.
			See format description at https://code.google.com/p/psimi/wiki/PsimiTab27Format
			
'meta_statistics[date].txt'
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
			FuncCatA		Functional category of the first interactor – see FuncCats tab
			DescA		Description of the first interactor
			FuncCatB		Functional category of the second interactor – see FuncCats tab
			DescB		Description of the second interactor

'cons_statistics[date].txt'
			Contains statistics relevant to the produced consensus meta-interactome.

At the moment, this only makes predictions based off presence of the same OGs as in the consensus network.
It needs to verify that both OGs in the predicted PPI are present in the target species.
Redundant predictions (the same interaction from the same taxon ID) are merged.

CHANGES COMPLETE:

*Downloads eggNOG map file (LUCA-level and bacteria specific) and IntAct interactions (just bacteria specific)
*Generates meta-interactome and rudimentary consensus meta-interactome.
*IntAct data cleaned before using (removes "intact" and "chebi" interactors)
*A few basic counts (interactors and interactions) are made for meta-interactome and consensus sets
*Counts for all consensus interactions are also made across the whole meta-interactome and provided in consensus network
*Downloads the eggNOG annotation file for all NOGs but doesn't do anything with it yet
*Gets taxon IDs, names, and parent taxon IDs. Adds them to interactions in consensus network but doesn't compare to eliminate redundant taxids
	Checks for parent and child relationships between taxon IDs to limit redundancy.
*Gets and maps FuncCat and description annotations (for both LUCA-level and bacteria) to OGs. Use them in the consensus network. 

IN PROGRESS:
Priorities:

**Some non-bacterial proteins are present within PPI in the input interaction set, or at least I found taxids for humans in the consensus. Check on why.
**Verify that the taxids in the consensus really correspond to the interaction.
	All taxids should have both interactors in their genomes.
**Filter by FuncCat and produce subsets.
**Assign methods to interactions (more general than original data, so we can detect spoke expansion)

*Get counts and statistics for input data and various interactomes.
*Do interactome prediction for a given proteome.
*Use protein and species count from eggNOG (it's in the annotation file).
*Output interaction sets, filtered by FuncCat (and especially OG UFs).
*Perform ANOVA between different FuncCats to see consensus interaction patterns.
*Download a proteome with a search query and set up OG mapping for it.

