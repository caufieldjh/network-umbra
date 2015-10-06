# network-umbra.py

Predicts interactions in a protein interaction network based off a meta-interactome network.
Uses eggNOG v.4.1.

**REQUIRES**: Biopython 1.65 or more recent. Also needs ~600 MB of available disk space to accomodate data files and output. More space may be necessary for proteome files.

**INPUT**: Downloads all available protein-protein interactions for bacteria from IntAct.
  Alternatively, uses a provided PPI data file in PSI-MI TAB 2.7 format.
  REMOVE THE HEADER ROW if it's present!
  Downloads highest-level (LUCA) and bacteria-specific Uniprot ID to NOG mappings from eggNOG v.4.1.
  Downloads highest-level (LUCA) bacteria-specific NOG annotations from eggNOG v.4.1.

**OUTPUT**: 

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
   
1. InteractorA  The first interactor. Usually an OG.
2. InteractorB  The second interactor. Usually an OG.
3. InteractionCount  Count of individual PROTEIN interactions contributing to this consensus interaction, as per the meta-interactome.
4. TaxonCount  Count of different taxons (here, a proxy for species) corresponding to the interaction.
   Similar taxons have been grouped together where possible, e.g. two different E. coli K-12 strains are just considered E. coli K-12.
5. Taxons  The taxons corresponding to this interaction.
6. FuncCatA  Functional category of the first interactor
7. DescA  Description of the first interactor
8. FuncCatB  Functional category of the second interactor
9. DescB  Description of the second interactor

'cons_statistics[date].txt'
   Contains statistics relevant to the produced consensus meta-interactome.
   
'subgraph_expansion_[FuncCat]_[date].txt'
   A set of subgraphs of the consensus meta-interactome, filtered by conservation of interactions and function of interactors.Each line is one interaction between a consensus interactor and a unique protein, accompanied by the source taxid of the unique protein.
   
'subgraph_expansion_[FuncCat]_nodes_[date].txt'
   Annotation file for the nodes in the expanded subgraphs.
   
Outputs proteome maps and predicted interactomes in respective folders and creates folders in working directory if not found.

Uses PSIQUIC service to retrieve IntAct data - see https://github.com/micommunity/psicquic

**CHANGES COMPLETE**:

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

Added unique proteins to predicted interactome output...but the counts seem far too small. 

**IN PROGRESS**:

Need to trace back interactome predictions to specific proteins, add to predicted interactome, and get counts.
This work is in progress but output needs to be verified.
For a group of interactome predictions, get counts for the following:
 (For Fig 4A)
 1. Unique proteins in proteome
 2. Proteins not in PPI
 3. Proteins with experimental PPI and direct predictions (either the original PPI, or PPI from shared OGs in the same species)
 4. Proteins with predicted PPI (interactions originally seen in different species)
 
 (For Fig 4B)
 5. Unique OGs in proteome
 6. OGs without interactions
 7. OGs with experimental interactions (see item 3)
 8. OGs with predicted interactions (see item 4)
 
 (For Fig 4C)
 9. Unique OG interactions in the predicted network
 10. Unique experimental OG interactions 

 Also want to know total count of PPI in predicted interactome vs. count of unique proteins, per species

Do interactome prediction for a given proteome. Download proteome on request.

Use protein and species count from eggNOG (it's in the annotation file).

Perform ANOVA between different FuncCats to see consensus interaction patterns.

 Or at least get interaction frequencies by FuncCat (as in filtering goal above)
 
Assign methods to interactions (more general than original data, so we can detect spoke expansion)

 Spoke expansion could actually be filtered out in the input set (use complex:"-") but would rather keep it.


