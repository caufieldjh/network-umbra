# network-umbra.py

Predicts interactions in a protein interaction network based off a meta-interactome network.
Intended for use with bacterial proteins.
Uses eggNOG v.4.1 for orthology assignments.
Can also retrieve proteomes from Uniprot, assign their proteins to orthologs, and predict interactions among those orthologs.

Uses Entrez services.
Change the Entrez.email value to your email address before use.

**REQUIRES**: 

**Biopython 1.65** or more recent. Try installing with **pip**:

    pip install numpy

    pip install biopython

or see http://biopython.org/DIST/docs/install/Installation.html

Also requires **BeautifulSoup 4**. Try:

    pip install beautifulsoup4

or see http://www.crummy.com/software/BeautifulSoup/bs4/doc/

For storage, network-umbra needs ~5 GB of available disk space to accomodate data files and output. More space may be necessary for downloaded proteome files.

**INPUT**: 

Downloads all available protein-protein interactions for bacteria from IntAct.

*NOTE: As of Feb 9, 2016, IntAct taxonomy filtering is not working properly.  Please use manually-downloaded interaction databases for now.*

Alternatively, uses a provided PPI data file in PSI-MI TAB 2.7 format.

REMOVE THE HEADER ROW if it's present!
  
Downloads highest-level (LUCA) and bacteria-specific Uniprot ID to NOG mappings from eggNOG v.4.1.

Uses PSICQUIC service to retrieve IntAct data - see https://github.com/micommunity/psicquic

**OUTPUT**: 

Outputs proteome maps and predicted interactomes in respective folders and creates folders in working directory if not found.

'metainteractome[date].txt'
- A meta-interactome composed of all available bacterial protein-protein interactions.
   Follows PSI-MI Tab27 format, with the addition of two ortholog identifiers per row.
   See format description at https://code.google.com/p/psimi/wiki/PsimiTab27Format
   
'meta_statistics[date].txt'
- Contains statistics relevant to the produced meta-interactome.

'taxid_context[date].txt'
- Contains NCBI taxonomy IDs, names, parent IDs, and domains for all input interactions.
   All domains should be Bacteria.
   Used as a reference if meta-interactome and consensus meta-interactome not
   built during the same session.
   
'consensus[date].txt'
- A consensus meta-interactome composed of all available bacterial protein-protein interactions.
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
- Contains statistics relevant to the produced consensus meta-interactome.
   
'subgraph_expansion_[FuncCat]_[date].txt'
- A set of subgraphs of the consensus meta-interactome, filtered by conservation of interactions and function of interactors.Each line is one interaction between a consensus interactor and a unique protein, accompanied by the source taxid of the unique protein.
   
'subgraph_expansion_[FuncCat]_nodes_[date].txt'
- Annotation file for the nodes in the expanded subgraphs.
   

**PLANNED ADDITIONS**:

Need to trace back interactome predictions to specific proteins, add to predicted interactome, and get counts.

Also want to know total count of PPI in predicted interactome vs. count of unique proteins, per species

Use protein and species count from eggNOG (it's in the annotation file).

Perform ANOVA between different FuncCats to see consensus interaction patterns.

Or at least get interaction frequencies by FuncCat (as in filtering goal above)
 
Assign methods to interactions (more general than original data, so we can detect spoke expansion)

Spoke expansion could actually be filtered out in the input set (use complex:"-") but would rather keep it.

Use eggNOG HMMs to extend orthology assignments when mappings are not available.
