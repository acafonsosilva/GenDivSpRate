# GenDivSpRate

This is a test


1. Folder_per_family.R:
* Extract GenBank sequence for a given family to get metadata per family using NCBI taxonomy database and prepare folders to parallelize SuperCrunch per family 
* Output:
    *  <family>.fasta,
    * Taxon_families.txt, - for SuperCrunch pipeline
    * Taxon_families_Nseq.txt, - number of sequences per family
    * genBank_metadata_complete.rds - to make the matching easier - GenBank organism name
    * cleanSeqLabels.txt - taxa list for all retrieved sequences but excluding names with problems (abbreviated genus, sp., other characters, and hybrid species names).

2. nameMatching 
* Scripts that improve how many sequences in GenBank can be matched into Upham et al. 2019 tree

2.1 - makeSynonymsCleanList.R
* Build up a synonyms list based on ITIS database
* input: df_Syns_mamm_update2016.txt - using Meyer et al. (2015) synonyms list to start
* output: synonymsClean.rds - new synonyms list
* Notes:
        - all valid itis species (in Accepted) are in the synonyms column with a valid flag in nameUsage
        - SynSpp column is an extra synonym matching column
        - nameUsage refers to the synonyms column, meaning that potential lumped or splited (problems) species in synSpp are not flagged
        - tsn columns have ids for itis but iucn and otl species names are flagged with just 'iucn' and 'otl'

2.2 - addUphamSynonymsToSynList.R
* Match tree species names to this new synonyms list and search for non-matched species names including iucn and tol database
* input: synonymsClean.rds and upham_NCBI_matchup.txt (metadata of Upham phylogenies matching IUCN with NCBI)
* output: synonymsClean_uphamUpdated.rds

2.3 - addNCBISynonymsToSynList.R
* Match sequence species names to the same synonyms list and search for non-matched species names including iucn and tol database
* input: synonymsClean_uphamUpdated.rds and genBank_metadata_complete.rds
* output: synonymsClean_uphamUpdatedNCBIupdated.rds

2.4 - makeTreeLabels_matchSyn.R
* Match every tree species name to a synonym ID to later better match directly to sequence species names. Only matches species names in DNAonly tree
* input: synonymsClean_uphamUpdatedNCBIupdated.rds and MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt (Upham phylogenies metadata)
* output: treelabels_synIDs.rds

2.5 - makeSeqLabels_matchSynUpham.R
* Match all sequence species names to a tree species names. Mostly through ID in synonyms list but also if sequence species name is in the tree but not in the list given certain conditions  
* input: genBank_metadata_complete.rds, synonymsClean_uphamUpdatedNCBIupdated.rds, treelabels_synIDs.rds
* output: genBank_metadata_complete_matchedUpham.txt, seqlabels_matchedUpham.txt
            
3. superCrunch_family - superCrunch.sh
* input: Taxon_families.txt (need to get the file from the cluster), <family>.fasta (need to get the files from the cluster), cleanSeqLabels.txt, genBank_metadata_complete_matchedUpham.txt
    1. Taxa_Assessment.py
    2. Reference_Blast_Extract.py
    3. Adjust_Direction.py
    4. Coding_Translation_Tests.py
    5. Align.py
    6. Fasta_file_per_species.R <- not a SuperCrunch script
    7. Trim_Alignments_Trimal.py
* output: Output_alignments (too many intermediate files from superCrunch so only the final alignments are shared)

4. genetic_diversity

4.1 - geneticDiversity.R
* Estimate pi and theta for species with at least five individuals using Ferretti et al. 2012 derived formulas to handle missing data, and for species with more than five individuals subsample 1000 times to five individuals.
* input: 3.superCrunch_family/Output_alignments
* output: GenDiv_subsample5Ind.txt, GenDiv_subsample5Ind1000rep.rds
* Notes: subsampling takes time but this script is without the parallelising calculations done in a cluster.

4.2 - geneticDiversityInFrame.R
* Estimate pi and theta for each of the three codon site positions for species with at least five individuals 
* input: 3.superCrunch_family/Output_alignments
* output: GenDiv_frame123.txt, speciesFrame.txt

5. speciation_rate 

5.1 - preparePhylogenies.R
* Prepare MCC and 100 posterior trees from Nathan Upham’s DNA only phylogenies, as well per clade sampling fraction
* input: MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, MamPhy_BDvr_DNAonly_topoFree_NDexp_4098sp_MCC_target.tree, MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_sample100_nexus.trees
* output: upham_4064sp_FR_MCCposterior100.rdata

5.2 - cluster_subClaDS2_Julia.sh > ClaDS2_upham_4064sp_MCCposterior100.jl
* To run Odile Maliet's ClaDS2 model implemented in Julia with data augmentation before it was available as a package the scrips in ClaDS_Julia_scripts.zip were used, so the names of the functions are a bit different from what is currently implemented in PANDA.jl (https://hmorlon.github.io/PANDA.jl/stable/#ClaDS). A cluster submission file to parallelise 100 posterior trees and MCC tree analyses was needed.
* input: upham_4064sp_FR_MCCposterior100.rdata
* output: upham_tree<tree1:MCC>.Rdata

5.3 - extractMAPS.R
* To extract the Maximum A Posteriori (MAPS) for each of the parameters, the tip rates and the hyperparameters from each of the analysed trees.
* input:  upham_tree<tree1:MCC>.Rdata (data was a bit too heavy so stayed in the cluster and just kept locally the output files)
* output: MCCposterior100_pars.txt (hyperparameters), MCCposterior100_tipRate.txt (branch and tip species-specific rates, MCCposterior100_treeMAPS.rds (MAPS from all the analyses)

* Can't find the code used to make the file Upham_FBD_clades_nodes.txt but it's mostly a summary of per clade sampling fraction and the node position for each clade for plotting in figure 1.

6.mutationRate (to be improved)

6.1 - prepareMutRate.R
* Input: upham_4064sp_FR_MCCposterior100.rdata, baseml.ctl, MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, FIN3_pruned1_zPrunedMisID_zNumOK_CYTB.fasta
* Output: clades/<mutclade>/<100trees+MCC>/<mutclade>_3rd.phy & <mutclade>.tree & baseml.ctl (edited for each clade), clades.txt (mutclade - adjusted for a few clades)

6.2 - mutation rate - subBaseml.sh > baseml.sh
* Runs PAML - baseml program for each tree for each mutclade with settings defined in baseml.ctl
* For clades with at least 4 species
* Easier to setup batches of 5 clades at a time and ran MCC trees in parallel independently

6.3 - mutation rate - extractMutRate.R
* input: upham_4064sp_FR_MCCposterior100.rdata,  clades/<mutclade>/<100trees+MCC>/<mutclade>_res.txt
* output: mutationRate_cytb3rdcodonPAML.rds

7. pgls 

7.1 - subpgls_All.sh > runPGLS.R
* 100 posterior trees + MCC tree
* analyses:
        * global analyses (estimated pi vs speciation rate, mean subsampled pi vs speciation rate, estimated theta vs speciation rate, mean subsampled theta vs speciation rate)
        * per clade analyses (estimated pi vs speciation rate, mean subsampled pi vs speciation rate, estimated theta vs speciation rate, mean subsampled theta vs speciation rate)
        * trait analyses: estimated pi vs speciation rate + traits (body size, range area, mean temperature across range area, litter size, generation length), estimated pi vs traits, speciation rate vs traits
        * mutation rate analyses: mutation rate vs speciation rate, Ne (from mutation rate and estimate pi) vs speciation rate.       
* input: GenDiv_subsample5Ind.txt, upham_4064sp_FR_MCCposterior100.rdata, MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, MCCposterior100_tipRate.txt, otherData/traits/matchedTraits.txt, mutationRate_cytb3rdcodonPAML.rds
* output stored in the cluster: global_gendivSpRate_results_<tree>.rds, clade_gendivSpRate_results_<tree>.rds,  global_gendivSpRateTraits_results_<tree>.rds, global_gendivSpRateMutRate_results_<tree>.rds

7.2 - subpgls_globalPerSite.R
* Without species that have no genetic diversity for each of the codon sites, just species which more than 10 individuals, and without species that are not in frame 
* analyses: estimated pi vs speciation rate, pi from first codon sites vs speciation rate, pi from second codon sites vs speciation rate, pi from third codon sites vs speciation rate
* input: GenDiv_subsample5Ind.txt, upham_4064sp_FR_MCCposterior100.rdata, MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, MCCposterior100_tipRate.txt, speciesFrame.txt
* output: global_gendivSpRatePerSite_results_<tree>.rds (stored in the cluster), extra/GendivSpRate10.rds


* extra has several files created for preliminary analyses and to help plotting

7.3 - extractPGLS.R

* input stored in the cluster: global_gendivSpRate_results_<tree>.rds, clade_gendivSpRate_results_<tree>.rds,  global_gendivSpRateTraits_results_<tree>.rds, global_gendivSpRateMutRate_results_<tree>.rds, global_gendivSpRatePerSite_results_<tree>.rds
* output: global_gendivSpRate_PGLSresults.rds, clade_gendivSpRate_PGLSresults.rds, global_gendivSpRateTraits_PGLSresults.rds, global_gendivSpRateMutRate_PGLSresults.rds, global_gendivSpRatePerSite_PGLSresults.rds


### still to finish fixing

8. bmlm

subBMLM.shrunBMLM_global_mutRate.RextractBMLM_globalAll.R

otherData

Figures.Rmd





