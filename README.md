# GenDivSpRate

This is a test


1. Folder_per_family.R:
* Extract GenBank sequence for a given family to get metadata per family using ncbi taxonomy database and prepare folders to parallelise SuperCrunch per family 
* Output:
    *  <family>.fasta,
    * Taxon_families.txt, - for SuperCrunch pipeline
    * Taxon_families_Nseq.txt, - number of sequences per family
    * genBank_metadata_complete.rds - to make the matching easier - genBank organism name
    * cleanSeqLabels.txt - taxa list for all retrieved sequences but excluding names with problems (abbreviated genus, sp., other characters, and hybrid species names).

2. nameMatching - makeSynonymsCleanList.R
* build up synonyms list based on ITIS database
* input: df_Syns_mamm_update2016.txt - using Meyer et al. (2015) synonyms list
* output: synonymsClean.rds - new synonyms list
* Notes:
        - all valid itis species (in Accepted) are in synonyms column with valid flag in nameUsage
        - SynSpp column are extra synonyms matching synonyms column
        - nameUsage refers to the synonyms column, meaning that potential lumped or splited (problems) species in synSpp are not flagged
        - tsn columns have ids for itis but iucn and otl species names are flagged with just 'iucn' and 'otl'

3. nameMatching - addUphamSynonymsToSynList.R
* Match tree species names to this new synonyms list and search for non matched species names including iucn and tol database
* input: synonymsClean.rds and upham_NCBI_matchup.txt (metadata of Upham phylogenies matching IUCN with NCBI)
* output: synonymsClean_uphamUpdated.rds

4. nameMatching - addNCBISynonymsToSynList.R
* match sequence species names to same synonyms list and search for non matched species names including iucn and tol database
* input: synonymsClean_uphamUpdated.rds and genBank_metadata_complete.rds
* output: synonymsClean_uphamUpdatedNCBIupdated.rds

5. nameMatching - makeTreeLabels_matchSyn.R
* match every tree species name to a synonym ID to later better match directly to sequence species names. Only matches species names in DNAonly tree
* input: synonymsClean_uphamUpdatedNCBIupdated.rds and MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt (Upham phylogenies metadata)
* output: treelabels_synIDs.rds

6. nameMatching - makeSeqLabels_matchSynUpham.R
* match every sequence species names to a tree species names. Mostly through ID in synonyms list but also if sequence species name in tree but not in list given certain conditions  
* input: genBank_metadata_complete.rds, synonymsClean_uphamUpdatedNCBIupdated.rds, treelabels_synIDs.rds
* output: genBank_metadata_complete_matchedUpham.txt, seqlabels_matchedUpham.txt
* Potential problems with different conditions for each
            - valid - maybe ok
            - subspecies - maybe ok
            - synonymised - maybe ok
            - lumped - maybe ok
            - invalid - to check
            - invalidspp - to check
            - iucn - to check
            - iucn_valid - to check
            - extinct - exclude
            - iucn_problSyn - likely to exclude
            - otl_problSyn - likely to exclude
            - problSyn - likely to exclude
            - neverFound - keep if sp in tree and if spp12 is but spp13 not or spp13 is but spp12 not
            
7. superCrunch_family - superCrunch.sh
* input: Taxon_families.txt, <family>.fasta, cleanSeqLabels.txt, genBank_metadata_complete_matchedUpham.txt
    1. Taxa_Assessment.py
    2. Reference_Blast_Extract.py
    3. Adjust_Direction.py
    4. Coding_Translation_Tests.py
    5. Align.py
    6. Fasta_file_per_species.R <- not a SuperCrunch script
    7. Trim_Alignments_Trimal.py
* output: Output_alignments

8. genetic_diversity - geneticDiversity.R
* Estimate pi and theta for species with at least five individuals using Ferretti et al. 2012 derived formulas to handle missing data, and for species with more than five individuals subsample 1000 times to five individuals.
* input: Output_alignments
* output: GenDiv_subsample5Ind.txt, GenDiv_subsample5Ind1000rep.rds
* Notes: subsampling takes time but this script is to do it without parallelising calculation.

9. speciation_rate - preparePhylogenies.R
* Prepare MCC and 100 posterior trees from Nathan Upham’s DNA only phylogenies, as well per clade sampling fraction
* input: MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, MamPhy_BDvr_DNAonly_topoFree_NDexp_4098sp_MCC_target.tree, MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_sample100_nexus.trees
* output: upham_4064sp_FR_MCCposterior100.rdata

10. speciation_rate - cluster_subClaDS2_Julia.sh > ClaDS2_upham_4064sp_MCCposterior100.jl
* To run Odile Maliet ClaDS2 model implemented in Julia with data augmentation. Needs a cluster submission file to parallelise 100 posterior trees and MCC tree analyses - 
* input: upham_4064sp_FR_MCCposterior100.rdata
* output: upham_tree<tree1:MCC>.Rdata
* Notes: Each tree took between 1 day to a bit over a week to complete.

11. speciation_rate - extractMAPS.R
* Extracts the Maximum A Posteriori (MAPS) for each of the parameters
* input:  upham_tree<tree1:MCC>.Rdata
* output: MCCposterior100_pars.txt (hyperparameters), MCCposterior100_tipRate.txt (branch and tip species-specific rates, MCCposterior100_treeMAPS.rds (MAPS from all the analyses)

12. pgls - preparePGLS.R
* input: GenDiv_subsample5Ind.txt, upham_4064sp_FR_MCCposterior100.rdata, MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, MCCposterior100_tipRate.txt
* output: global/datacomp_global_gendivSpRate.rds, clade/datacomp_clade_gendivSpRate.rds

13. pgls - subpgls_global.sh > runPGLS_global.R
* 101 global sets - 100 posterior trees + MCC tree
* analyses:
    * datacomp_global_gendivSpRate:
        * estimated pi vs speciation rate - running
        * mean subsampled pi vs speciation rate - running
        * estimated theta vs speciation rate - running
        * mean subsampled theta vs speciation rate - running
* input: global/datacomp_global_gendivSpRate.rds
* output: global/output/global_gendivSpRate_results_<tree>.rds

14. pgls - subpgls_clade.sh > runPGLS_clade.R
* Analysis for 14 clades that have data for at least 20 species (each with 101 sets)
* analyses:
    * datacomp_clades_gendivSpRate:
        * estimated pi vs speciation rate - running
        * mean subsampled pi vs speciation rate - running
        * estimated theta vs speciation rate - running 
        * mean subsampled  theta vs speciation rate - running
* input: clade/datacomp_clade_gendivSpRate.rds
* output: clade/output/clade_gendivSpRate_results_<clade>.rds

15. pgls - extractPGLS.R

16. bmlm - prepare bmlm
17. bmlm - sub > bmlm
18. bmlm - extract bmlm 

19. pgls - preparePGLS_traits.R
* input: GenDiv_subsample5Ind.txt, upham_4064sp_FR_MCCposterior100.rdata, MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, MCCposterior100_tipRate.txt, TRAITS
* output: datacomp_global_gendivSpRateTraits.rds

20. pgls - subpgls_global.sh > runPGLS_global_traits.R
* 101 global sets - 100 posterior trees + MCC tree
* analyses:
    * datacomp_global_gendivSpRateTraits:
        * estimated pi vs speciation rate + traits,
        * estimated pi vs traits, 
        * speciation rate vs traits
* input: datacomp_global_gendivSpRateTraits.rds
* output: global_gendivSpRateTraits_results_<tree>.rds

21. pgls - extractPGLS_traits.R

22. bmlm - prepare bmlm with traits
23. bmlm - sub > bmlm with traits
24. bmlm - extract bmlm with traits

25. mutation rate - prepareMutRate.R
* Input: upham_4064sp_FR_MCCposterior100.rdata, baseml.ctl (from Leandro), MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, FIN3_pruned1_zPrunedMisID_zNumOK_CYTB.fasta (from Nathan)
* Output: clades/<mutclade>/<100trees+MCC>/<mutclade>_3rd.phy & <mutclade>.tree & baseml.ctl (edited for each clade), clades.txt (mutclade - adjusted for a few clades)

26. mutation rate - subBaseml.sh > baseml.sh
* Runs PAML - baseml program for each tree for each mutclade with settings defined in baseml.ctl
* For clades with at least 4 species
* Easier to setup batches of 5 clades at a time and ran MCC trees in parallel independently

27. mutation rate - extractMutRate.R
* input: upham_4064sp_FR_MCCposterior100.rdata,  clades/<mutclade>/<100trees+MCC>/<mutclade>_res.txt
* output: mutationRate_cytb3rdcodonPAML.rds

