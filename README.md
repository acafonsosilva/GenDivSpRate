# Negative global-scale association between genetic diversity and speciation rates in mammals

This is a little overview of the scripts needed to run the analyses of the manuscript submitted entitled "Negative global-scale association between genetic diversity and speciation rates in mammals". Currently only scripts are stored on github, all inputs and outputs can be found at https://figshare.com/s/8917e270827adce72b89.

### Figures_v3.Rmd

Rmarkdown file with code used to make all the figures, including some extra models mentioned in the manuscript. It also includes the code used to export the Source Data files used for each figure, versions of the used packages and exports the R working environment. Output as HTML file Figures_v3.html that can be open in a browser.

## Description of the scripts used in the analyses

### 1. Folder_per_family.R:
* Extract GenBank sequence for a given family to get metadata per family using NCBI taxonomy database and prepare folders to parallelize SuperCrunch per family 
* Output:
    *  <family>.fasta,
    * Taxon_families.txt, - for SuperCrunch pipeline
    * Taxon_families_Nseq.txt, - number of sequences per family
    * genBank_metadata_complete.rds - to make the matching easier - GenBank organism name
    * cleanSeqLabels.txt - taxa list for all retrieved sequences but excluding names with problems (abbreviated genus, sp., other characters, and hybrid species names).

### 2.nameMatching 
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
* output: genBank_metadata_complete_matchedUpham.txt
            
### 3.superCrunch_family - superCrunch.sh
    1. Taxa_Assessment.py
    2. Reference_Blast_Extract.py <- customized helper script extractReferencefamily.R
    3. Adjust_Direction.py
    4. Coding_Translation_Tests.py
    5. Align.py
    6. Fasta_file_per_species.R <- not a SuperCrunch script
    7. Trim_Alignments_Trimal.py
    Further manual checking step with script fix_notInFrame_getAccession.R
* input: cleanSeqLabels.txt, genBank_metadata_complete_matchedUpham.txt
* output: Output_alignments (too many intermediate files from superCrunch so only the final alignments are shared as Cytb_alignments.zip); GenDiv_accessionIDsPerSpecies.txt (file produced later and added code to fix_notInFrame_getAccession.R script; it's the genBank_metadata_complete_matchedUpham.txt file filtered for sequences in the final output_alignments)

### 4.genetic_diversity

4.1 - resamplingSynNonSyn_faster.R
* Estimate pi and theta for species with at least five individuals using Ferretti et al. 2012 derived formulas to handle missing data, and for species with more than five individuals subsample 1000 times to five individuals.
* input: 3.superCrunch_family/Output_alignments
* output: GenDiv_SynNonSyn_resampled4ind.txt, GenDivSyn_resample4ind_output.rds
* Notes: subsampling takes time so the script shows the code used to parallelise calculations done in a cluster.

4.2 - geneticDiversity_GeographicLocationChecks.R
* Takes the coordinates used in Theodoridis et al. 2021 to get genetic diversity data for a smaller subset of species with available geographic information.
* input: 3.superCrunch_family/Output_alignments; otherData/Theodoridis2021/cytb_coordinates.csv'
* output: GenDiv_geo.txt

### 5.speciation_Rate 

5.1 - preparePhylogenies.R
* Prepare MCC and 100 posterior trees from Nathan Upham’s DNA only phylogenies, as well per clade sampling fraction
* input: MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, MamPhy_BDvr_DNAonly_topoFree_NDexp_4098sp_MCC_target.tree, MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_sample100_nexus.trees
* output: upham_4064sp_FR_MCCposterior100.rdata

5.2 - cluster_subClaDS2_Julia.sh > ClaDS2_upham_4064sp_MCCposterior100.jl; extractMAPS.R
* To run Odile Maliet's ClaDS2 model implemented in Julia with data augmentation before it was available as a package the scrips in ClaDS_Julia_scripts.zip were used, so the names of the functions are a bit different from what is currently implemented in PANDA.jl (https://hmorlon.github.io/PANDA.jl/stable/#ClaDS). A cluster submission file to parallelise 100 posterior trees and MCC tree analyses was needed. To extract the Maximum A Posteriori (MAPS) for each of the parameters, the tip rates and the hyperparameters from each of the analysed trees.
* input: upham_4064sp_FR_MCCposterior100.rdata
* output: upham_tree<tree1:MCC>.Rdata (data was a bit too heavy so stayed in the cluster and just kept locally the output files); MCCposterior100_pars.txt (hyperparameters), MCCposterior100_tipRate.txt (branch and tip species-specific rates, MCCposterior100_treeMAPS.rds (MAPS from all the analyses)

5.3 - match_edges.R
* Customized function to plot Figure 1 and run in the Figures_v3.Rmd rmarkdown file.

### 6.mutationRate 

6.1 - prepareMutRate.R
* Input: 5.speciation_Rate/outputs/upham_4064sp_FR_MCCposterior100.rdata, baseml.ctl, 5.speciation_Rate/inputs/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, FIN3_pruned1_zPrunedMisID_zNumOK_CYTB.fasta
* Output: For each mutclade and tree, a folder with: /clades/[mutclade]/[MCC_100trees]//[mutclade]3rd.phy; [mutclade].tree; baseml.ctl. clades.txt (file that shows the mutclades because some clades had to be adjusted for analysis to converge)

6.2 - mutation rate - subBaseml.sh > baseml.sh
* Runs PAML - baseml program for each tree for each mutclade with settings defined in baseml.ctl
* For clades with at least 4 species
* Easier to setup batches of 5 clades at a time and ran MCC trees in parallel independently in the cluster

6.3 - mutation rate - extractMutRate.R
* input: upham_4064sp_FR_MCCposterior100.rdata,  clades/<mutclade>/<100trees+MCC>/<mutclade>res.txt
* output: mutationRate_cytb3rdcodonPAML.rds

### 7.pgls 

7.1 - subpgls_All.sh > runPGLS_v2_REML.R > extractPGLS_v2.R
* 100 posterior trees + MCC tree
* analyses:
        * global analyses (estimated pi vs speciation rate, mean subsampled pi vs speciation rate, estimated theta vs speciation rate, mean subsampled theta vs speciation rate). Analyses done with measures estimated from all sites, synonymous and nonsynonymous sites.
        * per clade analyses (estimated synonymous pi vs speciation rate, mean subsampled synonymous pi vs speciation rate, estimated synonymous theta vs speciation rate, mean synonymous subsampled theta vs speciation rate)
        * trait analyses: estimated syn pi vs speciation rate + traits (body size, latitudinal midpoint, mean temperature across range area, litter size, generation length), estimated syn pi vs traits, speciation rate vs traits
        * mutation rate analyses: mutation rate vs speciation rate, Ne (from mutation rate and estimated syn pi) vs speciation rate. Done with both mutRate ((expNsub * GenerationLength_y) / timeyear) and mutRate_y (expNsub/timeyear), expNsub are the substitution rates estimated with PAML used as proxy for mutation rates.
        * an Extra analysis was run when extracting the pgls results to get R2 according to Ives et al. 2018.
* input: GenDiv_SynNonSyn_resampled4ind.txt, upham_4064sp_FR_MCCposterior100.rdata, MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, MCCposterior100_tipRate.txt, otherData/traits/matchedTraits_v3.txt, mutationRate_cytb3rdcodonPAML.rds
* output: gendivSpRate_PGLSresultsAll_REML.rds, clade_gendivSpRate_PGLSresults.rds, gendivSpRate_PGLS_R2EstPiSyn.rds

7.2 - runPGLS_extraTraits.R
* Added later to reply to reviewer comments. The script was run in the cluster and contains how it was extracted. In the end analyses with range area were not included in the manuscript.
* input: same as previous
* output: extraTraits_modelOutputs.rds

* Extra folder has several files created in Figures_v3.Rmd to complement analyses and help plotting.


### 8.bmlm

8.1 - subBMLM.sh > runBMLM_global_v2.R | runBMLM_global_traits_v2.R | runBMLM_global_mutRate_v2.R >extractBMLM_globalAll.R
* 100 posterior trees + MCC tree (each analysis accounts for clade data as random factor)
* analyses:
        * global analyses (estimated pi vs speciation rate, mean subsampled pi vs speciation rate, estimated theta vs speciation rate, mean subsampled theta vs speciation rate). Analyses were only done with synonymous sites measures.
        * trait analyses: estimated syn pi vs speciation rate + traits (body size, latitudinal midpoint, mean temperature across range area, litter size, generation length), estimated syn pi vs traits, speciation rate vs traits
        * mutation rate analyses: mutation rate vs speciation rate, Ne (from mutation rate and estimated syn pi) vs speciation rate. Done only with mutRate_y|Ne_y.
* input: GenDiv_SynNonSyn_resampled4ind.txt, upham_4064sp_FR_MCCposterior100.rdata, MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt, MCCposterior100_tipRate.txt, otherData/traits/matchedTraits_v3.txt, mutationRate_cytb3rdcodonPAML.rds
* output: global_allPosterior.rds, global_traits_allPosterior.rds, global_mutRate_allPosterior_year.rds




