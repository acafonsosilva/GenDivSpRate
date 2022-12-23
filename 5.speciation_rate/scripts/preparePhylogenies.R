library(ape)
library(tidyverse)
#### Prepare tree data from Nathan Upham to run ClaDS in julia
folder_path <- "/data/biodiv/asilva/5.speciation_rate/"

#### prepare species sampling fraction table
allData <- read.delim(paste0(folder_path,'inputs/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt'))

## not counting extinct species to get sampling fraction
noExtinct <- allData %>% 
  filter(extinct. == 0)

## Sampling fraction obtained per patched clade comparing DNAonly dataset vs complete dataset
final_clades <- data.frame(table(noExtinct[noExtinct$samp %in% 'sampled','PC']), table(noExtinct$PC))[-1]
final_clades$FR <- round(as.numeric(final_clades$Freq)/as.numeric(final_clades$Freq.1),3)
colnames(final_clades) <- c('sampled_notExtinct_PC', 'PC','complete_notExtinct_PC','sampling_fraction_PC')

species_SF <- allData %>% 
  filter(extinct. == 0, samp %in% "sampled") %>% 
  left_join(., final_clades, by = 'PC') %>% 
  mutate(species = word(tiplabel,1,2, sep = "_")) %>%
  select(PC, species, sampling_fraction_PC)

## prepare MCC tree
Upham <- read.tree(paste0(folder_path,'inputs/MamPhy_BDvr_DNAonly_topoFree_NDexp_4098sp_MCC_target.tree')) #includes extinct species that should be excluded
Upham <- ladderize(reorder.phylo(Upham))

Upham$tip.label <- word(Upham$tip.label,1,2, sep = "_")

treeMCC <- drop.tip(Upham, as.character(Upham$tip.label[which(!Upham$tip.label %in% 
                                                                  species_SF$species)]))

## prepare 100 posterior trees
tree100 <- read.nexus(paste0(folder_path,'inputs/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_sample100_nexus.trees'))

TreeSet <- list()
for (i in 1:length(tree100)){
  t <- ladderize(reorder.phylo(tree100[[i]]))
  t$tip.label <- word(t$tip.label,1,2, sep = "_")
  
  TreeSet[[i]] <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% species_SF$species)]))
}

names(TreeSet) <- paste0("tree",seq(1,100))
TreeSet[['treeMCC']] <- treeMCC

save(TreeSet, species_SF, file = paste0(folder_path,'outputs/upham_4064sp_FR_MCCposterior100.rdata'))

