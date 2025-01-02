library(ape)
library(tidyverse)

folder_path <- '/data/biodiv/asilva/'
wd <- paste0(folder_path, 'mutationRate/output/clades/')
dir.create(wd)
  
load(paste0(folder_path,'speciation_rate/output/upham_4064sp_FR_MCCposterior100.rdata')) #TreeSet

ctl <- readLines(paste0(folder_path, 'mutationRate/input/baseml.ctl'))

#### Prepare sequence for each clade
dna <- paste0(folder_path, 'mutationRate/input/FIN3_pruned1_zPrunedMisID_zNumOK_CYTB.fasta')

clades <- read.table(paste0(folder_path,'speciation_rate/input/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt'), header = T) %>% 
  drop_na(PC) %>% 
  mutate(species = word(tiplabel, 1,2, sep = "_"),
         clades = sub("^PC\\d+_",  "", PC)) %>%
  filter(species %in% TreeSet[['treeMCC']]$tip.label)

seqDNA <- seqinr::read.fasta(dna, seqtype = c("DNA"), 
                             seqonly = FALSE, forceDNAtolower = FALSE, set.attributes = FALSE)

seqNames <- data.frame(tips = names(seqDNA)) %>% 
  separate(tips, into= c("species", "family","order","ID1","ID2","ID3","seqSize"), sep="__") 

names(seqDNA) <- seqNames$species

### Some issue related with Carnivora and Muridae that had to be splitted in nodes (decided visualy)
extra <- data.frame(clade = c('Carnivora1', 'Carnivora2', 'Muridae1','Muridae2'), 
                    Node =  c(7365,7213,4976,4644))

for(i in 1:dim(extra)[1]){
 ext <- extract.clade(TreeSet[['treeMCC']], extra[i,'Node'])
 clades[clades$species %in% ext$tip.label,'clades'] <- as.character(extra[i,'clade'])
}

### Make clades folders and split DNA per clade with at least 4 species in phylogeny to estimate mutation rate
for (clade in levels(as.factor(clades$clades))){
  d <- paste0(wd,clade)
  setDNA <- seqDNA[names(seqDNA) %in% clades[clades$clades %in% clade,'species']]
  if(length(setDNA) > 3){  
    dir.create(d)
    write.dna(setDNA, paste0(d,"/",clade,".phy"), format = 'sequential', 
              nbcol = -1, colsep = "")
    ## Extract 3rd codon position for each alignment
    #fix(write.dna) to put two spaces between sp name and sequence
    # fmt <- paste("%-", max.nc + 2, "s", sep = "")  #replace max.nc + 1  by max.nc + 2
    
    setDNA <- read.dna(paste0(d,"/",clade,".phy"), as.character = T, as.matrix = T)
    setDNA_3rd <- setDNA[,seq(from=3, to=dim(setDNA)[2], by=3)]
    write.dna(setDNA_3rd, file = paste0(d,"/",clade,"_3rd.phy"), colsep ="",nbcol=-1)
  }
}

#### Prepare data for each posterior tree
### since for some species there is no cytb information from Nathan's alignment then each tree needs to have the same labels as the DNA alignment
mutClades <- list.dirs(wd, recursive = FALSE, full.names = FALSE)
for ( i in mutClades){
  for (j in names(TreeSet)){
    dir.create(paste0(wd,i,'/',j))
    dna <- read.dna(paste0(wd,i,'/',i,'_3rd.phy'), as.character = TRUE)
    file.copy(paste0(wd,i,'/',i,'_3rd.phy'), paste0(wd,i,'/',j))
    
    tips <- row.names(dna)
    
    tree <- TreeSet[[j]]
    subTree <- drop.tip(tree, tree$tip.label[which(!tree$tip.label %in% tips)])
    write.tree(subTree, paste0(wd,i,'/',j,'/',i,".tree"))

    ctl1 <- gsub("XXX", i,ctl)
    cat(ctl1, file = paste0(wd,i,'/',j,'/baseml.ctl'), sep="\n")
  }
}

write.table(data.frame(mutClades),paste0(folder_path, 'mutationRate/output/clades.txt'), row.names = F, col.names = F, quote = F, sep = '\t')

# for ( i in mutClades){
#   for (j in names(TreeSet)){
# file.copy(paste0('/Volumes/data/mutationRate/clades/',i,"/",j,"/",i, "_res.txt"), paste0('/Users/acas/Dropbox/Post-docs/Morlon_Lab/manuscript_projects_info/mammals_genDiv_SpRate/manuscript_scripts_data/mutationRate/output/clades/',i,"/",j,"/",i, "_res.txt"))
#   }
# }
