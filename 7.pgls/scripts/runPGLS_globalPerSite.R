library(ape)
library(tidyverse)
library(caper)
select = dplyr::select

args <- commandArgs()
print(args)
repi <- as.numeric(args[4])

#### set paths
folder_path <- '/data/biodiv/asilva/'

# discard species without genetic diversity for each of the frames, species which more than 10 individuals, and species not in frame 

checkFrame <- read.delim(paste0(folder_path,
                                '4.genetic_diversity/outputs/speciesFrame.txt')) 

gen.div <- read.delim(paste0(folder_path,
                                   '4.genetic_diversity/outputs/GenDiv_subsample5Ind.txt')) %>%
  filter(EstPi > 0, !outPi %in% 'out', !outTheta %in% 'out',
         species %in% checkFrame[checkFrame$check %in% 'inFrame1','species'])

clades <- read.delim(paste0(folder_path,
                            '5.speciation_rate/inputs/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt')) %>%
  drop_na(PC) %>%
  mutate(species = word(tiplabel, 1,2, sep = "_"),
         clades = sub("^PC\\d+_",  "", PC)) %>%
  dplyr::select(species, clades)

spRate <- read.delim(paste0(folder_path,'5.speciation_rate/outputs/MCCposterior100_tipRate.txt')) %>%
  rename(set = treeN ) %>%
  left_join(., clades, by = 'species')

GendivSpRate10 <- spRate %>%
  arrange(clades) %>%
  inner_join(.,select(gen.div, species, Nind, EstPi,EstPi_1:EstPi_3), by = 'species') %>%
  filter(Nind > 10, EstPi > 0, EstPi_1 > 0, EstPi_2 > 0, EstPi_3 > 0) %>%
  droplevels()

saveRDS(GendivSpRate10, paste0(folder_path,'7.pgls/extra/GendivSpRate10.rds'))



repiTree <- unique(GendivSpRate10$set)[repi]

DataSubset <- GendivSpRate10 %>% 
  filter(set %in% repiTree) %>% 
  drop_na() 

##load phylogenies and speciation rates from 100 posterior trees + MCC tree
load(paste0(folder_path,'5.speciation_rate/outputs/upham_4064sp_FR_MCCposterior100.rdata')) #loads TreeSet
tr <- TreeSet[[repiTree]]

TreeSubset <- drop.tip(tr, 
                       as.character(tr$tip.label[which(!tr$tip.label %in% DataSubset$species)]))
TreeSubset

datacomp <- comparative.data(data = DataSubset, 
                                          phy = TreeSubset, 
                                          names.col = "species", 
                                          vcv = TRUE,
                                          na.omit = TRUE, warn.dropped = TRUE)

tryCatch({
  modEstPi <- NULL
  modEstPi <- pgls(log(EstPi) ~ log(tipRate), data = datacomp, lambda = 'ML')  ## to compare the with same size data set
  
  modEstPi1 <- NULL
  modEstPi1 <- pgls(log(EstPi_1) ~ log(tipRate), data = datacomp, lambda = 'ML')
  
  modEstPi2 <- NULL
  modEstPi2 <- pgls(log(EstPi_2) ~ log(tipRate), data = datacomp, lambda = 'ML')
  
  modEstPi3 <- NULL
  modEstPi3 <- pgls(log(EstPi_3) ~ log(tipRate), data = datacomp, lambda = 'ML')
  
}, error=function(e){cat("Error in",conditionMessage(e), "\n")})

pgls <- list(EstPi = modEstPi,
             EstPi1 = modEstPi1,
             EstPi2 = modEstPi2,
             EstPi3 = modEstPi3)

saveRDS(pgls, paste0(wd, '7.pgls/outpus/global_gendivSpRatePerSite_results_',repiTree,'.rds')) ### kept in the cluster
