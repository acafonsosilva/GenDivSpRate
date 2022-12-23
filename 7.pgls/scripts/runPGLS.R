library(ape)
library(tidyverse)
library(caper)
select = dplyr::select

args <- commandArgs()
print(args)
repi <- as.numeric(args[4])

#### load data
folder_path <- '/data/biodiv/asilva/'

## discard species without genetic diversity and species which either mean subsampled pi or theta are outside the range of 1000 subsamples to 5 individuals
gen.divCIc_sc <- read.delim(paste0(folder_path, '4.genetic_diversity/GenDiv_subsample5Ind.txt')) %>%
  filter(EstPi > 0, !outPi %in% 'out', !outTheta %in% 'out') ### this way doesn't exclude species with 5 individuals


##load phylogenies and speciation rates from 100 posterior trees + MCC tree
load(paste0(folder_path,'5.speciation_rate/outputs/upham_4064sp_FR_MCCposterior100.rdata')) #loads TreeSet

clades <- read.delim(paste0(folder_path, 
                            '5.speciation_rate/inputs/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt')) %>%
  drop_na(PC) %>%
  mutate(species = word(tiplabel, 1,2, sep = "_"),
         clades = sub("^PC\\d+_",  "", PC)) %>%
  dplyr::select(species, clades) 

spRate <- read.delim(paste0(folder_path,'5.speciation_rate/outputs/MCCposterior100_tipRate.txt')) %>% 
  rename(set = treeN ) %>% 
  left_join(., clades, by = 'species') 

gendivSpRate <- dplyr::select(gen.divCIc_sc, species, EstPi, EstTheta) %>%
  left_join(., spRate, by = 'species') 

gendivSpRateSub <- dplyr::select(gen.divCIc_sc, species, subPi_mean, subTheta_mean) %>%
  left_join(., spRate, by = 'species')


repiTree <- levels(gendivSpRate$set)[repi]
t <- TreeSet[[repiTree]]

#### global ####

##Prepare global analyses sets - parallelise per tree
## Prepare Estimated dataset independently of subsampled dataset with species with more than 5 individuals and resampled 1000x to 5
DataSubset <- gendivSpRate %>% 
  filter(set %in% repiTree) %>% 
  droplevels() 
dim(DataSubset)

TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
TreeSubset 

datacompGlobal <- comparative.data(data = DataSubset, 
                                   phy = TreeSubset, 
                                   names.col = "species", 
                                   vcv = TRUE,
                                   na.omit = TRUE, warn.dropped = TRUE)

DataSubsetSub <- gendivSpRateSub %>% 
  filter(set %in% repiTree) %>% 
  droplevels() 

TreeSubsetSub <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubsetSub$species)]))
TreeSubset 

datacompGlobalSub <- comparative.data(data = DataSubsetSub, 
                                   phy = TreeSubsetSub, 
                                   names.col = "species", 
                                   vcv = TRUE,
                                   na.omit = TRUE, warn.dropped = TRUE)

if(!paste0('global_gendivSpRate_results_',repiTree,'.rds') %in% list.files(paste0(wd,'output'))){
  tryCatch({
    #estimated pi vs speciation rate at the tips
    modpglsEstPi <- NULL
    modpglsEstPi <- pgls(log(EstPi) ~ log(tipRate), data = datacompGlobal, lambda = 'ML')
    
    #mean subsampled pi vs speciation rate at the tips
    modpglssubPi <- NULL
    modpglssubPi <- pgls(log(subPi_mean) ~ log(tipRate), data = datacompGlobalSub, lambda = 'ML')
    
    #estimated theta vs speciation rate at the tips
    modpglsEstTheta <- NULL
    modpglsEstTheta <- pgls(log(EstTheta) ~ log(tipRate), data = datacompGlobal, lambda = 'ML')
    
    #mean subsampled theta vs speciation rate at the tips
    modpglssubtheta <- NULL
    modpglssubtheta <- pgls(log(subTheta_mean) ~ log(tipRate), data = datacompGlobalSub, lambda = 'ML')
  }, error=function(e){cat("Error in",conditionMessage(e), "\n")})
  
  pgls <- list(EstPi = modpglsEstPi, subPi = modpglssubPi,
               EstTheta = modpglsEstTheta, subtheta = modpglssubtheta)
  
  saveRDS(pgls, paste0(wd, '7.pgls/outputs/global_gendivSpRate_results_',repiTree,'.rds')) ### kept in the cluster
}

#### clades #### 
##Prepare per clade analyses sets without traits - parallelise per clade

freqClade <- table(gendivSpRate[gendivSpRate$set %in% 'treeMCC','clades'])
datacompClades <- list()
datacompCladesSub <- list()
for (clade in levels(as.factor(gendivSpRate$clades))){
  if(clade %in% names(freqClade[freqClade > 19])){ ##only for clades with at least 20 species
    
    DataSubset <- gendivSpRate %>% 
      filter(set %in% repiTree,  clades %in% clade) %>% 
      droplevels()   
    
    TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
    
    datacompClades[[clade]] <- comparative.data(data = DataSubset, 
                                                phy = TreeSubset, 
                                                names.col = "species", 
                                                vcv = TRUE, 
                                                na.omit = TRUE, warn.dropped = TRUE)
    
    DataSubsetSub <- gendivSpRateSub %>% 
      filter(set %in% repiTree,  clades %in% clade) %>% 
      droplevels() 
    
    TreeSubsetSub <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubsetSub$species)]))
    
    datacompCladesSub <- comparative.data(data = DataSubsetSub, 
                                          phy = TreeSubsetSub, 
                                          names.col = "species", 
                                          vcv = TRUE,
                                          na.omit = TRUE, warn.dropped = TRUE)
  }
}


pgls <- list()
for (clade in names(datacompClades)){
  tryCatch({
    
    #estimated pi vs speciation rate at the tips
    modpglsEstPi <- NULL
    modpglsEstPi <- pgls(log(EstPi) ~ log(tipRate), data = datacompClades[[clade]], lambda = 'ML')
    
    # #mean subsampled pi vs speciation rate at the tips
    modpglssubPi <- NULL
    modpglssubPi <- pgls(log(subPi_mean) ~ log(tipRate), data = datacompCladesSub[[clade]], lambda = 'ML')
    
    #estimated theta vs speciation rate at the tips
    modpglsEstTheta <- NULL
    modpglsEstTheta <- pgls(log(EstTheta) ~ log(tipRate), data = datacompClades[[clade]], lambda = 'ML')
    
    # #mean subsampled theta vs speciation rate at the tips
    modpglssubtheta <- NULL
    modpglssubtheta <- pgls(log(subTheta_mean) ~ log(tipRate), data = datacompCladesSub[[clade]], lambda = 'ML')
    
    pgls[[clade]] <- list(EstPi = modpglsEstPi, SubPi = modpglssubPi,
                      EstTheta = modpglsEstTheta, SubTheta = modpglssubtheta)
    
    print(clade)
  }, error=function(e){cat("Error in",conditionMessage(e), "\n")})
}

saveRDS(pgls, paste0(wd,'7.pgls/outputs/clade_gendivSpRate_results_',repiTree,'.rds'))  ### kept in the cluster

#### Traits ####

traitData <- read.delim(paste0(folder_path,'otherData/traits/matchedTraits.txt'), stringsAsFactors = FALSE) %>%
  mutate(mean_temp = mean_temp + abs(min(mean_temp, na.rm = T)) + 1) 

gendivSpRateTrait <- dplyr::select(gen.divCIc_sc, species, EstPi, EstTheta) %>%
  left_join(., spRate, by = 'species') %>%
  left_join(., traitData, by = 'species')

#### Prepare trait analyses sets only for global analysis - parallelise per tree
DataSubset <- gendivSpRateTrait %>% 
  filter(set %in% repiTree) %>% 
  drop_na() 

TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
TreeSubset 

datacompGlobalTraits <- comparative.data(data = DataSubset, 
                                         phy = TreeSubset, 
                                         names.col = "species", 
                                         vcv = TRUE,
                                         na.omit = TRUE, warn.dropped = TRUE)

#estimated pi vs speciation rate + traits
if(!paste0('global_gendivSpRateTraits_results_',repiTree,'.rds') %in% list.files(paste0(wd,'output'))){
  tryCatch({
    modpiSpRateTraits <- NULL
    modpiSpRateTraits <- pgls(log(EstPi) ~ log(tipRate) + log(BodyMassKg_notInputed) + log(geoArea_km2) + log(mean_temp) + log(litter_or_clutch_size_n) + log(GenerationLength_d), 
                              data = datacompGlobalTraits, lambda = 'ML')
    
    modpiTraits<- NULL
    modpiTraits <- pgls(log(EstPi) ~ log(BodyMassKg_notInputed) + log(geoArea_km2) + log(mean_temp) + log(litter_or_clutch_size_n) + log(GenerationLength_d), 
                        data = datacompGlobalTraits, lambda = 'ML')
    
    modSpRateTraits <- NULL
    modSpRateTraits <- pgls(log(tipRate) ~ log(BodyMassKg_notInputed) + log(geoArea_km2) + log(mean_temp) + log(litter_or_clutch_size_n) + log(GenerationLength_d), 
                            data = datacompGlobalTraits, lambda = 'ML')
  }, error=function(e){cat("Error in",conditionMessage(e), "\n")})
  
  pgls <- list(piSpRateTraits = modpiSpRateTraits, 
               piTraits = modpiTraits, 
               SpRateTraits = modSpRateTraits)
  
  saveRDS(pgls, paste0(wd, '7.pgls/outputs/global_gendivSpRateTraits_results_',repiTree,'.rds'))  ### kept in the cluster
}

#### Mutation rate ####
MutRateAll <- readRDS(paste0(folder_path, '6.mutationRate/output/mutationRate_cytb3rdcodonPAML.rds'))

## Join data 
gendivSpRateMutRate <- dplyr::select(gen.divCIc_sc, species, EstPi) %>%
  gendivSpRateTrait  %>%
  left_join(., MutRateAll, by = c('set','species')) %>%
  select(gen.divCIc_sc, species, set, EstPi, tipRate, time, GenerationLength_d, expNsub) %>% 
  mutate(timeyear = time * 1000000,
         GenerationLength_y = GenerationLength_d/365,
         timeGeneration = timeyear/GenerationLength_y,
         mutRate = expNsub / timeGeneration,
         Ne = EstPi / mutRate,
         mutRate_y = expNsub/timeyear) %>%
  drop_na()

## Prepare analyses sets only for global analysis - parallelise per tree

DataSubset <- gendivSpRateMutRate %>% 
  filter(set %in% repiTree) %>% 
  drop_na() 

TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
TreeSubset

datacompGlobalMutRate <- comparative.data(data = DataSubset, 
                                          phy = TreeSubset, 
                                          names.col = "species", 
                                          vcv = TRUE,
                                          na.omit = TRUE, warn.dropped = TRUE)

## Run mutation rate analyses

#mutation rate vs speciation rate
if(!paste0('global_gendivSpRateMutRate_results_',repiTree,'.rds') %in% list.files(paste0(wd,'output'))){
  tryCatch({
    modMutRateSpRate <- NULL
    modMutRateSpRate <- pgls(log(mutRate) ~ log(tipRate), data = datacompGlobalMutRate, lambda = 'ML')
    
    #Ne (pi/u) vs speciation rate
    modNeSpRate <- NULL
    modNeSpRate <- pgls(log(Ne) ~ log(tipRate), data = datacompGlobalMutRate, lambda = 'ML')
    
   }, error=function(e){cat("Error in",conditionMessage(e), "\n")})
  
  pgls <- list(MutRateSpRate = modMutRateSpRate, 
               NeSpRate = modNeSpRate)
  
  saveRDS(pgls, paste0(wd, '7.pgls/outputs/global_gendivSpRateMutRate_results_',repiTree,'.rds')) ### kept in the cluster
}

### in the cluster the script extractPGLS.R was used to get the summary of all PGLS analyses across all trees and tested models
