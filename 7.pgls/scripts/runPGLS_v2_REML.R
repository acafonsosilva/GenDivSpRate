library(ape)
library(tidyverse)
library(nlme)
select = dplyr::select

args <- commandArgs()
print(args)
repi <- as.numeric(gsub('--file=', '', args[6]))
print(repi)

#### load data
folder_path <- '/data/biodiv/asilva/rerunAnalyses/'

## discard species without genetic diversity and species which either mean subsampled pi or theta are outside the range of 1000 subsamples to 5 individuals
gen.divCIc_sc <- read.delim(paste0(folder_path, 'genetic_diversity/GenDiv_SynNonSyn_resampled4ind.txt')) %>%
  filter(EstPiSyn > 0, outPiSyn %in% 'in', outThetaSyn %in% 'in') ### this way doesn't exclude species with 5 individuals

##load phylogenies and speciation rates from 100 posterior trees + MCC tree

load(paste0(folder_path,'speciation_rate/output/upham_4064sp_FR_MCCposterior100.rdata')) #loads TreeSet
repiTree <- names(TreeSet)[repi]
print(repiTree)
t <- TreeSet[[repiTree]]

clades <- read.delim(paste0(folder_path, 
                            'speciation_rate/input/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt')) %>%
  drop_na(PC) %>%
  mutate(species = word(tiplabel, 1,2, sep = "_"),
         clades = sub("^PC\\d+_",  "", PC)) %>%
  dplyr::select(species, clades) 

spRate <- read.delim(paste0(folder_path,'speciation_rate/output/MCCposterior100_tipRate.txt')) %>% 
  rename(set = treeN ) %>% 
  filter(set %in% names(TreeSet)[repi]) %>%
  left_join(., clades, by = 'species') 

gendivSpRate <- dplyr::select(gen.divCIc_sc, species, EstPiSyn, EstThetaSyn) %>%
  left_join(., spRate, by = 'species') 

gendivSpRateOther <- dplyr::select(gen.divCIc_sc, species, EstPiTotal, EstPiNonSyn, 
                                   EstThetaTotal, EstThetaNonSyn) %>%
  filter(EstPiTotal > 0, EstPiNonSyn > 0, EstThetaTotal > 0, EstThetaNonSyn > 0) %>%
  left_join(., spRate, by = 'species') 

gendivSpRateSub <- dplyr::select(gen.divCIc_sc, species, subPiSyn_mean, subThetaSyn_mean) %>%
  left_join(., spRate, by = 'species')

#### global ####

## Prepare Estimated dataset independently of subsampled dataset with species with more than 5 individuals and resampled 1000x to 5
DataSubset <- gendivSpRate %>%
  filter(set %in% repiTree) %>%
  droplevels()
dim(DataSubset)

TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
TreeSubset

DataSubsetOther <- gendivSpRateOther %>%
  filter(set %in% repiTree) %>%
  droplevels()

TreeSubsetOther <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubsetOther$species)]))
TreeSubsetOther

DataSubsetSub <- gendivSpRateSub %>%
  filter(set %in% repiTree) %>%
  droplevels()
dim(DataSubsetSub)

TreeSubsetSub <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubsetSub$species)]))
TreeSubset


## Run global
modpglsEstPiTotal <- NULL
modpglsEstPiTotal <- gls(log(EstPiTotal) ~ log(tipRate), 
                         data = DataSubsetOther,
                         correlation = corPagel(1, phy = TreeSubsetOther, 
                                                form =~species), 
                         method = "REML")

modpglsEstPiNonSyn <- NULL
modpglsEstPiNonSyn <- gls(log(EstPiNonSyn) ~ log(tipRate), 
                          data = DataSubsetOther,
                          correlation = corPagel(1, phy = TreeSubsetOther, 
                                                 form =~species), 
                          method = "REML")

modpglsEstThetaTotal <- NULL
modpglsEstThetaTotal <- gls(log(EstThetaTotal) ~ log(tipRate), 
                            data = DataSubsetOther,
                            correlation = corPagel(1, phy = TreeSubsetOther, 
                                                   form =~species), 
                            method = "REML")

modpglsEstThetaNonSyn <- NULL
modpglsEstThetaNonSyn <- gls(log(EstThetaNonSyn) ~ log(tipRate), 
                             data = DataSubsetOther,
                             correlation = corPagel(1, phy = TreeSubsetOther, 
                                                    form =~species), 
                             method = "REML")

#estimated pi vs speciation rate at the tips
modpglsEstPiSyn <- NULL
modpglsEstPiSyn <- gls(log(EstPiSyn) ~ log(tipRate), 
                       data = DataSubset,
                       correlation = corPagel(1, phy = TreeSubset, 
                                              form =~species), 
                       method = "REML")

#mean subsampled pi vs speciation rate at the tips
modpglssubPiSyn_mean <- NULL
modpglssubPiSyn_mean <- gls(log(subPiSyn_mean) ~ log(tipRate), 
                            data = DataSubsetSub,
                            correlation = corPagel(1, phy = TreeSubsetSub, 
                                                   form =~species), 
                            method = "REML")

#estimated theta vs speciation rate at the tips
modpglsEstThetaSyn <- NULL
modpglsEstThetaSyn <- gls(log(EstThetaSyn) ~ log(tipRate), 
                          data = DataSubset,
                          correlation = corPagel(1, phy = TreeSubset, 
                                                 form =~species), 
                          method = "REML")

#mean subsampled theta vs speciation rate at the tips
modpglssubThetaSyn_mean <- NULL
modpglssubThetaSyn_mean <- gls(log(subThetaSyn_mean) ~ log(tipRate), 
                               data = DataSubsetSub,
                               correlation = corPagel(1, phy = TreeSubsetSub, 
                                                      form =~species), 
                               method = "REML")

## Run clade
##Prepare per clade analyses sets without traits - parallelise per clade

freqClade <- table(gendivSpRate[gendivSpRate$set %in% 'treeMCC','clades'])

pglsClade <- list()
for (clade in levels(as.factor(gendivSpRate$clades))){
  if(clade %in% names(freqClade[freqClade > 19])){ ##only for clades with at least 20 species
    
    dtset  <- gendivSpRate %>%
      filter(set %in% repiTree,  clades %in% clade) %>%
      droplevels()
    TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% dtset$species)]))
    
    dtsetSub <- gendivSpRateSub %>%
      filter(set %in% repiTree,  clades %in% clade) %>%
      droplevels()
    
    TreeSubsetSub <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% dtsetSub$species)]))
    
    #estimated pi vs speciation rate at the tips
    modpglsEstPiSync <- NULL
    modpglsEstPiSync <- gls(log(EstPiSyn) ~ log(tipRate), 
                            data = dtset,
                            correlation = corPagel(1, phy = TreeSubset, 
                                                   form =~species), 
                            method = "ML")
    
    # #mean subsampled pi vs speciation rate at the tips
    modpglssubPiSyn_meanc <- NULL
    modpglssubPiSyn_meanc <- gls(log(subPiSyn_mean) ~ log(tipRate), 
                                 data = dtsetSub,
                                 correlation = corPagel(1, phy = TreeSubsetSub, 
                                                        form =~species), 
                                 method = "ML")
    
    #estimated theta vs speciation rate at the tips
    modpglsEstThetaSync <- NULL
    modpglsEstThetaSyn <- gls(log(EstThetaSyn) ~ log(tipRate), 
                              data = dtset,
                              correlation = corPagel(1, phy = TreeSubset, 
                                                     form =~species), 
                              method = "ML")
    
    # #mean subsampled theta vs speciation rate at the tips
    modpglssubThetaSyn_meanc <- NULL
    modpglssubThetaSyn_meanc <- gls(log(subThetaSyn_mean) ~ log(tipRate), 
                                    data = dtsetSub,
                                    correlation = corPagel(1, phy = TreeSubsetSub, 
                                                           form =~species), 
                                    method = "ML")
    
    pglsClade[[clade]] <- list(EstPiSyn = modpglsEstPiSync,
                               SubPiSyn = modpglssubPiSyn_meanc,
                               EstThetaSyn = modpglsEstThetaSync,
                               SubThetaSyn = modpglssubThetaSyn_meanc)
    
    
  }
}


saveRDS(pglsClade, paste0(folder_path,'/pgls/output/clade_gendivSpRate_results_',repiTree,'.rds'))

#### Traits ####

traitData <- read.delim(paste0(folder_path,'traits/matchedTraits_v3.txt'), stringsAsFactors = FALSE) %>%
  mutate(mean_temp = mean_temp + abs(min(mean_temp, na.rm = T)) + 1,
         latitude_mean = abs(latitude_mean)) 

gendivSpRateTrait <- dplyr::select(gen.divCIc_sc, species, EstPiSyn) %>%
  left_join(., spRate, by = 'species') %>%
  left_join(., traitData, by = 'species')

#### Prepare trait analyses sets only for global analysis - parallelise per tree

DataSubset <- gendivSpRateTrait %>%
  select(species, set, EstPiSyn, tipRate, BodyMassKg_notInputed, mean_temp,
         litter_or_clutch_size_n, GenerationLength_d, latitude_mean) %>%
  filter(set %in% repiTree) %>%
  drop_na()

dim(DataSubset)

TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
TreeSubset

modpiSpRateTraits <- NULL
modpiSpRateTraits <- gls(log(EstPiSyn) ~ log(tipRate) + log(BodyMassKg_notInputed) + log(mean_temp) + log(latitude_mean) +
                           log(litter_or_clutch_size_n) + log(GenerationLength_d), 
                         data = DataSubset,
                         correlation = corPagel(1, phy = TreeSubset, 
                                                form =~species), 
                         method = "REML")

modpiTraits <- NULL
modpiTraits <- gls(log(EstPiSyn) ~ log(BodyMassKg_notInputed) + 
                     log(mean_temp) + log(latitude_mean) +
                     log(litter_or_clutch_size_n) +
                     log(GenerationLength_d), 
                   data = DataSubset ,
                   correlation = corPagel(1, phy = TreeSubset, 
                                          form =~species), 
                   method = "REML")

modSpRateTraits <- NULL
modSpRateTraits <- gls(log(tipRate) ~ log(BodyMassKg_notInputed) + log(mean_temp) + log(latitude_mean) +
                         log(litter_or_clutch_size_n) + log(GenerationLength_d), 
                       data = DataSubset ,
                       correlation = corPagel(1, phy = TreeSubset, 
                                              form =~species), 
                       method = "REML")


#### Mutation rate ####
MutRateAll <- readRDS(paste0(folder_path, 'mutation_rate/output/mutationRate_cytb3rdcodonPAML.rds'))

## Join data
gendivSpRateMutRate <- select(gendivSpRateTrait, species, set, EstPiSyn, tipRate, GenerationLength_d, geoArea_km2)   %>%
  #drop_na() %>%
  left_join(., MutRateAll, by = c('set','species')) %>%
  select(-mutrate) %>% ## mutation rate per Mya and not per year
  mutate(timeyear = time * 1000000,
         GenerationLength_y = GenerationLength_d/365,
         mutRate = (expNsub * GenerationLength_y) / timeyear,
         Ne = EstPiSyn / mutRate,
         mutRate_y = expNsub/timeyear,
         Ne_y = EstPiSyn / mutRate_y) %>%
  arrange(set)

## Prepare trait analyses sets only for global analysis - parallelise per tree
DataSubset1 <- gendivSpRateMutRate %>%
  select(species, set, mutRate, Ne,mutRate_y, Ne_y, tipRate) %>%
  filter(set %in% repiTree) %>%
  drop_na()


TreeSubset1 <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset1$species)]))
TreeSubset1

## Run mutation rate analyses

#mutation rate vs speciation rate
modMutRateSpRate <- NULL
modMutRateSpRate <- gls(log(mutRate) ~ log(tipRate), 
                        data = DataSubset1,
                        correlation = corPagel(1, phy = TreeSubset1, 
                                               form =~species), 
                        method = "REML")

#Ne (pi/u) vs speciation rate
modNeSpRate <- NULL
modNeSpRate <- gls(log(Ne) ~ log(tipRate), 
                   data = DataSubset1,
                   correlation = corPagel(1, phy = TreeSubset1, 
                                          form =~species), 
                   method = "REML")

#mutation rate vs speciation rate
modMutyRateSpRate <- NULL
modMutyRateSpRate <- gls(log(mutRate_y) ~ log(tipRate), 
                        data = DataSubset1,
                        correlation = corPagel(1, phy = TreeSubset1, 
                                               form =~species), 
                        method = "REML")

#Ne (pi/u) vs speciation rate
modNeySpRate <- NULL
modNeySpRate <- gls(log(Ne_y) ~ log(tipRate), 
                   data = DataSubset1,
                   correlation = corPagel(1, phy = TreeSubset1, 
                                          form =~species), 
                   method = "REML")

####
pgls <- list(EstPiTotal = modpglsEstPiTotal,
             EstPiNonSyn = modpglsEstPiNonSyn,
             EstThetaTotal = modpglsEstThetaTotal,
             EstThetaNonSyn = modpglsEstThetaNonSyn,
             EstPiSyn = modpglsEstPiSyn,
             SubPiSyn = modpglssubPiSyn_mean,
             EstThetaSyn = modpglsEstThetaSyn,
             SubThetaSyn = modpglssubThetaSyn_mean,
             piSpRateTraits = modpiSpRateTraits,
             piTraits = modpiTraits,
             SpRateTraits = modSpRateTraits,
             MutRateSpRate = modMutRateSpRate,
             NeSpRate = modNeSpRate,
             MutyRateSpRate = modMutyRateSpRate,
             NeySpRate = modNeySpRate)

saveRDS(pgls, paste0(folder_path,'/pgls/output/global_gendivSpRateALL_results_',repiTree,'.rds'))


# #####

extract_strings_in_parentheses <- function(text) {
  matches <- str_match_all(text, "\\((.*?)\\)")[[1]]
  if (is.null(matches))
    return(NA_character_)

  extracted_strings <- matches[, 2]
  return(paste(extracted_strings, collapse = " + "))
}

# Set the directory path
directory <- "/data/biodiv/asilva/rerunAnalyses/pgls/output/"

# Get the list of .rds files in the directory
file_list <- list.files(directory, pattern = "\\.rds$", full.names = TRUE)

# Initialize an empty list to store the results
results <- list()

# Loop through each .rds file
for (file_path in file_list) {
  # Read the .rds file as 'dt'
  dt <- readRDS(file_path)

  # Extract the 'TREENAME' from the file path
  TREENAME <- str_extract(file_path, "(?<=_results_)[^_]*?(?=\\.rds)")

  # Loop through each 'MODELS' and create a data frame
  model_data <- map_dfr(names(dt), function(model_name) {
    modelF <- as.character(dt[[model_name]]$call)[2]
    output <- data.frame(coef(summary(dt[[model_name]]))) %>%
      rownames_to_column('term')

    pgls_lambda <- attributes(dt[[model_name]]$apVar)$Pars['corStruct']

    tibble(set = TREENAME, analysis = model_name, modelF, df = dt[[model_name]]$dims$N,output, pgls_lambda)
  })

  # Append the model_data to the results list
  results[[TREENAME]] <- model_data
}

# Combine all data frames into a single data frame
final_df <- bind_rows(results)

# Print the final data frame
print(final_df)
saveRDS(final_df,'/data/biodiv/asilva/rerunAnalyses/pgls/gendivSpRate_PGLSresultsAll_REML.rds')
