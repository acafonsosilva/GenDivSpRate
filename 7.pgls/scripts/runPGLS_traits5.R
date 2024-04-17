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

traitData <- read.delim(paste0(folder_path,'traits/matchedTraits_v2.txt'), stringsAsFactors = FALSE) %>%
  mutate(mean_temp = mean_temp + abs(min(mean_temp, na.rm = T)) + 1,
         latitude_mean = abs(latitude_mean)) 

gendivSpRateTrait <- dplyr::select(gen.divCIc_sc, species, EstPiSyn) %>%
  left_join(., spRate, by = 'species') %>%
  left_join(., traitData, by = 'species')

#### Prepare trait analyses sets only for global analysis - parallelise per tree

DataSubset <- gendivSpRateTrait %>%
  select(species, set, EstPiSyn, tipRate, BodyMassKg_notInputed, geoArea_km2, mean_temp,
         litter_or_clutch_size_n, GenerationLength_d, latitude_mean) %>%
  filter(set %in% repiTree) %>%
  drop_na()

TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
TreeSubset

modpiSpRateTraits <- gls(log(EstPiSyn) ~ log(tipRate) + log(BodyMassKg_notInputed) + 
                           log(mean_temp) + log(latitude_mean) +
                           log(litter_or_clutch_size_n) + log(GenerationLength_d), 
                         data = DataSubset,
                         correlation = corPagel(1, phy = TreeSubset, 
                                                form =~species), 
                         method = "ML")

modpiTraits <- gls(log(EstPiSyn) ~ log(BodyMassKg_notInputed) + 
                     log(mean_temp) + log(latitude_mean) +
                     log(litter_or_clutch_size_n) +
                     log(GenerationLength_d), 
                   data = DataSubset ,
                   correlation = corPagel(1, phy = TreeSubset, 
                                          form =~species), 
                   method = "ML")

modSpRateTraits <- gls(log(tipRate) ~ log(BodyMassKg_notInputed) + 
                         log(mean_temp) + log(latitude_mean) +
                         log(litter_or_clutch_size_n) + 
                         log(GenerationLength_d), 
                       data = DataSubset ,
                       correlation = corPagel(1, phy = TreeSubset, 
                                              form =~species), 
                       method = "ML")


pglsList <- list(piSpRateTraits5 = modpiSpRateTraits,
                 piTraits5 = modpiTraits,
                 SpRateTraits5 = modSpRateTraits)


saveRDS(pglsList, paste0(folder_path,"/pgls/output_traits5/gendivSpRate_traits5_", repiTree,".rds"))

# # Set the directory path
# directory <- "/data/biodiv/asilva/rerunAnalyses/pgls/output_traits5/"
# 
# # Get the list of .rds files in the directory
# file_list <- list.files(directory, pattern = "\\.rds$", full.names = TRUE)
# 
# # Initialize an empty list to store the results
# results <- list()
# 
# # Loop through each .rds file
# for (file_path in file_list) {
#   # Read the .rds file as 'dt'
#   dt <- readRDS(file_path)
# 
#   # Extract the 'TREENAME' from the file path
#   TREENAME <- sub(".*/extraTraits_(.*?)\\.rds", "\\1", file_path)
# 
#   # Loop through each 'MODELS' and create a data frame
#   model_data <- map_dfr(names(dt), function(model_name) {
#     modelF <- as.character(dt[[model_name]]$call)[2]
#     output <- list(coef(summary(dt[[model_name]])))
# 
#     tibble(set = TREENAME, model = model_name, modelF, output)
#   })
# 
#   # Append the model_data to the results list
#   results[[TREENAME]] <- model_data
# }
# 
# # Combine all data frames into a single data frame
# final_df <- bind_rows(results)
# 
# # Print the final data frame
# print(final_df)
# saveRDS(final_df,'/data/biodiv/asilva/rerunAnalyses/pgls/traits5_modelOutputs.rds')