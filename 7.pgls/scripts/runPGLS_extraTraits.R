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
gen.divCIc_sc <- read.delim(paste0(folder_path, '4.genetic_diversity/outputs/GenDiv_SynNonSyn_resampled4ind.txt')) %>%
  filter(EstPiSyn > 0, outPiSyn %in% 'in', outThetaSyn %in% 'in') ### this way doesn't exclude species with 5 individuals

##load phylogenies and speciation rates from 100 posterior trees + MCC tree
load(paste0(folder_path,'5.speciation_rate/outputs/upham_4064sp_FR_MCCposterior100.rdata')) #loads TreeSet
repiTree <- names(TreeSet)[repi]
print(repiTree)
t <- TreeSet[[repiTree]]

clades <- read.delim(paste0(folder_path, 
                            '5.speciation_rate/inputs/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt')) %>%
  drop_na(PC) %>%
  mutate(species = word(tiplabel, 1,2, sep = "_"),
         clades = sub("^PC\\d+_",  "", PC)) %>%
  dplyr::select(species, clades) 

spRate <- read.delim(paste0(folder_path,'5.speciation_rate/outputs/MCCposterior100_tipRate.txt')) %>% 
  rename(set = treeN ) %>% 
  filter(set %in% names(TreeSet)[repi]) %>%
  left_join(., clades, by = 'species') 

traitData <- read.delim(paste0(folder_path,'otherData/traits/matchedTraits_v3.txt'), stringsAsFactors = FALSE) %>%
  mutate(mean_temp = mean_temp + abs(min(mean_temp, na.rm = T)) + 1,
         latitude_mean = abs(latitude_mean)) 

gendivSpRateTrait <- dplyr::select(gen.divCIc_sc, species, EstPiSyn) %>%
  left_join(., spRate, by = 'species') %>%
  left_join(., traitData, by = 'species')

#### Prepare trait analyses sets only for global analysis - parallelise per tree

## per trait analysis was done with the subset of species that includes geoArea (1119 species) but it 
##should have been done with the subset without (1131 species) since in the end this was excluded from the article given 
##we wanted to avoid confounding variables and we were requested to add latitude as a covariate, which was highly correlated with geographic area
## since this is just for supplementary figure and it's not a big difference in number of species these analyses were not repeated

DataSubset <- gendivSpRateTrait %>%
  select(species, set, EstPiSyn, tipRate, BodyMassKg_notInputed, #geoArea_km2, 
         mean_temp,
         litter_or_clutch_size_n, GenerationLength_d, latitude_mean) %>%
  filter(set %in% repiTree) %>%
  drop_na()

TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
TreeSubset

mod21 <- gls(log(tipRate) ~ log(BodyMassKg_notInputed), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML")  
print("mod21 done")

mod22 <- gls(log(tipRate) ~ log(geoArea_km2), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod22 done")

mod23 <- gls(log(tipRate) ~ log(mean_temp), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod23 done")

mod24 <- gls(log(tipRate) ~ log(latitude_mean), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod24 done")

mod25 <- gls(log(tipRate) ~ log(litter_or_clutch_size_n), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod25 done")

mod26 <- gls(log(tipRate) ~ log(GenerationLength_d), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod26 done")

mod10 <- gls(log(EstPiSyn) ~ 1, 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML")  
print("mod10 done")

mod11 <- gls(log(EstPiSyn) ~ log(BodyMassKg_notInputed), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML")  
print("mod11 done")

mod12 <- gls(log(EstPiSyn) ~ log(geoArea_km2), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod12 done")

mod13 <- gls(log(EstPiSyn) ~ log(mean_temp), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod13 done")

mod14 <- gls(log(EstPiSyn) ~ log(latitude_mean), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod14 done")

mod15 <- gls(log(EstPiSyn) ~ log(litter_or_clutch_size_n), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod15 done")

mod16 <- gls(log(EstPiSyn) ~ log(GenerationLength_d), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML")

print("mod16 done")

mod31 <- gls(log(EstPiSyn) ~ log(tipRate) + log(BodyMassKg_notInputed), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML")  
print("mod31 done")

mod32 <- gls(log(EstPiSyn) ~ log(tipRate)  + log(geoArea_km2), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 

print("mod32 done")

mod33 <- gls(log(EstPiSyn) ~ log(tipRate)  + log(mean_temp), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod33 done")

mod34 <- gls(log(EstPiSyn) ~ log(tipRate)  + log(latitude_mean), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod34 done")

mod35 <- gls(log(EstPiSyn) ~ log(tipRate)  + log(litter_or_clutch_size_n), 
             data = DataSubset ,
             correlation = corPagel(1, phy = TreeSubset, 
                                    form =~species), method = "ML") 
print("mod35 done")

pglsList <- list(mod10 = mod10,
                 mod11 = mod11,
                 mod12 = mod12,
                 mod13 = mod13,
                 mod14 = mod14,
                 mod15 = mod15,
                 mod16 = mod16,
                 mod21 = mod21,
                 mod22 = mod22,
                 mod23 = mod23,
                 mod24 = mod24,
                 mod25 = mod25,
                 mod26 = mod26,
                 mod31 = mod31,
                 mod32 = mod32,
                 mod33 = mod33,
                 mod34 = mod34,
                 mod35 = mod35)


saveRDS(pglsList, paste0(folder_path,"/pgls/output_extra/extraTraits_", repiTree,".rds"))

# # Set the directory path
# directory <- "/data/biodiv/asilva/rerunAnalyses/pgls/output_extra/"
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
# saveRDS(final_df,'/data/biodiv/asilva/rerunAnalyses/pgls/extraTraits_modelOutputs.rds')