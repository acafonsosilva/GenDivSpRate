unlink(".RData")
library(tidyverse)
library(ape)

select = dplyr::select
filter=dplyr::filter
mutate=dplyr::mutate
Sys.setenv(USE_CXX14 = 1)

args <- commandArgs()
print(args)
a <- as.numeric(args[4])
#a <- as.numeric(gsub('--file=', '', args[6]))
Sys.sleep(a*10)

#### load data
folder_path <- '/data/biodiv/asilva/rerunAnalyses/'

wd <- paste0(folder_path, 'bmlm/global_traits/output/')

traitData <- read.delim(paste0(folder_path,'traits/matchedTraits_v2.txt'), stringsAsFactors = FALSE) %>%
  mutate(mean_temp = mean_temp + abs(min(mean_temp, na.rm = T)) + 1) %>% ### to not exclude species with negative mean temperatures when log transform 
  select(-DispDist_notInputed)

## discard species without genetic diversity and species which either mean subsampled pi or theta are outside the range of 1000 subsamples to 5 individuals
gen.divCIc_sc <- read.delim(paste0(folder_path, 'genetic_diversity/GenDiv_SynNonSyn_resampled4ind.txt')) %>%
  filter(EstPiSyn > 0, outPiSyn %in% 'in', outThetaSyn %in% 'in') %>%
  select(species, EstPiSyn, Nind)

##load phylogenies and speciation rates from 100 posterior trees + MCC tree

load(paste0(folder_path,'speciation_rate/output/upham_4064sp_FR_MCCposterior100.rdata')) #loads TreeSet
set <- names(TreeSet)[a]

clades <- read.delim(paste0(folder_path, 
                            'speciation_rate/input/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt')) %>%
  drop_na(PC) %>%
  mutate(species = word(tiplabel, 1,2, sep = "_"),
         clades = sub("^PC\\d+_",  "", PC)) %>%
  dplyr::select(species, clades) 

spRate <- read.delim(paste0(folder_path,'speciation_rate/output/MCCposterior100_tipRate.txt')) %>% 
  left_join(., clades, by = 'species') 

gendivSpRateTrait <- gen.divCIc_sc %>%
  left_join(., spRate, by = 'species', multiple = "all") %>%
  left_join(., traitData, by = 'species')

### For each set
DataSubset <- gendivSpRateTrait %>%
  filter(treeN %in% set) %>%
  drop_na() 
dim(DataSubset)

tree <- TreeSet[[set]]
TreeSubset <- drop.tip(tree, as.character(tree$tip.label[which(!tree$tip.label %in% DataSubset$species)]))
TreeSubset
InvTree <- vcv.phylo(TreeSubset)
colnames(DataSubset)[1] <- "phylo"

library(brms)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(javascript=FALSE)
library(cmdstanr)
set_cmdstan_path('/users/biodiv/asilva/.cmdstan/cmdstan-2.31.0')
library(V8)

#### https://github.com/paul-buerkner/brms/issues/114
# dir.create(paste0(wd,'piSpRateTraits'))
if(!paste0('brms_piSpRateTraits_randomSlope_',set,".rds") %in% list.files(paste0(wd,'piSpRateTraits/'))){
brm(log(EstPiSyn) ~ log(tipRate) + log(BodyMassKg_notInputed) + log(geoArea_km2) + log(mean_temp) + log(latitude_mean) + log(litter_or_clutch_size_n) + log(GenerationLength_d) + (1|clades) + (1|gr(phylo, cov = InvTree)),  # (log(tipRate)|clades) for better comparison with other models should only use patch clade for random slopes so not correlated with tipRate
    data = DataSubset,
    family = gaussian(),
    save_mevars = TRUE,
    data2 = list(InvTree = InvTree),
    file = paste0(wd,'piSpRateTraits/brms_piSpRateTraits_randomSlope_',set,".rds"), #refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
    chains = 2, cores = 1,
    iter = 4000, warmup = 1000, thin = 2)
}

# if(!paste0('brms_piTraits_randomSlope_',set,".rds") %in% list.files(paste0(wd,'piTraits/'))){
# dir.create(paste0(wd,'piTraits'))
# brm(log(EstPiSyn) ~ log(BodyMassKg_notInputed) + log(geoArea_km2) + log(mean_temp) + log(latitude_mean) + log(litter_or_clutch_size_n) + log(GenerationLength_d) + (1|clades) + (1|gr(phylo, cov = InvTree)),
#     data = DataSubset,
#     family = gaussian(),
#     save_mevars = TRUE,
#     data2 = list(InvTree = InvTree),
#     file = paste0(wd,'piTraits/brms_piTraits_randomSlope_',set,".rds"), refresh = 0,
#     control = list(adapt_delta = 0.99, max_treedepth = 20),
#     chains = 2, cores = 1,
#     iter = 4000, warmup = 1000, thin = 2)
# }

# if(!paste0('brms_SpRateTraits_randomSlope_',set,".rds") %in% list.files(paste0(wd,'SpRateTraits/'))){
# dir.create(paste0(wd,'SpRateTraits'))
# brm(log(tipRate) ~ log(BodyMassKg_notInputed) + log(geoArea_km2) + log(mean_temp) + log(latitude_mean) + log(litter_or_clutch_size_n) + log(GenerationLength_d) + (1|clades) + (1|gr(phylo, cov = InvTree)),
#     data = DataSubset,
#     family = gaussian(),
#     save_mevars = TRUE,
#     data2 = list(InvTree = InvTree),
#     file = paste0(wd,'SpRateTraits/brms_SpRateTraits_randomSlope_',set,".rds"), refresh = 0,
#     control = list(adapt_delta = 0.99, max_treedepth = 20),
#     chains = 2, cores = 1,
#     iter = 4000, warmup = 1000, thin = 2)
# }
