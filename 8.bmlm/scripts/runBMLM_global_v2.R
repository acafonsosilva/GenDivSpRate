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
wd <- paste0(folder_path, 'bmlm/global/output/')
dir.create(wd)

## discard species without genetic diversity and species which either mean subsampled pi or theta are outside the range of 1000 subsamples to 5 individuals
gen.divCIc_sc <- read.delim(paste0(folder_path, 'genetic_diversity/GenDiv_SynNonSyn_resampled4ind.txt')) %>%
  filter(EstPiSyn > 0, outPiSyn %in% 'in', outThetaSyn %in% 'in') 

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

gendivSpRate <- gen.divCIc_sc %>%
  left_join(., spRate, by = 'species', multiple = "all") %>%
  select(treeN, clades, species, EstPiSyn, subPiSyn_mean, subPiSyn_sd, EstThetaSyn, subThetaSyn_mean, subThetaSyn_sd, tipRate)

### For each set
DataSubset <- gendivSpRate %>%
  filter(treeN %in% set) %>%
  droplevels() 
dim(DataSubset)

t <- TreeSet[[set]]
TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
TreeSubset
InvTree <- vcv.phylo(TreeSubset)
colnames(DataSubset)[3] <- "phylo"

library(brms)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(javascript=FALSE)
library(cmdstanr)
set_cmdstan_path('/users/biodiv/asilva/.cmdstan/cmdstan-2.31.0')
library(V8)

#### https://github.com/paul-buerkner/brms/issues/114
if(!paste0('brms_EstPiSpRate_correlatedInterSlope_',set,".rds") %in% list.files(paste0(wd,'estPiSyn/'))){
dir.create(paste0(wd,'EstPiSyn'))
brm(log(EstPiSyn) ~ log(tipRate) + (log(tipRate)|clades) + (1|gr(phylo, cov = InvTree)),
    data = DataSubset,
    family = gaussian(),
    save_mevars = TRUE,
    data2 = list(InvTree = InvTree),
    file = paste0(wd,'estPiSyn/brms_EstPiSpRate_correlatedInterSlope_',set,".rds"), refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
    chains = 2, cores = 2,
    iter = 4000, warmup = 1000, thin = 2)
}

if(!paste0('brms_SubPiSpRate_correlatedInterSlope_',set,".rds") %in% list.files(paste0(wd,'subPiSyn/'))){
dir.create(paste0(wd,'SubPiSyn'))
brm(log(subPiSyn_mean) | se(subPiSyn_sd, sigma = TRUE) ~ log(tipRate) + (log(tipRate)|clades) + (1|gr(phylo, cov = InvTree)),
    data = DataSubset,
    family = gaussian(),
    save_mevars = TRUE,
    data2 = list(InvTree = InvTree),
    file = paste0(wd,'subPiSyn/brms_SubPiSpRate_correlatedInterSlope_',set,".rds"), refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
    chains = 2, cores = 2,
    iter = 4000, warmup = 1000, thin = 2)
}

if(!paste0('brms_EstThetaSpRate_correlatedInterSlope_',set,".rds") %in% list.files(paste0(wd,'estThetaSyn/'))){
dir.create(paste0(wd,'EstThetaSyn'))
brm(log(EstThetaSyn) ~ log(tipRate) + (log(tipRate)|clades) + (1|gr(phylo, cov = InvTree)),
    data = DataSubset,
    family = gaussian(),
    save_mevars = TRUE,
    data2 = list(InvTree = InvTree),
    file = paste0(wd,'estThetaSyn/brms_EstThetaSpRate_correlatedInterSlope_',set,".rds"), refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
    chains = 2, cores = 2,
    iter = 4000, warmup = 1000, thin = 2)
}

if(!paste0('brms_SubThetaSpRate_correlatedInterSlope_',set,".rds") %in% list.files(paste0(wd,'subThetaSyn/'))){
dir.create(paste0(wd,'SubThetaSyn'))
brm(log(subThetaSyn_mean) | se(subThetaSyn_sd, sigma = TRUE) ~ log(tipRate) + (log(tipRate)|clades) + (1|gr(phylo, cov = InvTree)),
    data = DataSubset,
    family = gaussian(),
    save_mevars = TRUE,
    data2 = list(InvTree = InvTree),
    file = paste0(wd,'subThetaSyn/brms_SubThetaSpRate_correlatedInterSlope_',set,".rds"), refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
    chains = 2, cores = 2,
    iter = 4000, warmup = 1000, thin = 2)
}

