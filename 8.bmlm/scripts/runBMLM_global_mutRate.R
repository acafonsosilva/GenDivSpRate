unlink(".RData")
library(tidyverse)
library(ape)
library(Rmisc)
Sys.setenv(USE_CXX14 = 1)

args <- commandArgs()
print(args)
a <- as.numeric(args[4])
Sys.sleep(a*10)

#### load data
folder_path <- '/data/biodiv/asilva/SuperCrunchClean/'
wd <- paste0(folder_path, 'bmlm/global_mutRate/output/')

MutRateAll <- readRDS(paste0(folder_path, 'mutationRate/output/mutationRate_cytb3rdcodonPAML.rds'))

## discard species without genetic diversity and species which either mean subsampled pi or theta are outside the range of 1000 subsamples to 5 individuals
gen.divCIc_sc <- read.delim(paste0(folder_path, 'genetic_diversity/GenDiv_subsample5Ind.txt')) %>%
  filter(EstPi > 0, !outPi %in% 'out', !outTheta %in% 'out') %>% ### this way doesn't exclude species with 5 individuals
  select(species, EstPi, Nind)

##load phylogenies and speciation rates from 100 posterior trees + MCC tree

load(paste0(folder_path,'speciation_rate/output/upham_4064sp_FR_MCCposterior100.rdata')) #loads TreeSet
treeN <- names(TreeSet)[a]

clades <- read.delim(paste0(folder_path, 
                            'speciation_rate/input/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt')) %>%
  drop_na(PC) %>%
  mutate(species = word(tiplabel, 1,2, sep = "_"),
         clades = sub("^PC\\d+_",  "", PC)) %>%
  dplyr::select(species, clades) 

spRate <- read.delim(paste0(folder_path,'speciation_rate/output/MCCposterior100_tipRate.txt'),
                     stringsAsFactors = FALSE) %>% 
  dplyr::rename(set = treeN) %>% 
  left_join(., clades, by = 'species') 

traitData <- read.delim(paste0(folder_path,'traits/matchedTraits.txt'), 
                        stringsAsFactors = FALSE) %>%
  select(species, GenerationLength_d) %>%
  drop_na()

#### Join data ####
gendivSpRateMutRate <- dplyr::select(gen.divCIc_sc, species, EstPi) %>%
  left_join(., spRate, by = 'species') %>%
  left_join(., MutRateAll, by = c('set','species')) %>%
  left_join(., traitData, by = 'species') %>%
  mutate(timeyear = time * 1000000,
         GenerationLength_y = GenerationLength_d/365,
         timeGeneration = timeyear/GenerationLength_y,
         mutRate = expNsub / timeGeneration,
         Ne = EstPi / mutRate) %>%
  drop_na()

### For each set
DataSubset <- gendivSpRateMutRate %>%
  dplyr::filter(set %in% treeN) %>%
  drop_na() %>%
  droplevels()
dim(DataSubset)

t <- TreeSet[[treeN]]
TreeSubset <- drop.tip(t, as.character(t$tip.label[which(!t$tip.label %in% DataSubset$species)]))
TreeSubset
InvTree <- vcv.phylo(TreeSubset)
colnames(DataSubset)[1] <- "phylo"

#### https://github.com/paul-buerkner/brms/issues/114
library(brms)
library(rstan)
library(cmdstanr)
rstan_options(javascript=FALSE)
library(V8)

# mutation rate vs speciation rate
dir.create(paste0(wd,'MutRateSpRate'))
brm(log(mutrate) ~ log(tipRate) + (1|mutclade) + (1|gr(phylo, cov = InvTree)),
    data = DataSubset,
    family = gaussian(),
    data2 = list(InvTree = InvTree),
    file = paste0(wd,'MutRateSpRate/brms_MutRateSpRate_randomSlope_',treeN,".rds"), refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
    chains = 2, cores = 2, backend = "cmdstanr", threads = threading(2),
    iter = 4000, warmup = 1000, thin = 1)

#Ne vs speciation rate
dir.create(paste0(wd,'NeSpRate'))
brm(log(Ne) ~ log(tipRate) + (1|mutclade) + (1|gr(phylo, cov = InvTree)),
    data = DataSubset,
    family = gaussian(),
    data2 = list(InvTree = InvTree),
    file = paste0(wd,'NeSpRate/brms_NeSpRate_randomSlope_',treeN,".rds"), refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
    chains = 2, cores = 2, backend = "cmdstanr", threads = threading(2),
    iter = 4000, warmup = 1000, thin = 1)

