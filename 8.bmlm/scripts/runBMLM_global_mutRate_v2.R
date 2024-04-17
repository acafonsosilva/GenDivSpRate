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
#Sys.sleep(a*30)

#### load data
folder_path <- '/data/biodiv/asilva/rerunAnalyses/'
wd <- paste0(folder_path, 'bmlm/global_mutRate/output/')

MutRateAll <- readRDS(paste0(folder_path, 'mutation_rate/output/mutationRate_cytb3rdcodonPAML.rds')) %>%
  select(-mutrate)

head(MutRateAll)
## discard species without genetic diversity and species which either mean subsampled pi or theta are outside the range of 1000 subsamples to 5 individuals
gen.divCIc_sc <- read.delim(paste0(folder_path, 'genetic_diversity/GenDiv_SynNonSyn_resampled4ind.txt')) %>%
  filter(EstPiSyn > 0, outPiSyn %in% 'in', outThetaSyn %in% 'in') %>%
  select(species, EstPiSyn, Nind)
head(gen.divCIc_sc)
##load phylogenies and speciation rates from 100 posterior trees + MCC tree

load(paste0(folder_path,'speciation_rate/output/upham_4064sp_FR_MCCposterior100.rdata')) #loads TreeSet
print(paste0(folder_path,'speciation_rate/output/upham_4064sp_FR_MCCposterior100.rdata'))
print(names(TreeSet))
treeN <- names(TreeSet)[a]
print(treeN)
spRate <- read.delim(paste0(folder_path,'speciation_rate/output/MCCposterior100_tipRate.txt'),
                     stringsAsFactors = FALSE) %>% 
  dplyr::rename(set = treeN) 
head(spRate)
traitData <- read.delim(paste0(folder_path,'traits/matchedTraits_v2.txt'), stringsAsFactors = FALSE) %>%
  select(species, GenerationLength_d) 
head(traitData)
DataSubset <- gen.divCIc_sc %>%
  left_join(., spRate, by = 'species', multiple = "all") %>%
  left_join(., traitData, by = 'species') %>%
  left_join(., MutRateAll, by = c('set','species')) %>%
  mutate(timeyear = time * 1000000,
         GenerationLength_y = GenerationLength_d/365,
         mutRate = (expNsub * GenerationLength_y) / timeyear,
         Ne = EstPiSyn / mutRate,
         mutRate_y = expNsub/timeyear) %>%
  select(species, set, mutclade, mutRate, Ne, tipRate) %>%
  filter(set %in% treeN) %>%
  drop_na() %>%
  droplevels()

dim(DataSubset)

tree <- TreeSet[[treeN]]
tree
TreeSubset <- drop.tip(tree, as.character(tree$tip.label[which(!tree$tip.label %in% DataSubset$species)]))
TreeSubset
InvTree <- vcv.phylo(TreeSubset)
colnames(DataSubset)[1] <- "phylo"

#### https://github.com/paul-buerkner/brms/issues/114
library(brms)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(javascript=FALSE)
library(cmdstanr)
set_cmdstan_path('/users/biodiv/asilva/.cmdstan/cmdstan-2.31.0')

library(V8)

# mutation rate vs speciation rate
dir.create(paste0(wd,'MutRateSpRate'))
if(!paste0('brms_MutRateSpRate_randomSlope_',treeN,".rds") %in% list.files(paste0(wd,'MutRateSpRate/'))){
brm(log(mutRate) ~ log(tipRate) + (1|mutclade) + (1|gr(phylo, cov = InvTree)),
    data = DataSubset,
    family = gaussian(),
    data2 = list(InvTree = InvTree),
    file = paste0(wd,'MutRateSpRate/brms_MutRateSpRate_randomSlope_',treeN,".rds"), refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
    chains = 2, cores = 2, backend = "cmdstanr", threads = threading(2),
    iter = 4000, warmup = 1000, thin = 2)

}
#Ne vs speciation rate

dir.create(paste0(wd,'NeSpRate'))
if(!paste0('brms_NeSpRate_randomSlope_',treeN,".rds") %in% list.files(paste0(wd,'NeSpRate/'))){
brm(log(Ne) ~ log(tipRate) + (1|mutclade) + (1|gr(phylo, cov = InvTree)),
    data = DataSubset,
    family = gaussian(),
    data2 = list(InvTree = InvTree),
    file = paste0(wd,'NeSpRate/brms_NeSpRate_randomSlope_',treeN,".rds"), refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
    chains = 2, cores = 2, backend = "cmdstanr", threads = threading(2),
    iter = 4000, warmup = 1000, thin = 2)
}