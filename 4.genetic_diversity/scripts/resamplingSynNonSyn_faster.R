library(seqinr)
library(bioseq)
library(tidyverse)
library(ape)
library(Rmisc)

args <- commandArgs()
print(args)
i <- as.numeric(gsub('--file=', '', args[6]))

###running in the cluster
setwd('/data/biodiv/asilva/genetic_diversity/')

seq_list <- list.files("/data/biodiv/asilva/genetic_diversity/cleanInFrameAlignments/", 
                       recursive = F, full.names = T)

###### functions ######

pi_per_site <- function(seqDNA){
  alleles<-table(seqDNA)
  alleles<-alleles[names(alleles) != "-"]
  n_alleles<-sum(alleles)
  return(1-sum(alleles*(alleles-1))/(n_alleles*(n_alleles-1)))
}
pi_estimator_fastaSites <- function(sequences, segSites) {
  pi_sites <- sapply(1:ncol(sequences), function(i) pi_per_site(sequences[,i]))
  pi_sites[is.na(pi_sites)] <- 0
  sum(pi_sites[segSites])/ncol(sequences)
}
theta_estimator_fastaSites <- function (seqDNA, segSites) {
  S <- length(segSites)
  sum_harmonic <- 0
  for (j in 1:ncol(seqDNA)){
    if (length(which(seqDNA[,j] %in% "-"))<nrow(seqDNA)-1){
      for (i in 1:(length(which(!seqDNA[,j] %in% "-"))-1)){
        sum_harmonic <- sum_harmonic + 1/i
      }
    } 
  }
  return(S/sum_harmonic)
}
seg_sites_bioseq <- function (bioseqdna) {
  segSitess <- c()
  for (i in 1:unique(seq_nchar(bioseqdna))){
    pp <- unlist(unique(lapply(strsplit(bioseqdna, ""), function(x) x[[i]])))
    
    if(length(pp[!pp %in% "-" & !pp %in% "N"]) > 1) 
      segSitess <- c(segSitess, i)
  }
  return(segSitess)
} ### too slow, use the other function but keep the code
nonSynSites_bioseq <- function (bioseqdna, segSites) { ###function for vertebrate mtDNA and all sequences set in frame 1
  trans <- seq_translate(bioseqdna, code = 2, codon_frame = 1)
  aaPoli <- c() ### get aminoacids positions with substitutions
  for (i in 1:unique(seq_nchar(trans))){
    pp <- unlist(unique(lapply(strsplit(trans, ""), function(x) x[[i]])))
    
    if(length(pp[!pp %in% "X"]) > 1)
      aaPoli <- c(aaPoli, i)
  }
  aaPoli3range <- c(aaPoli *3-2, aaPoli *3-1, aaPoli *3) ### get all positions that are in codons with substitutions
  nonSynSites <- segSites[which(segSites %in% aaPoli3range)] ### get only the segregating sites in the DNA from all the nonsyn codon sites
  return(nonSynSites)
}
seg_sites_fasta <- function (x) {
  S <- which(sapply(1:ncol(x), function(i) length(which(!unique(x[,i]) %in% c("-","N","n"))))>1)
  return(S)
}

###### Get resampled genetic diversity estimates ######
ReData <- data.frame()
seqDNA <- read.dna(seq_list[i],format="fasta", as.character = T)

sp <- tools::file_path_sans_ext(basename(seq_list[i]))
ind <- nrow(seqDNA)

if(ind > 4){  
  
  for(j in 1:100){
    setseqDNA <- seqDNA[sample(ind, 4),]
    
    setdnabioseq <- read_fasta(seq_list[i])[rownames( setseqDNA)]
    
    allsegSites <- seg_sites_fasta(setseqDNA)
    nonSynSites <- nonSynSites_bioseq(setdnabioseq, allsegSites)
    SynSites <- allsegSites[!allsegSites %in% nonSynSites]
    
    ReData <- bind_rows(ReData, data.frame(species = sp, rep = j,
                                           EstPiTotal = pi_estimator_fastaSites(setseqDNA,allsegSites),
                                           EstPiSyn = pi_estimator_fastaSites(setseqDNA,SynSites),
                                           EstThetaTotal = theta_estimator_fastaSites(setseqDNA,allsegSites),
                                           EstThetaSyn = theta_estimator_fastaSites(setseqDNA,SynSites)))
  }
}

saveRDS(ReData, paste0("/data/biodiv/asilva/genetic_diversity/output/GenDivSyn_resample4ind_",sp,".rds"))

#### after sending calculations per species:
library(tidyverse)
dir_path <- '/data/biodiv/asilva/rerunAnalyses/genetic_diversity/output'
df <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE) %>% 
  map(readRDS) %>% 
  bind_rows()
saveRDS(df,'/data/biodiv/asilva/genetic_diversity/GenDivSyn_resample4ind_output.rds')

dfSum <- df %>% 
  group_by(species) %>%
  dplyr::summarise(subPiTotal_mean = mean(EstPiTotal,na.rm = TRUE),
                   subPiTotal_min = min(EstPiTotal,na.rm = TRUE),
                   subPiTotal_max = max(EstPiTotal,na.rm = TRUE),
                   subPiSyn_mean = mean(EstPiSyn,na.rm = TRUE),
                   subPiSyn_sd = sd(EstPiSyn, na.rm = TRUE),
                   subPiSyn_min = min(EstPiSyn,na.rm = TRUE),
                   subPiSyn_max = max(EstPiSyn,na.rm = TRUE),
                   subThetaTotal_mean = mean(EstThetaTotal,na.rm = TRUE),
                   subThetaTotal_min = min(EstThetaTotal,na.rm = TRUE),
                   subThetaTotal_max = max(EstThetaTotal,na.rm = TRUE),
                   subThetaSyn_mean = mean(EstThetaSyn,na.rm = TRUE),
                   subThetaSyn_sd = sd(EstThetaSyn, na.rm = TRUE),
                   subThetaSyn_min = min(EstThetaSyn,na.rm = TRUE),
                   subThetaSyn_max = max(EstThetaSyn,na.rm = TRUE)) 
#write.table(dfSum, '/data/biodiv/asilva/genetic_diversity/GenDiv_resampled4ind.txt', row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
#dfSum <- read.delim('/Users/acas/Dropbox/Post-docs/Morlon_Lab/analyses/manuscript_scripts_data/genetic_diversity/GenDiv_resampled4ind.txt', header = TRUE, sep = "\t")

gen.divAll <- left_join(gen.div, dfSum, by = 'species') %>%
  mutate(outPiTotal = case_when(EstPiTotal > subPiTotal_max | EstPiTotal < subPiTotal_min ~ 'out',
                           TRUE ~ 'in'),
         outPiSyn = case_when(EstPiSyn > subPiSyn_max | EstPiSyn < subPiSyn_min ~ 'out',
                                TRUE ~ 'in'),
         outThetaTotal = case_when(EstThetaTotal > subThetaTotal_max | EstThetaTotal < subThetaTotal_min ~ 'out',
                                TRUE ~ 'in'),
         outThetaSyn = case_when(EstThetaSyn > subThetaSyn_max | EstThetaSyn < subThetaSyn_min ~ 'out',
                                TRUE ~ 'in'))

write.table(gen.divAll, '4.genetic_diversity/outputs/GenDiv_SynNonSyn_resampled4ind.txt', row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

