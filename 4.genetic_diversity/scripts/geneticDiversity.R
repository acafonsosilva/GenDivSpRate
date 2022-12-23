library(tidyverse)
library(ape)
library(Rmisc)

wdpath <- '/data/biodiv/asilva/'

#list to per species alignments
seq_list <- list.files(paste0(wdpath,"3.supercrunch_family/Output_alignments"),
                       recursive = F, full.names = T)

###### functions to estimate genetic diversity sensu Ferretti et al. 2012######
pi_per_site <- function(seqDNA){
  alleles<-table(seqDNA)
  alleles<-alleles[names(alleles) != "-"]
  n_alleles<-sum(alleles)
  return(1-sum(alleles*(alleles-1))/(n_alleles*(n_alleles-1)))
}

##Luca Ferretti considers that including sites with just Ns or gaps in total sequence length could bias the estimation
cleanDNA <- function(sequences){ 
  sequences <- sequences[,sapply(1:ncol(sequences), 
                                 function(i) length(which(!sequences[,i] %in% c("-","N","n")))>1),drop=F]
  return(sequences)
}

pi_estimator_fasta <- function(sequences) {
  pi_sites <- sapply(1:ncol(sequences), function(i) pi_per_site(sequences[,i]))
  pi_sites[is.na(pi_sites)] <- 0
  sum(pi_sites)/ncol(sequences)
}

seg_sites_fasta <- function (x) {
  S <- which(sapply(1:ncol(x), function(i) length(which(!unique(x[,i]) %in% c("-","N","n"))))>1)
  return(S)
}

theta_estimator_fasta <- function (seqDNA) {
  S <- length(seg_sites_fasta(seqDNA))
  sum_harmonic <- 0
  for (j in 1:ncol(seqDNA)){if (length(which(seqDNA[,j] %in% "-"))<nrow(seqDNA)-1){
    for (i in 1:(length(which(!seqDNA[,j] %in% "-"))-1)){
      sum_harmonic <- sum_harmonic + 1/i}}}
  return(S/sum_harmonic)
}

###### Get genetic diversity estimates from raw data ######
gen.div <- data.frame()
ReList <- list()
for(i in 1:length(seq_list)){
  seqDNA <- cleanDNA(read.dna(seq_list[i],format="fasta",as.character = T))
  
  ### make all ambiguous sites into gaps so they are not accounted for
  seqDNA[which(!seqDNA %in% c("a","c","g","t","-"))] <- "-"
  
  sp <- gsub('_cytb','',tools::file_path_sans_ext(basename(seq_list[i])))
  ind <- nrow(seqDNA)
  
  ### only run for species with more than 4 individuals
  if(nrow(seqDNA) > 4){
    
    gen.div0 <- data.frame(species = sp,
                           Nind = ind,
                           size = ncol(seqDNA),
                           EstPi = pi_estimator_fasta(seqDNA),
                           EstTheta = theta_estimator_fasta(seqDNA))
    
    gen.divf <- data.frame()
    for (j in 1:3){
      subseqDNA <- seqDNA[,seq(j, ncol(seqDNA), 3)]
      
      gen.divf <- bind_rows(gen.divf, data.frame(species = sp,
                                                 frame = j,
                                                 S = length(seg_sites_fasta(subseqDNA)),
                                                 EstPi = pi_estimator_fasta(subseqDNA),
                                                 EstTheta = theta_estimator_fasta(subseqDNA)))
    }
    
    gen.div0 <- gen.divf %>% 
      pivot_wider(names_from = frame, values_from = c(S, EstPi, EstTheta)) %>%
      right_join(., gen.div0, by = 'species')
    
    ## sub-sample 1000x to 5 for species that have more than 5 individuals
    pi_resample1000 <- c()
    theta_resample1000 <- c()
    if(ind > 5){      
      for(j in 1:1000){

        setseqDNA <- cleanDNA(seqDNA[sample(ind, 5), ])
        pi_resample1000 <- c(pi_resample1000, pi_resample = pi_estimator_fasta(setseqDNA))
        theta_resample1000 <- c(theta_resample1000, theta_resample = theta_estimator_fasta(setseqDNA))
      }
      ReList[[sp]] <- list(pi_resample = pi_resample1000,
                           theta_resample = theta_resample1000)

      sub <- tibble(pi = pi_resample1000,
                        theta = theta_resample1000)

        resData <- sub %>%
        summarise(species = sp,
                  subPi_min = min(pi,na.rm = TRUE),
                  subPi_mean = mean(pi,na.rm = TRUE),
                  subPi_max = max(pi,na.rm = TRUE),
                  subPi_upper.CI = CI(pi, ci = 0.95)[1],
                  subPi_lower.CI = CI(pi, ci = 0.95)[3],
                  subPi_sd = sd(pi, na.rm = TRUE),
                  subPi_se = subPi_sd / sqrt(sum(!is.na(pi))),
                  subTheta_min = min(theta,na.rm = TRUE),
                  subTheta_mean = mean(theta,na.rm = TRUE),
                  subTheta_max = max(theta,na.rm = TRUE),
                  subTheta_upper.CI = CI(theta, ci = 0.95)[1],
                  subTheta_lower.CI = CI(theta, ci = 0.95)[3],
                  subTheta_sd = sd(theta, na.rm = TRUE),
                  subTheta_se  = subTheta_sd / sqrt(sum(!is.na(theta))))
    
      ### To check estimated genetic diversity is within range
      gen.div0 <- left_join(gen.div0, resData, by = 'species') %>%
        mutate(outPi = case_when(EstPi > subPi_max | EstPi < subPi_min ~ 'out',
                                 TRUE ~ 'in'),
               outTheta = case_when(EstTheta > subTheta_max | 
                                      EstTheta < subTheta_min ~ 'out',
                                    TRUE ~ 'in'))
    }
    gen.div <- bind_rows(gen.div,gen.div0)
  }
  print(paste0(i,"/",length(seq_list)))
}

saveRDS(ReList, paste0(wdpath,'4.genetic_diversity/outputs/GenDiv_subsample5Ind1000rep.rds'))
write.table(gen.div, paste0(wdpath,"4.genetic_diversity/outputs/GenDiv_subsample5Ind.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)