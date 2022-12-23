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


gen.div <- data.frame()
for (i in 1:length(seq_list)){
  sp <- tools::file_path_sans_ext(basename(seq_list[i]))
  seqDNA <- cleanDNA(read.dna(seq_list[i],format="fasta",as.character = T))
  
  if(nrow(seqDNA) > 4){
    seqDNA[which(!seqDNA %in% c("a","c","g","t","-"))] <- "-"
    
    ind <- nrow(seqDNA)
    
    gen.divf <- data.frame()
    for (j in 1:3){
      subseqDNA <- seqDNA[,seq(j, ncol(seqDNA), 3)]
      
      gen.divf <- bind_rows(gen.divf, data.frame(species = sp,
                             frame = j,
                             S = length(seg_sites_fasta(subseqDNA)),
                             EstPi = pi_estimator_fasta(subseqDNA),
                             EstTheta = theta_estimator_fasta(subseqDNA)))
    }
    gen.div <- gen.divf %>% 
      pivot_wider(names_from = frame, values_from = c(S, EstPi, EstTheta)) %>%
      bind_rows(gen.div, .)
  }
}

write.table(gen.div, paste0(wdpath,"4.genetic_diversity/outputs/GenDiv_frame123.txt"), col.names = T, row.names = F, quote = F)



#### detect which species are out of frame to remove them from data set before analyses
library(bioseq)

checkFrame <- data.frame()
for (i in seq_list){
  dna <- read_fasta(i)
  
  if(length(dna) > 4){
    species <- tools::file_path_sans_ext(basename(i))
    trans <- seq_translate(dna, code = 2, codon_frame = 1)
    checkFrame <- bind_rows(checkFrame, data.frame(species = species,
                                                   check = ifelse(any(str_detect(trans, "\\*")),
                                                                  "notInFrame1", "inFrame1")))
    
    if(checkFrame[checkFrame$species %in% species,'check'] %in% 'notInFrame1'){
      dna <- read_fasta(i)
      trans <- seq_translate(dna, code = 2, codon_frame = 2)
      checkFrame[checkFrame$species %in% species,'check'] <- ifelse(any(str_detect(trans, "\\*")),
                                                                    "notInFrame12", "inFrame2")
    }
    if(checkFrame[checkFrame$species %in% species,'check'] %in% 'notInFrame12'){
      trans <- seq_translate(dna, code = 2, codon_frame = 3)
      checkFrame[checkFrame$species %in% species,'check'] <- ifelse(any(str_detect(trans, "\\*")),
                                                                    "notInFrame123", "inFrame3")
    }
  }
}

write.table(checkFrame, row.names = FALSE, quote = FALSE, sep = '\t',
            paste0(folder_path, '4.genetic_diversity/outputs/speciesFrame.txt'))