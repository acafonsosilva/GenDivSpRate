library(seqinr)
library(bioseq)
library(tidyverse)
library(ape)
library(Rmisc)


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

### to avoid any potential issues because of missing data and ambiguous sites any ambiguous codes were converted to gaps and alignments were exported so no changes happened between both reading dna functions
seq_list <- list.files('/Users/acas/Dropbox/Post-docs/Morlon_Lab/analyses/manuscript_scripts_data/superCrunch_family/Output_alignments/', full.names = TRUE)
for (i in 1:length(seq_list)){
  sp <- tools::file_path_sans_ext(basename(seq_list[i]))
  seqDNA <- read.dna(seq_list[i],format="fasta",as.character = T)
  ind <- nrow(seqDNA)
  if(ind > 4){
    seqDNA[which(!seqDNA %in% c("a","c","g","t","-"))] <- "-"
    
    segs <- length(seg_sites_fasta(seqDNA))
    outFasta <-paste0('/Users/acas/Dropbox/Post-docs/Morlon_Lab/analyses/manuscript_scripts_data/superCrunch_family/tempAlignments/', sp,".fasta")
    write.dna(seqDNA, outFasta, format = 'fasta')
    
    dnabioseq <- read_fasta(outFasta)
    segs2 <- length(seg_sites_bioseq(dnabioseq))
    
    if(segs != segs2){
      stop(paste0(sp, " has issue with cleaning sites"))
    }
    
    trans <- seq_translate(dnabioseq, code = 2, codon_frame = 1)
    if(any(str_detect(trans, "\\*"))){
      print(paste(i, "still not in frame"))
    } 
  }
}

#####
seq_list <- list.files('/Users/acas/Dropbox/Post-docs/Morlon_Lab/analyses/manuscript_scripts_data/superCrunch_family/tempAlignments/', full.names = TRUE)


gen.div <- data.frame()
for (i in 1:length(seq_list)){
  
  seqDNA <- read.dna(seq_list[i],format="fasta",as.character = T)
  
  ind <- nrow(seqDNA)
  if(ind > 4){
  sp <- tools::file_path_sans_ext(basename(seq_list[i]))
  dnabioseq <- read_fasta(seq_list[i])

  allsegSites <- seg_sites_fasta(seqDNA)
  nonSynSites <- nonSynSites_bioseq(dnabioseq, allsegSites)
  SynSites <- allsegSites[!allsegSites %in% nonSynSites]

    gen.div <- bind_rows(gen.div, data.frame(species = sp,
                                             Nind = ind,
                                             size = ncol(seqDNA),
                                             segSitestotal = length(allsegSites),
                                             synSites = length(SynSites),
                                             nonSynSites = length(nonSynSites),
                                             EstPiTotal = pi_estimator_fastaSites(seqDNA,allsegSites),
                                             EstPiSyn = pi_estimator_fastaSites(seqDNA,SynSites),
                                             EstPiNonSyn = pi_estimator_fastaSites(seqDNA,nonSynSites),
                                             EstThetaTotal = theta_estimator_fastaSites(seqDNA,allsegSites),
                                             EstThetaSyn = theta_estimator_fastaSites(seqDNA,SynSites),
                                             EstThetaNonSyn = theta_estimator_fastaSites(seqDNA,nonSynSites)))
    if(any(is.na(filter(gen.div2, species %in% sp)))){
      stop(paste0("problem with species ",sp))
    }
  }
}

write.table(gen.div, '/Users/acas/Dropbox/Post-docs/Morlon_Lab/analyses/manuscript_scripts_data/genetic_diversity/GenDiv_SynNonSyn.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
