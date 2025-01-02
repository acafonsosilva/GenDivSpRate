library(seqinr)
library(bioseq)
library(Rmisc)

seq_list <-  list.files('3.superCrunch_family/Output_alignments/', full.names = TRUE)

theo2021Match2Loc <- read_csv('otherData/Theodoridis2021/cytb_coordinates.csv') %>% 
  filter(seq %in% SeqData$seq) %>% 
  arrange(species) %>%
  group_by(species) %>%
  filter(n_distinct(x, y) >= 2)

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

gen.divG <- data.frame()
for (i in 1:length(seq_list)){
  
  sp <- tools::file_path_sans_ext(basename(seq_list[i]))
  if(sp %in% unique(theo2021Match2Loc$species)){
    
    seqDNA <- read.dna(seq_list[i], format="fasta",as.character = T)
    setIDs <- rownames(seqDNA)[str_detect(rownames(seqDNA), 
                                          paste(pull(filter(theo2021Match2Loc, species %in% sp), seq), collapse = "|"))]
    #  if(length(setIDs) > 4){
    dnabioseq <- read_fasta(seq_list[i])
    
    allsegSites <- seg_sites_fasta(seqDNA[setIDs,])
    nonSynSites <- nonSynSites_bioseq(dnabioseq[setIDs], allsegSites)
    SynSites <- allsegSites[!allsegSites %in% nonSynSites]
    
    gen.divG <- bind_rows(gen.divG, data.frame(species = sp,
                                               Nind = length(setIDs),
                                               segSitestotalG = length(allsegSites),
                                               synSitesG = length(SynSites),
                                               EstPiTotalG = pi_estimator_fastaSites(seqDNA[setIDs,],allsegSites),
                                               EstPiSynG = pi_estimator_fastaSites(seqDNA[setIDs,],SynSites)))
    #  }                                      
  }
}
write.table(gen.divG, '4.genetic_diversity/outuputs/GenDiv_geo.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)