### Ana C. Afonso Silva script and not a SuperCrunch file

##### PACKAGES ####
library(seqinr)
library(tidyverse)
library(taxize)
library(rentrez)

set_entrez_key("78915277c9890b626a8b082c889ceca86a08")

#### read argument to define clade ####
args <- commandArgs()
clade <- args[6]
print(clade)

folder_family_path <- args[7]

# List matching extracted GenBank sequence to Upham tree
inputSeq <- read.delim(paste0(
  str_remove(folder_family_path, 
             paste0('/families/',clade)),'genBank_metadata_complete_matchedUpham.txt'))

big_fasta_path <- list.files(paste0(folder_family_path, "/5_Align/Alignments-MAFFT/"), pattern = "MAFFT_Aligned.fasta", full.names = TRUE)

big_fasta <- read.fasta(big_fasta_path, seqtype = c("DNA"),
                        seqonly = FALSE, forceDNAtolower = FALSE, set.attributes = FALSE)

new_folder_species <- paste0(folder_family_path, "/6_Per_species_fasta_files")
dir.create(new_folder_species)

inputSeqFam <- inputSeq %>%
  filter(accession %in% names(big_fasta))

seqTips <- unique(inputSeqFam$uphamTips)[!unique(inputSeqFam$uphamTips) %in% 'noTip' & !is.na(unique(inputSeqFam$uphamTips))]
### subset order fasta alignment into each species
for (sp in seqTips){
  path_species_file = paste0(new_folder_species, "/", sp, ".fasta")
  sPaccession <- inputSeqFam %>% 
    filter(uphamTips %in% sp) %>% 
    pull(accession)
  
  subfasta <- big_fasta[names(big_fasta)[which(names(big_fasta) %in% sPaccession)]] 
  write.fasta(subfasta, names = names(subfasta), file.out = path_species_file, open = "w" )
} 

sp5 <- inputSeqFam[!is.na(inputSeqFam$uphamTips),] %>% group_by(uphamTips) %>% tally() %>% filter(n > 4)

print(paste(length(unique(inputSeqFam$uphamTips)), 'species alignments were produced and',dim(sp5)[1], 'species have at least 5 individuals'))
print(paste(table(is.na(inputSeqFam$uphamTips))[2], 'sequences are from species that are not in phylogeny'))

