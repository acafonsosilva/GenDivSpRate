library(seqinr)
library(bioseq)
library(tidyverse)
library(ape)

#list to per species alignments
seq_list <- list.files("3.superCrunch_family/Output_alignments",
                       recursive = F, full.names = T)

for (i in seq_list){
  dna <- read_fasta(i)
  
  if(length(dna) > 4){
    species <- tools::file_path_sans_ext(basename(i))
    trans <- seq_translate(dna, code = 2, codon_frame = 1)
    
    if(any(str_detect(trans, "\\*"))){    
      argsFinal <- c('-Xmx2000M',
                     '-jar macse_v2.06.jar',
                     "-prog alignSequences",
                     "-gc_def 2",
                     paste0("-seq ", i, ".fasta"))
      system2('/usr/bin/java', argsFinal) ## this set of sequences were further manually reviewed
    } 
  }
}

### get accession for the final set of sequences

seq_list <- list.files('3.superCrunch_family/Output_alignments/', full.names = TRUE)

seqInfo <- data.frame()
for (i in 1:length(seq_list)){
  
  seqDNA <- read.dna(seq_list[i],format="fasta",as.character = T)
  
  ind <- nrow(seqDNA)
  if(ind > 4){
    sp <- tools::file_path_sans_ext(basename(seq_list[i]))
    dnabioseq <- read_fasta(seq_list[i])
    
    seqInfo <- bind_rows(seqInfo, data.frame(species = sp, 
                                             accession = word(names( dnabioseq), sep = "_",3)))
    
  }
}

seqData <- read_delim('3.superCrunch_family/inputs/genBank_metadata_complete_matchedUpham.txt') %>% 
  filter(accession %in% seqInfo$accession) %>% 
  rename(species = uphamTips)
write.table(seqData, '3.superCrunch_family/GenDiv_accessionIDsPerSpecies.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

  
  