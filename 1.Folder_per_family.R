## Script to retrieve Genbank sequence data per family and setup folders to parallelize SuperCrunch
library(tidyverse)
library(reutils)
library(taxize)
library(seqinr)

rentrez::set_entrez_key("78915277c9890b626a8b082c889ceca86a08")
key = '78915277c9890b626a8b082c889ceca86a08'

#### DATA input####
folder_path <- "/data/biodiv/asilva/supercrunch_family/"

Taxon_family <- downstream('mammalia', db = 'ncbi', 
                           downto = 'family', intermediate = FALSE)[['mammalia']] %>%
  arrange(childtaxa_name) %>%
  pull(childtaxa_name)

# Creating a folder per family
dir.create(paste0(folder_path, '/families'))


Taxon_family <- data.frame(families = Taxon_family, Nseq = NA)
for (family in Taxon_family$families) {
  print(family)
  new_folder <- paste0(folder_path, "/families/", family)
  
  dir.create(new_folder)
  
  #### get sequence from ncbi per family
  taxa_CYTB <- paste0(family,"[Organism] AND CYTB NOT Homo sapiens[Organism]")
  
  temp <- esearch(term = taxa_CYTB, db = 'nuccore', usehistory = TRUE)
  if (is.null(getError(temp)$wrnmsg)) {
  temp
  accessions <- efetch(temp, rettype = "fasta", retmode = "text",
                       outfile = paste0(new_folder,"/", family,"_cytb.fasta"))
  
  seqDNA <- read.fasta(paste0(new_folder,"/", family,"_cytb.fasta"), 
                               seqtype = c("DNA"),seqonly = FALSE, forceDNAtolower = FALSE, 
                               set.attributes = TRUE, whole.header = TRUE)
  
  data <- data.frame(accession = word(names(seqDNA)),
             genBank_seqName = sub(".*? ", "", names(seqDNA)),
             family = family,
             stringsAsFactors = FALSE)

  temp <- data.frame()
  for (i in data$accession){
    temp <- bind_rows(temp,data.frame(accession = i, suppressWarnings(
      ncbi_get_taxon_summary(genbank2uid(i)[[1]][1])),stringsAsFactors = FALSE))
  }
  
  data <- left_join(temp,data, by = 'accession') %>% 
    rename(seqSp = name)
  inputSeq <- bind_rows(inputSeq, data)
  
  if(length(data$seqSp) == length(seqDNA)){
    names(seqDNA) <- paste(data$accession, data$seqSp)
    write.fasta(seqDNA, names = names(seqDNA), 
                file.out = paste0(new_folder,"/", family,"_cytb.fasta"), open = "w" )
  }

  Taxon_family[which(Taxon_family$families %in% family),'Nseq'] <- length(seqDNA)
  }
}

saveRDS(inputSeq, paste0(folder_path, "/genBank_metadata_complete.rds"))

seqCleanSp <- inputSeq %>%
  filter(
    !str_detect(seqSp, "\\."), #exclude sp. and abreviated genus
         !str_detect(seqSp, " x "), #exclude hybrids
         !str_detect(seqSp, "\\?"), 
         !str_detect(seqSp, "-"),
         !str_detect(seqSp, "\\("),
         !str_detect(seqSp, "\\["),
         !str_detect(seqSp, "\\<i>"),
         !str_detect(seqSp, "\\'")) %>%
  arrange(seqSp) %>% pull(seqSp) %>% word(.,1,2) %>% unique()

write.table(seqCleanSp,paste0(folder_path, "/cleanSeqLabels.txt"), quote = F, col.names = F, row.names = F, sep = '\t')

write.table(Taxon_family$families,  file = paste0(folder_path, "/Taxon_families.txt"), 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(Taxon_family,  file = paste0(folder_path, "/Taxon_families_Nseq.txt"), 
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

dir.create(paste0(folder_path, '/Output_alignments'))