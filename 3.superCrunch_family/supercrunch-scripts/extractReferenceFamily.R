#### script to extract the closest sequence to the expected size (~1140) to become reference within the group
### Ana C. Afonso Silva script and not a SuperCrunch file

suppressWarnings(library(seqinr))

args <- commandArgs()
clade <- args[6]
print(paste("copying reference for",clade))

folder_family_path <- args[7]

big_fasta_path <- paste0(folder_family_path, "/1_Taxa_Assessment/Matched_Taxa.fasta")

big_fasta <- read.fasta(big_fasta_path, seqtype = c("DNA"),
                        seqonly = FALSE, forceDNAtolower = FALSE, set.attributes = TRUE)

ref <- big_fasta[names(big_fasta) %in% names(which.min(abs(lengths(big_fasta) - 1140)))]

if(length(big_fasta) > 1) {
write.fasta(ref, names = paste0('reference_',names(ref)), file.out = paste0(folder_family_path,"/reference.fasta"), open = "w")
  } else {
    file.copy(paste0(folder_path, "/reference.fasta"), paste0(folder_family_path,"/reference.fasta"), overwrite = TRUE)
}