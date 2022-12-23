#!/bin/bash
##### INPUT #####
# Prepared in Folder_per_family.R
# 1) cleanSeqLabels.txt - species list for all retrieved sequences but excluding names with problems
# 2) For each family - A folder with family.fasta

# Prepared with nameMatching scripts
## 3) genBank_metadata_complete_matchedUpham.txt - list matching sequence data to used Upham Phylogeny
## This is needed for the customized R script Fasta_file_per_species.R and both pipelines take time to run
## Instead of waiting for the nameMatching scripts to finish, it's easier to run them while running steps 1 to 5 of this script

# Enter the full path of this folder 
supercrunch_folder_path="/data/biodiv/asilva/supercrunch_family/"

# Enter the full path of this folder
supercrunch_scripts_path="/data/biodiv/asilva/supercrunch_family/supercrunch-scripts"

# extract family from submission file  when parallelizing analysis
clade=$(sed -n $@p $supercrunch_folder_path"/Taxon_families.txt")
echo $@p
echo $clade

family_folder=$supercrunch_folder_path/families/$clade

cd $family_folder

taxa_file=$supercrunch_folder_path/"cleanSeqLabels.txt"


### 1 Taxa Assessment - Taxa_Assessment.py
             
mkdir $family_folder/1_Taxa_Assessment

python $supercrunch_scripts_path/Taxa_Assessment.py -i $family_folder"/"$clade"_cytb.fasta" -t $taxa_file -o $family_folder/1_Taxa_Assessment


### 2 Blast Extract - Reference_Blast_Extract.py

/usr/bin/Rscript $supercrunch_scripts_path/extractReferencefamily.R  $clade $supercrunch_folder_path ### to use a sequence from within each family as outgroup

mkdir $family_folder/2_Blast_Extract

cp $family_folder/reference.fasta $family_folder/2_Blast_Extract
cp $family_folder/1_Taxa_Assessment/Matched_Taxa.fasta $family_folder/2_Blast_Extract

python $supercrunch_scripts_path/Reference_Blast_Extract.py -i $family_folder/2_Blast_Extract -d reference.fasta -e Matched_Taxa.fasta -b dc-megablast -m span -o $family_folder/3_Blast_Extract --threads 4 


### 3 Adjust Direction - Adjust_Direction.py

mkdir $family_folder/3_Adjust_Direction

cp $family_folder/2_Blast_Extract/Filtered-Fasta-Files/Matched_Taxa_extracted.fasta $family_folder/3_Adjust_Direction

python $supercrunch_scripts_path/Adjust_Direction.py -i $family_folder/3_Adjust_Direction --threads 4 -o $family_folder/3_Adjust_Direction


### 4 Coding Translation Checks - Coding_Translation_Tests.py

mkdir $family_folder/4_Coding_Translation_Checks

cp $family_folder/3_Adjust_Direction/Adjusted-Fasta-Files/*.fasta $family_folder/4_Coding_Translation_Checks

python $supercrunch_scripts_path/Coding_Translation_Tests.py -i $family_folder/4_Coding_Translation_Checks --table vertmtdna -o $family_folder/4_Coding_Translation_Checks


## 5 Align - Align.py - use muscle but script was edited to include --leavegappyregion with --accurate flag

mkdir $family_folder/5_Align

cp $family_folder/4_Coding_Translation_Checks/Translation-Passed-Seqs/*.fasta $family_folder/5_Align

python $supercrunch_scripts_path/Align.py -i $family_folder/5_Align -a mafft --threads 4 -o $family_folder/5_Align --accurate 


## 6 Per Species fasta files and copy to genetic diversity folder = not a SuperCrunch script

/usr/bin/Rscript $supercrunch_scripts_path/Fasta_file_per_species.R  $clade $supercrunch_folder_path


## 7 Removes any site with just gaps - mainly an issue in bigger alignments

python $supercrunch_scripts_path/Trim_Alignments_Trimal.py -i $family_folder/6_Per_species_fasta_files -a noallgaps -f fasta  -o $family_folder

cd $family_folder/Trimmed-fasta-Files/

for file in *; do mv "$family_folder/Trimmed-fasta-Files/$file" "$family_folder/Trimmed-fasta-Files/$(basename "$file" _trimmed.fasta).fasta"; done 
 
cp -r $family_folder/Trimmed-fasta-Files/* $supercrunch_folder_path/Output_alignments/

