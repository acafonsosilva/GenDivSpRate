#!/bin/bash
folder_path="/data/biodiv/asilva/mutationRate/"
echo $folder_path

while (( "$#" )); do 
  t=$1
  shift 
done
echo $t ##not clear to me why $@ is not sufficient 

clades=$folder_path"/clades.txt"
cat $clades

for i in {1..25};do

clade=$(sed -n "$i"p $clades)
cd $folder_path"clades/"$clade"/tree"$t
pwd
/users/biodiv/asilva/paml4.9j/src/baseml

done

# To parallelise sets of 5 clades (makes it faster without taking over all available cores 25 * 100)
#
# for i in {1..5};do  ##{10..15}  ## {16..20}  ## {21..25}
# 
# clade=$(sed -n "$i"p $clades)
# cd $folder_path"clades/"$clade"/tree"$t
# pwd
# /users/biodiv/asilva/paml4.9j/src/baseml
# 
# done

# To run MCC using argument of submission file for clade number
#
# clade=$(sed -n "$t"p $clades)
# cd $folder_path"clades/"$clade"/treeMCC"$t
# pwd
# /users/biodiv/asilva/paml4.9j/src/baseml
