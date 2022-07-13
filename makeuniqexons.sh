#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --time=0-2
#SBATCH --output=makeuniqexons.out
#SBATCH --mail-user=catherineyaozhou@gmail.com
#SBATCH --mail-type=ALL

out1="/oak/stanford/groups/ebutcher/catherine/PLN_EC/rmats/output/merged_output/2843-55622568"
out2="/oak/stanford/groups/ebutcher/catherine/PLN_EC/rmats/output/merged_output/2965-81446365"
spath="/oak/stanford/groups/ebutcher/catherine/PLN_EC/csv/merged_csv"

cd $out1

for cell in $(ls | sort -n)
do
        awk '{print $2, $3, $4, $5, $6, $7}' $out1/$cell/"SE.MATS.JC.txt" | sed '1d' | sed 's/\"//g' | cut -c 20- | sed -e 's/ /_/g' >> $spath/"allexons.txt"
done

cd $out2

for cell in $(ls | sort -n)
do
        awk '{print $2, $3, $4, $5, $6, $7}' $out2/$cell/"SE.MATS.JC.txt" | sed '1d' | sed 's/\"//g' | cut -c 20- | sed -e 's/ /_/g' >> $spath/"allexons.txt"
done

cat $spath/"allexons.txt" | sort | uniq >> $spath/"uniqexons.txt"
