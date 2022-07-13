#!/bin/bash

spath="/oak/stanford/groups/ebutcher/catherine/PLN_EC/input/2843-55622568"
bam="/oak/stanford/groups/ebutcher/scRNAseq/tmp/2018_TU_HEV_LN_HEV_TU_EC/2843-55622568/FASTQ_Generation_2017-12-05_04_34_16Z-65510628"

cd $bam
for cell in $(ls | cut -f1 -d'-' | cut -f2 -d'_' | uniq | sort)
do      
        for (( i = 1 ; i <= 4; i++ ))
        do
                d=($(ls -d *"_"$cell"_L00"$i*))
                f=($(cd $d; ls *".fastq.gz"))
                echo -n $bam/$d/$f, >> $spath/$cell.txt
        done
        truncate -s -1 $spath/$cell".txt" 
done
