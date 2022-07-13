#!/bin/bash

#SBATCH --array=0-383:10%10
#SBATCH --nodes=4
#SBATCH --mem=50GB
#SBATCH --time=0-24
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --mail-user=catherineyaozhou@gmail.com
#SBATCH --mail-type=ALL

out="/oak/stanford/groups/ebutcher/catherine/PLN_EC/output/2965-81446365"
gtf="/oak/stanford/groups/ebutcher/catherine/PLN_EC/files/Mus_musculus.GRCm39.105.gtf"
star="/oak/stanford/groups/ebutcher/catherine/PLN_EC/files/GRCm39_STAR_index"
bam="/oak/stanford/groups/ebutcher/scRNAseq/tmp/2018_TU_HEV_LN_HEV_TU_EC/2965-81446365/FASTQ_Generation_2018-06-13_02_24_30Z-103228135"
inputs="/oak/stanford/groups/ebutcher/catherine/PLN_EC/input/2965-81446365"
tmp="/oak/stanford/groups/ebutcher/catherine/PLN_EC/tmp/2965-81446365"

cd $bam

ID=($(ls | cut -f1 -d'-' | cut -f2 -d'_' | uniq | sort -n))

if [[ $SLURM_ARRAY_TASK_ID -lt 380 ]]
then
        for ((i=0;i<=9;i++)); do
                rmats.py --s1 $inputs/${ID[$SLURM_ARRAY_TASK_ID + $i]}.txt --gtf $gtf --bi $star -t single --readLength 76 --variable-read-length --nthread 4 --od $out/${ID[$SLURM_ARRAY_TASK_ID + $i]} --tmp $tmp/tmp_${ID[$SLURM_ARRAY_TASK_ID + $i]} --statoff
        done
else
        for ((i=0;i<=3;i++)); do
                rmats.py --s1 $inputs/${ID[$SLURM_ARRAY_TASK_ID + $i]}.txt --gtf $gtf --bi $star -t single --readLength 76 --variable-read-length --nthread 4 --od $out/${ID[$SLURM_ARRAY_TASK_ID + $i]} --tmp $tmp/tmp_${ID[$SLURM_ARRAY_TASK_ID + $i]} --statoff
        done
fi
