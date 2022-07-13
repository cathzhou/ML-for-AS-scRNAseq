#!/bin/bash

#SBATCH --time=0-12
#SBATCH --array=0-383%10
#SBATCH --nodes=2
#SBATCH --mem=100GB
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --mail-user=catherineyaozhou@gmail.com
#SBATCH --mail-type=ALL

bam="/oak/stanford/groups/ebutcher/scRNAseq/tmp/2018_TU_HEV_LN_HEV_TU_EC/2843-55622568/FASTQ_Generation_2017-12-05_04_34_16Z-65510628"
in="/oak/stanford/groups/ebutcher/catherine/PLN_EC/rmats/input/2843-55622568"
out="/oak/stanford/groups/ebutcher/catherine/PLN_EC/star/output/2843-55622568"
STARindex="/oak/stanford/groups/ebutcher/catherine/PLN_EC/files/GRCm39_STAR_index"
gtf="/oak/stanford/groups/ebutcher/catherine/PLN_EC/files/Mus_musculus.GRCm39.105.gtf"

cd $bam

ID=($(ls | cut -f1 -d'-' | cut -f2 -d'_' | uniq | sort -n))

slurm_id_new=$SLURM_ARRAY_TASK_ID

STAR --runThreadN 8 --outSAMtype BAM SortedByCoordinate --genomeDir $STARindex --sjdbGTFfile $gtf --outFileNamePrefix $out/${ID[$SLURM_ARRAY_TASK_ID]} --readFilesIn $(cat $in/${ID[$SLURM_ARRAY_TASK_ID]}.txt) --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts

