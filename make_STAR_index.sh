#!/bin/bash
#SBATCH --time=0-12
#SBATCH --mem=50GB
#SBATCH --output=make_STAR_index.out
#SBATCH --mail-user=catherineyaozhou@gmail.com
#SBATCH --mail-type=ALL
genome="/oak/stanford/groups/ebutcher/catherine/PLN_EC/files"

STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $genome/GRCm39_STAR_index --genomeFastaFiles $genome/Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile $genome/Mus_musculus.GRCm39.105.gtf
