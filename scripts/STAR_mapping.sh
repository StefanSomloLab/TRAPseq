#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=STAR
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxx.xxx@xxx.xxx


STAR --genomeDir /path-to-STAR-reference/mice --readFilesIn $1 $2 \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
--runThreadN 10 --outFileNamePrefix $3 --readFilesCommand zcat

