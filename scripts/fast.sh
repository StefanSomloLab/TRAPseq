#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=fastp
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxx.xxx@xxx.xxx


fastp -i $1 \
	  -I $2 \
	  -o $1".clean.fastq"\
	  -O $2".clean.fastq"\
	  -j $3".json" \
	  -h $3".html" \
	  -w 8

gzip $1".clean.fastq"
gzip $2".clean.fastq"

	  

