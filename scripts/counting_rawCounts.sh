#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=count
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxx.xxx@xxx.xxx


featureCounts -T 10 -p -t exon -g gene_id \
		-a  /path-to/gencode.vM25.annotation.gtf \
		-o featureCounts.gene_level_exon_counts_PE_samples \
		-J *.bam
