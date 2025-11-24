#!/bin/bash
#SBATCH -n 8 #Number of cores
#SBATCH -N 1
#SBATCH -t 1-00:00 # Runtime in minutes
#SBATCH --mem=100G #Memory per node in MB (see also --mem-per-cpu)
#SBATCH --job-name=featurecounts
#SBATCH -o ./log/featurecounts.log
#SBATCH -e ./log/featurecounts.err

{{FEATURECOUNTS_PATH}} -T 8 -p -t exon -g gene_name -B -C \
-a {{GTF}} \
-o ./feature_count/gene_counts.txt \
{{BAMS}}
