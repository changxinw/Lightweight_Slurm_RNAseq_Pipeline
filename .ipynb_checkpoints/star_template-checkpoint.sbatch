#!/bin/bash
#SBATCH -n 4 #Number of cores
#SBATCH -N 1
#SBATCH -t 3-00:00 # Runtime in minutes
#SBATCH --mem=100G #Memory per node in MB (see also --mem-per-cpu)
#SBATCH --job-name={{SAMPLE_ID}}
#SBATCH -o ./log/{{SAMPLE_ID}}.log
#SBATCH -e ./log/{{SAMPLE_ID}}.err

fastp --thread 4 --detect_adapter_for_pe --in1 {{FASTQ1}} --in2 {{FASTQ2}} --out1 ./fastp/{{SAMPLE_ID}}_R1_001.fastq.gz --out2 ./{{SAMPLE_ID}}_R2_001.fastq.gz -h ./fastp/{{SAMPLE_ID}}_fastp.html
{{STAR_PATH}} \
--genomeDir {{STAR_INDEX}} \
--runThreadN 4 \
--readFilesIn ./fastp/{{SAMPLE_ID}}_R1_001.fastq.gz ./fastp/{{SAMPLE_ID}}_R2_001.fastq.gz \
--outFileNamePrefix ./star/{{SAMPLE_ID}}_ \
--readFilesCommand zcat \
--genomeLoad NoSharedMemory \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMheaderHD @HD VN:1.4 SO:coordinate \
--outSAMunmapped Within \
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--sjdbScore 1