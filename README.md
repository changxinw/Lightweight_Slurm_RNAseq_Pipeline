# Light-weight RNA-seq Data Processing Pipeline for Slurm-based HPC

## Steps
### 1. Copy the folder and rename it as your project name
`cp -r /hpc/group/mchenlab/Changxin/RNAseqPipeline /hpc/group/mchenlab/Changxin/gene_KO_RNAseq`
### 2. Change the working directory to your project folder
`cd /hpc/group/mchenlab/Changxin/gene_KO_RNAseq`
### 3. Edit the `config.yaml` file, specify the software path, index path, genome, and fastq files
One trick here is that you can make a soft link of your file path
`ln -s /hpc/group/mchenlab/Changxin/geneKORNAseq/25257-02-11142025_151637 /hpc/group/mchenlab/Changxin/gene_KO_RNAseq/fastq`
Then you can edit the `config.yaml` file to specify the fastq files and genome.
### 4. Run the `make_sbatch` to generate sbtach files and submit
`Rscript make_sbatch --submit`
### 5. Monitor the job status
`squeue -u {netid}`
You should be able to find your count matrix at `feature_count` folder.
## 6. Perform differential expression analysis
This step requires some manual efforts. I have provided `rna.ipynb`, which includes code for differential expression analysis using DESeq2. You can run the notebook in https://dcc-ondemand-01.oit.duke.edu/pun/sys/dashboard when you are a member of Duke.
