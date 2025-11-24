# Light-weight RNA-seq Data Processing Pipeline for Chen Lab at DCC

## Steps
### 1. Copy the folder and rename it as your project name
`cp -r /hpc/group/mchenlab/Changxin/RNAseqPipeline /hpc/group/mchenlab/Changxin/CDK12_KO_RNAseq_V2`
### 2. Change the working directory to your project folder
`cd /hpc/group/mchenlab/Changxin/CDK12_KO_RNAseq_V2`
### 3. Edi the `config.yaml` file, specify the genome and fastq files
One trick here is that you can make a soft link of your file path
`ln -s /hpc/group/mchenlab/weiling/CDK12\ RNAseq/25257-02-11142025_151637 /hpc/group/mchenlab/Changxin/CDK12_KO_RNAseq_V2/fastq`
Then you can edit the `config.yaml` file to specify the fastq files and genome.

### 4. Run the `make_sbatch` to generate sbtach files and submit
`Rscript make_sbatch --submit`
### 5. Monitor the job status
`squeue -u {netid}`
You should be able to find your count matrix at `feature_count` folder.
## 6. Perform differential expression analysis
This step requires some manual efforts. I have provided `rna.ipynb`, which includes code for differential expression analysis using DESeq2. You can run the notebook in https://dcc-ondemand-01.oit.duke.edu/pun/sys/dashboard.
