#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=1 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=10G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=10G # < if you are threading, ask for whole memory and not memory per CPU
#SBATCH --time=1-00:00:00     # 1 day
#SBATCH --output=FastQC_QualCheck.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="FastQC_QualCheck"
#SBATCH -p aronsonlab

module load fastqc/0.11.9

if [[ ! -d ./FastQC_Results ]]; then
    mkdir FastQC_Results
fi

if [[ ! -d ./FastQC_Results/16S_FastQC ]]; then
    mkdir FastQC_Results/16S_FastQC
fi

if [[ ! -d ./FastQC_Results/ITS2_FastQC ]]; then
    mkdir FastQC_Results/ITS2_FastQC
fi

for FILE in 16S_Seqs/*.fastq.gz;
do
    f=$(basename $FILE)
    SAMPLE=${f%.fastq*}
    
    fastqc $FILE --outdir=./FastQC_Results/16S_FastQC
    
done

for FILE in ITS2_Seqs/*.fastq.gz;
do
    f=$(basename $FILE)
    SAMPLE=${f%.fastq*}
    fastqc $FILE --outdir=./FastQC_Results/ITS2_FastQC
    
done
