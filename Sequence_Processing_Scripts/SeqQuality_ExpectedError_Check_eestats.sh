#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=1 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=10G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=30G # < if you are threading, ask for whole memory and not memory per CPU
#SBATCH --time=1-00:00:00     # 1 day
#SBATCH --output=EEstats_Seq_Qual_Check.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="EEstats_Seq_Qual_Check.21"
#SBATCH -p batch

module load usearch/11

if [[ ! -d ./EEstats_Results ]]; then
    mkdir EEstats_Results
fi

if [[ ! -d ./EEstats_Results/16S_EEstats ]]; then
    mkdir EEstats_Results/16S_EEstats
fi

if [[ ! -d ./EEstats_Results/ITS2_EEstats ]]; then
    mkdir EEstats_Results/ITS2_EEstats
fi

## * be sure to gunzip files before running script

for FILE in 16S_Seqs/*.fastq;
do
    f=$(basename $FILE)
    SAMPLE=${f%.fastq*}
    
    #usearch -fastq_eestats $FILE -output ${SAMPLE}_eestats.txt
    usearch -fastq_eestats2 $FILE -output ${SAMPLE}_eestats2.txt
    
    mv ${SAMPLE}_eestats2.txt EEstats_Results/16S_EEstats
done

for FILE in ITS2_Seqs/*.fastq;
do
    f=$(basename $FILE)
    SAMPLE=${f%.fastq*}
    
    #usearch -fastq_eestats $FILE -output ${SAMPLE}_eestats.txt
    usearch -fastq_eestats2 $FILE -output ${SAMPLE}_eestats2.txt
    
    mv ${SAMPLE}_eestats2.txt EEstats_Results/ITS2_EEstats
    
done

# More info here: https://www.drive5.com/usearch/manual/cmd_fastq_eestats.html
# https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html

