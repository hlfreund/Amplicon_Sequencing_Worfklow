#!/bin/bash -l

for name in *.fa
do
    newname="${name#*EA_Pool_}"
    echo $newname
    mv $name $newname
done

for filename in *.fastq.gz
do
    mv "$filename" "${filename//_001/}"

done

for filename in *.fastq
do
    mv "$filename" "${filename/_S*_/_}"

done

for filename in *.fastq.gz
do
    mv "$filename" "${filename/_S[0-9]_/_}"

done


for i in *_R1_001.fastq;
do
    file=$(basename $i)
    SAMPLE=${file%_S*}
    echo ${SAMPLE} ${i} "Checking quality of sequences with FastQC" >> New_and_Old_SampleNames.txt
    
done

for dir in *metabat-bins1500-*
do
    file=$(basename $dir)
    SAMPLE=${file%_contigs*}
    mv $dir ${SAMPLE}_bins

done


for filename in *.fastq
do
    mv "$filename" "${filename/_S*_/_}"
    echo ${filename} ${filename/_S*_/_} "Checking quality of sequences with FastQC" >> New_and_Old_SampleNames.txt

done
