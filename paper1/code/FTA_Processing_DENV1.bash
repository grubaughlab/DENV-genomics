#!/bin/bash

#Unique name given to each library. This will be provided in a list file for each library.  
file_base=$1

#Creates a log for each processing step for each library 
log=${file_base}.FTA.pipeline.log

{

echo "***********************************" 
echo "begin consensus generation for sample: $file_base" 
echo "***********************************" 

#variables set of each read
f1=${file_base}_R1.fastq.gz
f2=${file_base}_R2.fastq.gz

#direct path to indexed reference genomes. (DENV1 Accession JX669465.1, DENV2 Accession KP188569.1)
DENV2_Brazil_Ref=/home/jrf69/project/DENV_FTA_Test/References/KP188569.1.fasta
DENV1_Brazil_Ref=/home/jrf69/project/DENV_FTA_Test/References/JX669465.1.fasta

#Clean reads using Trimmomatic
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 $f1 $f2 ${file_base}_R1_paired.fq.gz ${file_base}_R1_unpaired.fq.gz ${file_base}_R2_paired.fq.gz ${file_base}_R2_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#new variables set of each read
f1=${file_base}_R1_paired.fq.gz
f2=${file_base}_R2_paired.fq.gz

#Map reads to reference (DENV1 here)
bwa mem $DENV1_Brazil_Ref $f1 $f2 > ${file_base}_aln.sam

#.bam to .sample
samtools view -bS ${file_base}_aln.sam > ${file_base}_aln.bam 

#keep only mapped reads
samtools view -b -F 4 ${file_base}_aln.bam > ${file_base}_mapped.bam

#sort final mapped reads 
samtools sort ${file_base}_mapped.bam -o ${file_base}_mapped_sorted.bam 

echo "***********************************" 
echo "finished consensus generation for sample: $file_base" 
echo "***********************************" 

#finish log file for pipeline 
} 2>&1  | tee -a $log

