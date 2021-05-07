#!/bin/bash

# Instructions:
#	1. Create a folder run/ containing folders whose names are descriptive of each sample. 
#   2. Each sub-folder must contain one/two (single/paired end) fq.gz.
#	3. Resulting files will appear in each of these sub-folders.

#SBATCH -J STAR
#SBATCH -N 1
threads=6
#SBATCH -n $threads 		

## Change manually
indexDirectory=/home/areyes/notFur
gtf_path=/home/ska/jtena/genomes/Nfu_20150522.genes_20140630.gtf

## Don't change (except cores)
SCRATCH=/scratch/areyes_$SLURM_JOB_ID
mkdir -p $SCRATCH || exit $?
mkdir $SCRATCH/starIndex
mkdir $SCRATCH/results
wd=$PWD

folder=$1
namePrefix=${folder##*/}

cp $indexDirectory/* $SCRATCH/starIndex/

file1=$folder/*.fq.gz
#file2=$folder/*_2.fq.gz

cp $file1 $SCRATCH/1.fq.gz
#cp $file2 $SCRATCH/2.fq.gz
cp $gtf_path $SCRATCH/model.gtf

cd $SCRATCH
gunzip *.fq.gz

#Add if this option if no strand specific and wanna use cufflinks
#--outSAMstrandField intronMotif
STAR  --genomeDir ./starIndex --outSAMtype BAM SortedByCoordinate --sjdbGTFfile model.gtf --runThreadN 10 --readFilesIn 1.fq 2.fq --outFileNamePrefix results/$namePrefix --quantMode GeneCounts

dest=${wd}/${folder}
mv results/* $folder

rm -rf $SCRATCH
