#!/bin/bash

#SBATCH -J htseq-count
#SBATCH -N 1
threads=6
#SBATCH -n $threads 		

## Change manually
gtf_path=/home/ska/jtena/genomes/Nfu_20150522.genes_20140630.gtf

## Don't change (except cores)
SCRATCH=/scratch/areyes_$SLURM_JOB_ID
mkdir -p $SCRATCH || exit $?
mkdir $SCRATCH/results
wd=$PWD

folder=$1
namePrefix=${folder##*/}

file1=$folder/*.bam

cp $file1 $SCRATCH/alignment_file.bam
cp $gtf_path $SCRATCH/model.gtf

cd $SCRATCH

htseq-count -f bam --nonunique all -s no alignment_file.bam model.gtf > results/$namePrefix

mv results/* $folder

rm -rf $SCRATCH
