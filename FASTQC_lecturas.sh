#!/bin/bash

for folder in ./*; 
do
	cd $folder; file=./*.bam; 
	fastqc $file --outdir=/home/areyes/nfur_rnaseq/quality_control/ ;
	cd ..; 
done 
