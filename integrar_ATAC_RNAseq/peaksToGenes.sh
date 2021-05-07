##### USAGE: ./peaksToGenes.sh folder_to_bedfiles (don't use the dot path {./} to specify folder)

# Activate environment where closestBed is located

conda activate master 

tss=/home/areyes/nfur/genome/Nfu_20150522.genes_20150922.TSS.sorted.bed #tss file sorted (chr/start/end/strand)

for f in $1*bed
do

	name=$( echo $f | cut -d '.' -f 1 )
	output1=$name.upstream.txt
        output2=$name.downstream.txt

	# Change coordenates format (from "chr:start-end" to "chr\tstart\tend")
	
	sed -i 's/[:-]/\t/g' $f	
	
	# Obtain the closest gene to each peak
	
	echo "$f	Processing..."
	closestBed -id -D ref -k 1 -a $f -b $tss > $output1
	closestBed -iu -D ref -k 1 -a $f -b $tss > $output2
	
	# Extract gene name column and give Danre homologous
	
	output3=$name.killigenes.txt
	cat $output1 $output2 | cut -f8 | grep -v "\." > $output3
	Rscript --vanilla /home/areyes/nfur/ATACseq/scripts/mapNames.R $output3 $name

done
