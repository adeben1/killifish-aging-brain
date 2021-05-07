# Merge IDR files for differential peak analysis

# Generate merge file from IDRs peaks

tissue=MB

echo "Tissue --> $tissue"
echo "IDR files to be merged...
"

for folder in ./*${tissue};
do
	cp ${folder}/TRACKS/*idrConsPeaks.bed ./;
	ls ${folder}/TRACKS/*idrConsPeaks.bed;
done

cat *idrConsPeaks.bed | cut -f1,2,3 | sort -k1,1 -k2,2n > allPeaks1_${tissue}.bed
mergeBed -i allPeaks1_${tissue}.bed > allPeaks_${tissue}.bed

echo "IDR files merged"

rm *idrConsPeaks.bed
rm allPeaks1_${tissue}.bed

# Process nucfree.bed files and intersect with merged file

echo "Nucfree files to be intersected with merge file...
"

for file in ./*${tissue}/*;
do
	if [[ $file == *nucfree.bed ]];
	then 
		echo $file
		cut -f1,2,3 $file | sort -k1,1 -k2n > ${file}.sorted
		intersectBed -a allPeaks_${tissue}.bed \
			     -b ${file}.sorted -c \
			     > ${file}.count
		rm ${file}.sorted
	fi; 
done

echo "Count files generated"

### Next code lines were added after running the script (run them line by line)

# Generate 2 columns count files

# awk -v OFS='\t' '{$1=$1":"$2"-"$3; $2=$4; $3=""; $4=""}1' 4columns.counts > 2columns.count
# use it with a for loop that select *count files

# Make count matrices
# use sth similar to:

# paste 8sem_FB/N-fur_ATAC_8sem_FB1_nucfree.bed.counts2 
#	<(cut -f2 8sem_FB/N-fur_ATAC_8sem_FB2_nucfree.bed.counts2) 
#	<(cut -f2 16sem_FB/N-fur_ATAC_16sem_FB1_nucfree.bed.counts2) 
#	<(cut -f2 16sem_FB/N-fur_ATAC_16sem_FB2_nucfree.bed.counts2) 
#	<(cut -f2 24sem_FB/N-fur_ATAC_24sem_FB1_nucfree.bed.counts2) 
#	<(cut -f2 24sem_FB/N-fur_ATAC_24sem_FB2_nucfree.bed.counts2) 
#	> atacpeaks_counts_FB.txt
