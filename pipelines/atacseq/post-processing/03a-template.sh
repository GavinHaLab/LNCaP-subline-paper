
inFile1=intermediate/02.txt
inFile2=input/manifest.txt
outDir1=intermediate/03/results
bedCol=8

sample_ids=$(grep $sample $inFile1 | cut -f1 | grep -v $sample_id1)

peak_bed1=$(grep $sample_id1 $inFile2 | cut -f$bedCol)

for sample_id2 in $sample_ids
do
  peak_bed2=$(grep $sample_id2 $inFile2 | cut -f$bedCol)    
  outFile1=$outDir1/$sample_id1--$sample_id2.bed

  # require at least 50% overlap
  cmd="bedtools intersect -f 0.90 -r -a $peak_bed1 -b $peak_bed2 > $outFile1"
  echo $cmd;
  eval $cmd
done

