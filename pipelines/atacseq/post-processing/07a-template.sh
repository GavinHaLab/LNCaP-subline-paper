
inDir=intermediate/06
inFile2=input/manifest.txt
outDir1=intermediate/07/results
bamCol=9

bam=$(grep $sample_id $inFile2 | cut -f$bamCol)
outFile=$outDir1/$sample_id.txt

cmd="featureCounts -p -a $saf -F SAF -T 1 -o $outFile $bam"
echo $cmd; eval $cmd

