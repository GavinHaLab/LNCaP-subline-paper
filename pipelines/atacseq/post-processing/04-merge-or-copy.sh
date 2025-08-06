#!/bin/bash


inFile1=intermediate/02.txt
inFile2=input/manifest.txt
inDir1=intermediate/03/results
outDir=intermediate/04
peakCol=8

mkdir -p $outDir

samples=$(tail -n+2 $inFile1 | cut -f2 | sort | uniq)

for sample in $samples
do
    echo $sample
    num=$(grep $sample $inFile1 | wc -l)
    outFile=$outDir/$sample.bed

    # figure out what to search for
    prefix=$(echo $sample | sed 's/__.*//')
    suffix=$(echo $sample | sed 's/.*__//')
    search="__$suffix\__[0-9]"
    search="${prefix}__${suffix}__[0-9]"
    
    if [[ $num == 1 ]]; then
#	echo 'no replicates'
	inFile=$(grep $search $inFile2 | cut -f$peakCol)
	cat $inFile | grep -P '^chr[XY0-9]+\t' | cut -f1-3 | bedtools sort > $outFile
    else
#      echo 'replicates'
      find $inDir1 | grep bed | grep $search | xargs cat | grep -P '^chr[XY0-9]+\t' | bedtools sort | bedtools merge > $outFile		
    fi
done
		




