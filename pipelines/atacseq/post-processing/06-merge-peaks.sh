#!/bin/bash
# get union of intersecting peaks
# 20240502 arb - updated to turn off merging

inFile1=input/manifest.txt
inDir1=intermediate/04
outDir=intermediate/06

mkdir -p $outDir

datasets=$(find $inDir1 -type f | sed 's/.*__//' | sed 's/.bed//' | sort | uniq)
						       
for dataset in $datasets
do
  echo $dataset
  outFile=$outDir/$dataset.bed

  find $inDir1 | grep bed | grep $dataset | xargs cat | grep -P '^chr[XY0-9]+\t' | sort | uniq | bedtools sort | bedtools merge > $outFile
#  find $inDir1 | grep bed | grep $dataset | xargs cat | grep -P '^chr[XY0-9]+\t' | sort | uniq | bedtools sort > $outFile  
done
		




