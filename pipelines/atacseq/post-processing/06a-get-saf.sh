#!/bin/bash
# get saf file

inDir=intermediate/06
outDir=intermediate/06

beds=$(find $inDir | grep bed)
for bed in $beds
do
  echo $bed
  saf=$(echo $bed | sed 's/.bed/.saf/')
  awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' $bed > $saf
  echo $saf
done
