#!/bin/bash
# for each sample merge sv calls
# source from original manta sv vcf
# 20240209 arb - updated to run in parallel

svMerge=sv-merge/sv-merge.py
inFile1=intermediate/01a.txt
outDir1=intermediate/02
outDir2=$outDir1/scripts
outDir3=$outDir1/results
outDir4=$outDir1/logs

# set up directories for cluster output
mkdir -p $outDir1
mkdir -p $outDir2
mkdir -p $outDir3
mkdir -p $outDir4

samples=$(tail -n+2 $inFile1 | cut -f1 | sort | uniq) # | grep APIPC) 

for sample in $samples
do
  echo $sample

  searchString="$sample/"
  vcfs=$(awk -v sample=$sample '$1 == sample { print $2 }' $inFile1)

  script=$outDir2/$sample.sh
  std_out=$outDir4/$sample.out
  std_err=$outDir4/$sample.err

  echo "#!/bin/bash" > $script
  echo "#SBATCH --job-name=$sample" >> $script
  echo "#SBATCH --output=$std_out" >> $script
  echo "#SBATCH --error=$std_err" >> $script
  echo "#SBATCH --time=01:00:00" >> $script
  echo "#SBATCH --partition=campus-new" >> $script
  echo "#SBATCH --nodes=1" >> $script

  echo >> $script
  echo >> $script

  echo "time $svMerge -v $vcfs -s $sample -o $outDir3 -ro 0.8 --verbose T" >> $script

  cmd="sbatch -c 2 $script"
#  cmd="bash $script"
  echo $cmd;
  eval $cmd;

done
