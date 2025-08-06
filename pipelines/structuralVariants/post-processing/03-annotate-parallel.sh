#!/bin/bash
# combine with titan calls
# updated to run in parallel

inDir1=intermediate/02/results
inDir2=input/titan/hmm/optimalClusterSolution
outDir1=intermediate/03
outDir2=$outDir1/scripts
outDir3=$outDir1/results
outDir4=$outDir1/logs

# set up directories for cluster output
mkdir -p $outDir1
mkdir -p $outDir2
mkdir -p $outDir3
mkdir -p $outDir4

samples=$(ls -1 $inDir1 | grep vcf | sed 's/-sv-merge.vcf//') # | grep FGC | grep -v clone | grep -v SRR)
for sample in $samples
do
  file=$inDir1/$sample-sv-merge.vcf
    
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
  echo "" >> $script

  echo "sample=$sample" >> $script
  echo "file=$file" >> $script
  echo "inDir1=$inDir1" >> $script
  echo "inDir2=$inDir2" >> $script
  echo "outDir=$outDir3" >> $script
  echo >> $script
  echo "# ** template begin **" >> $script  

  cat 03a-template.sh >> $script

  cmd="sbatch -c 2 $script"
  echo $cmd;
  eval $cmd;
done


