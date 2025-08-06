#!/bin/bash
# get peaks that intersect with other peaks in other samples

inFile1=intermediate/02.txt
inFile2=input/manifest.txt
template_sh=03a-template.sh
outDir1=intermediate/03
outDir2=$outDir1/scripts
outDir3=$outDir1/results
outDir4=$outDir1/logs

mkdir -p $outDir1
mkdir -p $outDir2
mkdir -p $outDir3
mkdir -p $outDir4

sample_ids=$(tail -n+2 $inFile1 | cut -f1 | sort | uniq)
for sample_id1 in $sample_ids
do
  echo $sample_id1

  sample=$(grep $sample_id1 $inFile1 | cut -f2 | sort | uniq)
  
  script=$outDir2/$sample_id1.sh
  std_out=$outDir4/$sample_id1.out
  std_err=$outDir4/$sample_id1.err

  echo "#!/bin/bash" > $script
  echo "#SBATCH --job-name=$sample_id1" >> $script
  echo "#SBATCH --output=$std_out" >> $script
  echo "#SBATCH --error=$std_err" >> $script
  echo "#SBATCH --time=01:00:00" >> $script
  echo "#SBATCH --partition=campus-new" >> $script
  echo "#SBATCH --nodes=1" >> $script
  echo "" >> $script

  echo "sample=$sample" >> $script
  echo "sample_id1=$sample_id1" >> $script

  cat $template_sh >> $script

  cmd="sbatch -c 2 $script"
  echo $cmd;
  eval $cmd;
done

