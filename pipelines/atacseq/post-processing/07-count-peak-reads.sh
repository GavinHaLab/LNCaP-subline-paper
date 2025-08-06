#!/bin/bash
# get peaks that intersect with other peaks in other samples

inDir=intermediate/06
inFile1=intermediate/02.txt
template_sh=07a-template.sh
outDir1=intermediate/07
outDir2=$outDir1/scripts
outDir3=$outDir1/results
outDir4=$outDir1/logs

mkdir -p $outDir1
mkdir -p $outDir2
mkdir -p $outDir3
mkdir -p $outDir4

sample_ids=$(tail -n+2 $inFile1 | cut -f1 | sort | uniq)
for sample_id in $sample_ids
do
  echo $sample_id
  dataset=$(echo $sample_id | sed 's/__[0-9]$//' | sed 's/.*__//')
  saf="$inDir/$dataset.saf"

  script=$outDir2/$sample_id.sh
  std_out=$outDir4/$sample_id.out
  std_err=$outDir4/$sample_id.err

  echo "#!/bin/bash" > $script
  echo "#SBATCH --job-name=$sample_id" >> $script
  echo "#SBATCH --output=$std_out" >> $script
  echo "#SBATCH --error=$std_err" >> $script
  echo "#SBATCH --time=01:00:00" >> $script
  echo "#SBATCH --partition=campus-new" >> $script
  echo "#SBATCH --nodes=1" >> $script
  echo "" >> $script

  echo "sample_id=$sample_id" >> $script
  echo "dataset=$dataset" >> $script
  echo "saf=$saf" >> $script    
  cat $template_sh >> $script
  
  cmd="sbatch -c 1 $script"
  echo $cmd; eval $cmd;
done

