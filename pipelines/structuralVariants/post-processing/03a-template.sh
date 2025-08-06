script=titan-code/combineSVABAandTITAN.R
utils=titan-code/svaba_utils.R
genomeBuild=hg38
genomeStyle=UCSC
minSPAN=1
minInversionSPAN=1000
script2=03b-annotate-genes.R

echo $sample
vcf=$file

# check if lncap_fcg then add _ to make unique
label=$sample
if [ $sample == "LNCaP_FGC" ]; then
    label=$sample\_\_recalibrated
fi

bin=$(ls $inDir2/$sample\_recalibrated\_cluster?.titan.ichor.cna.txt 2> /dev/null)
bin=$(ls $inDir2/$label\_*cluster?.titan.ichor.cna.txt 2> /dev/null)
seg=$(ls $inDir2/$sample\_recalibrated\_cluster?.titan.ichor.seg.noSNPs.txt 2> /dev/null)
seg=$(ls $inDir2/$label\_*cluster?.titan.ichor.seg.noSNPs.txt 2> /dev/null)
outDir=$outDir
outputSVFile=$outDir/$sample.svabaTitan.sv.txt
outputSVFile2=$outDir/$sample.svabaTitan.genes.sv.txt  
outputCNFile=$outDir/$sample.svabaTitan.cn.txt
outputBedpeFile=$outDir/$sample.svabaTitan.sv.bedpe  

# combine with copy number data and annotate
cmd="time Rscript $script --id $sample --svaba_funcs $utils --svabaVCF $vcf --titanBinFile $bin --titanSegFile $seg --genomeBuild $genomeBuild --genomeStyle $genomeStyle --minSPAN $minSPAN --minInversionSPAN $minInversionSPAN --outDir $outDir --outputSVFile $outputSVFile --outputCNFile $outputCNFile --outputBedpeFile $outputBedpeFile"
echo $cmd; eval $cmd
    
# add intersecting gene names
cmd="time Rscript $script2 --inputSVFile $outputSVFile --outputSVFile $outputSVFile2"
echo $cmd; eval $cmd
../utils/qc.sh $outputSVFile2 1-2    
	
