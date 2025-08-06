# 20241125 arb
# generate 2 tables - one with peaks reproducible merged peaks per line

options(stringsAsFactors=F)

library(data.table)
setDTthreads(2)

inDir = 'intermediate/04'
outFile1 = 'intermediate/04a.txt'

files = dir(inDir)

out = NULL
for(file in files) {

  sample = sub('__.*','',file)
  dataset = sub('.*__','',file)
  dataset = sub('.bed','',dataset)
  print(sample)


  inFile = paste0(inDir,'/',file)
  data1 = fread(inFile,sep="\t",col.names=c('chrom','start','end'))
  data1 = data.frame(data1)

  n = nrow(data1)
  data1a = data.frame(sample = rep(sample,n),dataset = rep(dataset,n),data1)

  if(is.null(out)) { out = data1a; }
  else { out = rbind(out,data1a); } 
}

fwrite(out,outFile1,sep="\t")

