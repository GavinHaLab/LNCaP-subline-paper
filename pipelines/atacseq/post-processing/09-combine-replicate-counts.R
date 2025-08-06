# 20231222 arb
# combine replicate so we have line counts

options(stringsAsFactors=F)

inDir = 'intermediate/08/'
inFile1 = 'intermediate/02.txt'
outDir = 'intermediate/09/'

if(!file.exists(outDir)) { dir.create(outDir); } 

read = T

data1 = read.delim(inFile1)
data1 = data1[order(data1$sampleID2),]

files=dir(inDir)
#files = c(files[3])
for(file in files) {
  antibody = sub('.txt','',file)
  print(antibody)
  data1a = data1[grep(antibody,data1$sampleID),]

  if(read) { 
    inFile2 = paste0(inDir,file)
    data2 = read.delim(inFile2)
    colnames(data2) = sub('^X','',colnames(data2))
  }

  # get new matrix
  sampleID2s = unique(data1a$sampleID2)
  new_mat = matrix(NA,nrow=nrow(data2),ncol=(length(sampleID2s) + 1))
  colnames(new_mat) = c('peak',sampleID2s)
  new_mat[,1] = data2$peak

  for(sampleID2 in sampleID2s) {
    print(paste0('.',sampleID2))
    sampleIDs = data1a$sampleID[data1a$sampleID2 == sampleID2]

    # single replicate 
    if(length(sampleIDs) == 1) { new_mat[,sampleID2] = data2[,sampleIDs]; }
    # multiple replicates 
    else { new_mat[,sampleID2] = apply(data2[,sampleIDs],1,sum); } 
  }

  colnames(new_mat) = sub(paste0('__',antibody),'',colnames(new_mat))

  outFile1 = paste0(outDir,file)
  write.table(new_mat, outFile1, quote=F,row.names=F,sep="\t")
}
