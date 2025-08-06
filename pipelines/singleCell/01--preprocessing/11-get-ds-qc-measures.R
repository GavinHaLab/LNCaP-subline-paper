# 20240925 arb
# get table of sample qc values to plot

library(Seurat)
library(dplyr)

inDir = 'intermediate/10/'
outFile1 = 'intermediate/11.txt'

# check/create dir

files = dir(inDir)
# files = files[1:3]
out = NULL
for(file in files) {
  inFile = paste0(inDir,file)
  so = readRDS(inFile)
  print(inFile)

  so@meta.data$cellID = rownames(so@meta.data)
  tmp = select(so@meta.data,sample=orig.ident,cellID,numReads=nCount_RNA,numGenes=nFeature_RNA,percentMito=percent_mito)  
  out = rbind(out,tmp)
}

# write it yo
out = data.frame(out)
write.table(out,outFile1,quote=F,row.names=F,sep="\t")
