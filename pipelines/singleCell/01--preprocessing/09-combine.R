# 20241009 arb
# combine together qc metrics into a single table

options(stringsAsFactors=F)

library(dplyr)

inFile1 = 'intermediate/04.txt'  # old parsed manifest with paths
inFile2 = 'intermediate/06.txt'  # qc file
inFile3 = 'intermediate/07.txt'  # new metrics
outFile1 = 'intermediate/09a.txt' # qc
outFile2 = 'intermediate/09b.txt' # manifest

data1 = read.delim(inFile1)
data2 = read.delim(inFile2)
data3 = read.delim(inFile3)

# aggregate averages per sample
data3a = aggregate(numReads ~ sample, data = data3, FUN = mean)
data3b = aggregate(numGenes ~ sample, data = data3, FUN = mean)

# join and rename for qc
data1a = select(data1,sample,sampleID,line,SCSN,folder,path)
data2a = select(data2,sample,numCells = num_cells_after,startingNumCells = num_cells_before)
data1b = merge(data1a,data2a,by='sample')
data1c = merge(data1b,data3a,by='sample')
data1d = merge(data1c,data3b,by='sample')
stopifnot(nrow(data1d) == nrow(data1a))

# polish and write manifest
data1d = rename(data1d,avgNumReads = numReads,avgNumGenes = numGenes)
cols = c('sample','line','sampleID','SCSN','numCells','startingNumCells','avgNumReads','avgNumGenes','folder','path')
data1e = data1d[,cols]
write.table(data1e,outFile2,quote=F,row.names=F,sep="\t")

# wrangle and write
data2a = merge(data2,data3a,by='sample')
data2b = merge(data2a,data3b,by='sample')
write.table(data2b,outFile1,quote=F,row.names=F,sep="\t")
