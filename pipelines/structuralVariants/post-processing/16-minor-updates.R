# 20241220 arb
# minor added updates

options(stringsAsFactors=F)

inFile1 = 'intermediate/13.bedpe'
outFile1 = 'intermediate/16.bedpe'

data1 = read.delim(inFile1)

data1$genePair = paste0(data1$gene1,',',data1$gene2)

data1$genePair = sapply(data1$genePair, function(x) {
  uniqueValues = unique(sort(strsplit(x,',')[[1]]))
  if(length(uniqueValues) == 2) { uniqueValues = uniqueValues[!is.na(uniqueValues) & uniqueValues != 'NA']; }
  paste(uniqueValues,collapse=',')  
})

data1$isFusion = ifelse(!is.na(data1$gene1) & !is.na(data1$gene2) & data1$gene1 != data1$gene2, 'Y','N')

# if sv is not common with all lines and is not unique then it must be shared in other lines
data1$isNonUnique = ifelse(data1$isCommon == 'N' & data1$isUnique == 'N', 'Y','N')

# decided to exclude
data1a = data1[data1$isSpanningSupport == 'Y',]

# write yo data yo
write.table(data1a,outFile1,quote=F,row.names=F,sep="\t")
