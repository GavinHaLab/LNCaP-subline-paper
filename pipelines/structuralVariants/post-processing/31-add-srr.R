# 20250123 arb
# add srr back in

options(stringsAsFactors=F)

inFile1 = 'intermediate/16.bedpe'
inFile2 = 'intermediate/12.bedpe'
outFile1 = 'intermediate/31.bedpe'

data1 = read.delim(inFile1)
data2 = read.delim(inFile2)

data2a = data2[data2$sample == 'LNCaP_FGC_SRR7943697',]

data2a$svID = paste(data2a$chrom1,data2a$end1,data2a$chrom2,data2a$end2,sep="_")
data2a$isCommon = NA
data2a$isUnique = NA
data2a$isNonUnique = NA
data2a$isSpanningSupport = ifelse(data2a$svID %in% data1$svID | data2a$SR != 0,'Y','N')  # based on whether we have seen it before (exact) in other samples
data2a$isFusion = ifelse(!is.na(data2a$gene1) & !is.na(data2a$gene2) & data2a$gene1 != data2a$gene2, 'Y','N')

# create genePair column
data2a$genePair = paste0(data2a$gene1,',',data2a$gene2)
data2a$genePair = sapply(data2a$genePair, function(x) {
  uniqueValues = unique(sort(strsplit(x,',')[[1]]))
  if(length(uniqueValues) == 2) { uniqueValues = uniqueValues[!is.na(uniqueValues) & uniqueValues != 'NA']; }
  paste(uniqueValues,collapse=',')  
})

# get make sure we can 

missing = colnames(data1)[!colnames(data1) %in% colnames(data2a)]
stopifnot(length(missing) == 0)

# filter for svs not in black list
# ** ignore black list region with 1mb of AR on chrx **
start = 67544032
start = start - 1000000
end = 67730619
end = end + 1000000
idx1 = data2a$isENCODEBLv2 == 'N' | (data2a$chrom1 == 'chrX' & data2a$end1 >= start & data2a$end1 <= end) | (data2a$chrom2 == 'chrX' & data2a$end2 >= start & data2a$end2 <= end)
idx2 = data2a$numSVCallers >= 2
idx3 = data2a$isSpanningSupport == 'Y'

data2b = data2a[idx1 & idx2 & idx3,]

cols = colnames(data1)
data2c = data2b[,cols]

data1a = rbind(data1,data2c)

write.table(data1a,outFile1,quote=F,row.names=F,sep="\t")
