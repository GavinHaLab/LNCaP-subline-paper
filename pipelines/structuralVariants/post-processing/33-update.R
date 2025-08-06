# 20250219 arb
# last minute changes:
# 1. remove AR907/909 AR SVs - these are the result of AR vector

options(stringsAsFactors=F)

inFile1 = 'intermediate/31.bedpe'
outFile1 = 'intermediate/33.bedpe'

# read it
data1 = read.delim(inFile1)

# filter
idx = data1$sample %in% c('LNCaP_AR907','LNCaP_AR909') & !is.na(data1$genePair) & data1$genePair == 'AR'
data1a = data1[!idx,]

# write it yo
write.table(data1a,outFile1,quote=F,row.names=F,sep="\t")