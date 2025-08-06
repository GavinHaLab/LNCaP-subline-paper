# 20241219 arb
# get table of svs after de-duplication and find out which samples have them
# https://github.com/GavinHaLab/CDK12-CRPC-paper/blob/main/TDPmultiSamplePhylogenetic/TAN_TDPheterogeneity_v1.R
# de-duping does not make sense here because we need each samples sv to be counted so the totals add up
# updated to exclude FGC_SRR from comparison

options(stringsAsFactors=F)
library(dplyr)
library(GenomicRanges)
library(stringr)
source('utils/utils.R')

inFile1 = 'intermediate/12.bedpe'
outFile1 = 'intermediate/13a.txt'
outFile2 = 'intermediate/13b.txt'
outFile3 = 'intermediate/13.bedpe'

data1 = read.delim(inFile1)
data1$isENCODEBLv2 = NULL   # added to previous file to enable addition of _SRR SVs

read. = T
match. = T

if(read.) { 
#  data1a = data1[data1$isLargeSpan == 'Y',]
  data1$svID = paste(data1$chrom1,data1$end1,data1$chrom2,data1$end2,sep="_")
  data1a = data1[data1$sample != 'LNCaP_FGC_SRR7943697',]
  samples = sort(unique(data1a$sample))
}

if(match.) { 
# ** create sample specific matching table **
matchDF1 = NULL
data1b = unique(select(data1a,svID,chrom1,start1,end1,chrom2,start2,end2))
data1b = data1b[order(data1b$svID),]

# samples = c("LNCaP_FGC_SRR7943697")
#samples = c("LNCaP_FGC_SRR7943697")
for(sample in samples) {
  print(sample)
  data1c = data1a[data1a$sample == sample,]	   

  svIDs = data1c$svID

  matches = overlapSVs(data1b,data1c)

  tmp = data.frame(svID=data1b$svID,sample = rep(sample,nrow(data1b)))
  matches$sampleSVID = data1c$svID[matches$subjectHits]
  matches$svID = data1b$svID[matches$queryHits]

  matches2 = select(matches,svID,sampleSVID)

  result1 = merge(tmp,matches2,by='svID')
  stopifnot(nrow(result1) == nrow(matches2))
  result1$matchType = ifelse(result1$svID == result1$sampleSVID,'exact','notExact')
  result1 = result1[,c('sample','sampleSVID','svID','matchType')]

  # add in read support
  readSupport = select(data1c,sampleSVID = svID,SR,DR)
  result1a = merge(result1,readSupport,by='sampleSVID')
  stopifnot(nrow(result1a) == nrow(result1))

  matchDF1 = rbind(matchDF1,result1a)
}
# **
}

# determine if svID has been observed with spanning read support
matchDF1a = matchDF1[matchDF1$SR > 0,]
svIDsSpanningSupport = sort(unique(matchDF1a$svID))
matchDF1$isSpanningSupport = ifelse(matchDF1$svID %in% svIDsSpanningSupport,'Y','N')

# write details table
write.table(matchDF1,outFile1,quote=F,row.names=F,sep="\t")

# create svID specific matching table
matchDF2 = matchDF1
matchDF2$sampleSVID = matchDF2$matchType = NULL
counts = table(matchDF2$svID)
counts = data.frame(svID = names(counts),n=as.numeric(counts))

df1 = matchDF2 %>% group_by(svID) %>% summarize(samples = paste(sample,collapse=','))
df1 = as.data.frame(df1)
# updated samples to be unique

# deal with duplicate samples due to multi-mapping sv matches
df1$samples2 = sapply(df1$samples, function(x) {
  uniqueValues = unique(sort(strsplit(x,',')[[1]]))
  paste(uniqueValues,collapse=',')  
})
df1$n = str_count(df1$samples2,',') + 1

# remove old sample version with dups
df1$samples = df1$samples2
df1$samples2 = NULL

write.table(df1,outFile2,quote=F,row.names=F,sep="\t")

commonIDs = df1$svID[df1$n == length(samples)]
uniqueIDs = df1$svID[df1$n == 1]

data1a$isCommon = ifelse(data1a$svID %in% commonIDs, 'Y','N')
data1a$isUnique = ifelse(data1a$svID %in% uniqueIDs, 'Y','N')
data1a$isSpanningSupport = ifelse(data1a$svID %in% svIDsSpanningSupport, 'Y','N')

write.table(data1a,outFile3,quote=F,row.names=F,sep="\t")
