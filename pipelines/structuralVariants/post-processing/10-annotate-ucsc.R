# 20240927 arb
# annotate variants with mappability tracks from ucsc browser

options(stringsAsFactors=F)

library(data.table)
library(rtracklayer)
library(GenomicRanges)

inFile1 = 'intermediate/09.txt'
outFile1 = 'intermediate/10.txt'

read. = T

if(read.) {
  data1 = fread(inFile1);
  data1 = data.frame(data1)

  beds = list()
  beds[['UMAP100']] = 'input/umap100_20240927.bed'
  beds[['ENCODEBLv2']] = 'input/encodeBlacklistV2_20240927.bed'
  beds[['RepeatMasker']] = 'input/repeatMasker_20240927.bed'

  tracks = list()
  for(bed in names(beds)) {
    print(bed)
    bedFile = beds[[bed]]
    tracks[[bed]] = import(bedFile,format='BED')
  }
  names(tracks) = names(beds)
}

data1a = data1[data1$isSelect == 'Y',]
#data1a = data1
for(track in names(tracks)) {
  print(track)

  field = paste0('is',track)
  data1[[field]] = 'N'

  allMatchedIds = c()
  samples = sort(unique(data1$sample))
  for(sample in samples) {
    print(sample)
    data1b = data1a[data1a$sample == sample,]

    # look for overlap with breakpoints and track
    query1 = GRanges(data1b$chrom1,IRanges(data1b$pos1,data1b$pos1+1))
    query2 = GRanges(data1b$chrom2,IRanges(data1b$pos2,data1b$pos2+1))    
 
    matches1 = findOverlaps(query1,tracks[[track]],type = 'any')
    matches2 = findOverlaps(query2,tracks[[track]],type = 'any')    
    
    matchedQueries1 = unique(queryHits(matches1))
    matchedQueries2 = unique(queryHits(matches2))    
    matchedIds1 = data1b$id[matchedQueries1]
    matchedIds2 = data1b$id[matchedQueries2]    
    allMatchedIds = unique(c(allMatchedIds,matchedIds1,matchedIds2))
  }
  data1[[field]][data1$id %in% allMatchedIds] = 'Y'
}

write.table(data1,outFile1,quote=F,row.names=F,sep="\t")
