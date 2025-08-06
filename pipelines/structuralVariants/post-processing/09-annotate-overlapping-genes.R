# 20240718 arb
# annotate intra chrom SVs based on what genes they overlap

options(stringsAsFactors=F)
library(GenomicRanges)
library(dplyr)
library(data.table)
setDTthreads(2)
library(readxl)
library(dplyr)

inFile1 = 'intermediate/08.txt'
inFile2 = 'input/20230821_CompositeGeneList_PGcurate2.xlsx'
gene_annotation = "input/GRCh38.p12.ensembl.gene.annotations.sorted.include.AREnhancer.txt"
outFile1 = 'intermediate/09.txt'

read = T
if(read) {
  data1 = fread(inFile1)
  data1 = data.frame(data1)
  data2 = read_excel(inFile2,skip = 1)
  data2 = data.frame(data2)
  pcGenes = data2[,1]
}

# ** read gene annotation **
genomeStyle = "UCSC"
chrs = c(1:22,"X",'Y')
seqlevelsStyle(chrs) = genomeStyle

genes = read.delim(gene_annotation)
genes$Chromosome=paste0('chr',genes$Chr)
genes$Start = genes$cdsStart
genes$End = genes$cdsEnd
genes$length= genes$End - genes$Start

# keep longest transript
genes = genes[order(genes$length,decreasing=T),]
genes = genes[!duplicated(genes$Ensg_ID) & genes$Chromosome %in% chrs,]
genes = select(genes,ensembl_gene = Ensg_ID,Chromosome,Start,End,gene=HGNC_symbol)

# ** sourced from /fh/fast/ha_g/user/gha/software/git/TitanCNA/R/utils.R **
getOverlap <- function(svs, genes) {
  gene_list = rep(NA,nrow(svs))	   

  x <- as(svs, "GRanges")
  y <- as(genes, "GRanges")
  hits <- findOverlaps(query = x, subject = y, type = 'any')

  # polish a data frame
  hits = as.data.frame(hits)
  hits$sv_idx = hits$queryHits

  # get genes and remove the blanks
  hits$gene = genes$gene[hits$subjectHits]
  hits = hits[hits$gene != '',]

  hits = hits[order(hits$sv_idx,hits$gene),]
  return(hits)
#  hits = hits[!duplicated(hits$queryHits),]

#  gene_list[hits$queryHits] = hits$gene
#  return(gene_list)
}
# **

# get intra chrom events
data1a = data1[data1$chrom1 == data1$chrom2 & data1$chrom1 != -1 & data1$chrom2 != -1,]
segs = select(data1a,id,Chromosome=chrom1,Start=pos1,End=pos2)
result = getOverlap(segs,genes)
result = data.frame(result)
result$isPCGene = ifelse(result$gene %in% pcGenes, 'Y','N')

# ** add cols isPCGeneOverlap overlapGenes overlapPCGenes **
result$id = data1a$id[result$sv_idx]
result1 = sort(unique(result$id[result$isPCGene == 'Y']))
data1$isPCGeneOverlap = ifelse(data1$id %in% result1, 'Y','N')

result2 = aggregate(gene ~ id, data = result, FUN = paste,collapse=',')
result2 = data.frame(result2)
colnames(result2) = c('id','overlapGenes')

result3 = result[result$isPCGene == 'Y',]
result3a = aggregate(gene ~ id, data = result3, FUN = paste,collapse=',')
result3a = data.frame(result3a)
colnames(result3a) = c('id','overlapPCGenes')

data1a = merge(data1,result2,by='id',all.x = T)
data1a = merge(data1a,result3a,by='id',all.x = T)

# update cols
cols = colnames(data1a)
prefix = c('sample','event_id','variant_id','isSelect','isTransectPCGene','isPCGeneOverlap')
cols = c(prefix,cols[!cols %in% prefix])
data1a = data1a[,cols]

# write it yo
fwrite(data1a,outFile1,quote=F,row.names=F,sep="\t",na="NA")
