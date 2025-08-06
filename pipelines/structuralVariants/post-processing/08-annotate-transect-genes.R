# 20240718 arb
# annotate gene transections for each pair of breakpoints

options(stringsAsFactors=F)

library(GenomicRanges)
library(dplyr)
library(readxl)
library(data.table)
setDTthreads(2)

inFile1 = 'intermediate/07.txt'
inFile2 = 'input/20230821_CompositeGeneList_PGcurate2.xlsx'
gene_annotation = "input/GRCh38.p12.ensembl.gene.annotations.sorted.include.AREnhancer.txt"
outFile1 = 'intermediate/08.txt'

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

  # get rid of dups by ordering
  hits = hits[order(hits$gene),]
  hits = hits[!duplicated(hits$queryHits),]

  gene_list[hits$queryHits] = hits$gene
  return(gene_list)
}
# **

# create overlap queries for first and second partner
svs1 = select(data1,id=id,Chromosome = chrom1,Start = pos1)
svs1$End=svs1$Start+1
svs2 = select(data1,id=id,Chromosome = chrom2,Start = pos2)
svs2$End=svs2$Start+1

# get overlapping genes
data1$transectGene1 = getOverlap(svs1,genes)
data1$transectGene2 = getOverlap(svs2,genes)
data1$isTransectPCGene = ifelse(data1$transectGene1 %in% pcGenes | data1$transectGene2 %in% pcGenes, 'Y','N')

# update cols
cols = colnames(data1)
prefix = c('sample','event_id','variant_id','isSelect','isTransectPCGene')
cols = c(prefix,cols[!cols %in% prefix])
data1 = data1[,cols]
data1$is_pc_gene = NULL

# write it yo
fwrite(data1,outFile1,quote=F,row.names=F,sep="\t",na = 'NA')
