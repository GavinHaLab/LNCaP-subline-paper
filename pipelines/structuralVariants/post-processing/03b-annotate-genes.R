# 20231123 arb
# annotate breakpoints with gene information

options(stringsAsFactors=F)

library(optparse)
library(GenomicRanges)
library(dplyr)

# ** deal with arguments **
option_list <- list(
	make_option(c("--inputSVFile"), type = "character", help = "output file generated from combineSVABAandTitan.R script"),
	make_option(c("--outputSVFile"), type = "character", help = "name of updated file with gene annotations")
)
parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)

inFile1 <- opt$inputSVFile
#inFile1 = "intermediate/04/PCFHX_2001B_P29_CL.svabaTitan.sv.txt"
outFile1 <- opt$outputSVFile
#outFile1 = "intermediate/04/PCFHX_2001B_P29_CL.svabaTitan.genes.sv.txt"
gene_annotation = "/fh/fast/ha_g/WCDT_CRPC_WGS/Analysis/matrices_CNA/gene/GRCh38.p12.ensembl.gene.annotations.sorted.include.AREnhancer.txt"
# **

data1 = read.delim(inFile1)

# ** read gene annotation **
genomeStyle = "UCSC"
chrs = c(1:22,"X")
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
svs1 = select(data1,id=SV.id,Chromosome = chromosome_1,Start = start_1)
svs1$End=svs1$Start+1
svs2 = select(data1,id=SV.id,Chromosome = chromosome_2,Start = start_2)
svs2$End=svs2$Start+1

# get overlapping genes
data1$gene1 = getOverlap(svs1,genes)
data1$gene2 = getOverlap(svs2,genes)

print(paste0('writing ', outFile1, ' ...'))
write.table(data1,outFile1,quote=F,row.names=F,sep="\t")
