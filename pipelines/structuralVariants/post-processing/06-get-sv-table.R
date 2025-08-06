# 20231125 arb
# combine together merged calls with those that are annotated - we focus on those that are annotatable by combineTitanSvaba script

options(stringsAsFactors=F)
library(dplyr)
library(readxl)
library(data.table)
setDTthreads(2)

inFile1 = 'intermediate/05.txt'
inFile2 = 'input/20230821_CompositeGeneList_PGcurate2.xlsx' # gene list
inFile3 = 'input/robinson.txt'
inDir = 'intermediate/03/results/'
outFile1 = 'intermediate/06.txt'

data1 = read.delim(inFile1)
data2 = read_excel(inFile2)
data2 = data.frame(data2)
data3 = read.delim(inFile3)

pc_genes = data2[,1]
robinson_genes = data3[,1]
robinson_genes = robinson_genes[!grepl('fusions',robinson_genes)]
ets_genes = c('ELF1','ELF2','ELF4','GABPA','ERG','FLI1','FEV','ERF','ETV3','ELF3','ELF5','ESE3','ETS1','ETS2','SDEF','ETV4','ETV5','ETV1','ETV2','SPI1','SPIB','SPIC','ELK1','ELK4','ELK3','ETV6','ETV7')
stopifnot(all(robinson_genes %in% pc_genes))

data1$id = paste0(data1$sample,'__',data1$event_id)

out = NULL
files = dir(inDir,pattern='genes')
# files = files[grep('_170.1',files)]
for(file in files) {
  inFile = paste0(inDir,file)
  print(inFile)  
  data2 = fread(inFile)
  data2 = data.frame(data2)

  data2a = select(data2,sample=Sample,event_id=event,support,orient_1,orient_2,SPAN,type,CN_overlap_type,gene1,gene2,support)
  data2a$id = paste0(data2a$sample,'__',data2a$event_id)
  data2a$sample = NULL  # sample will come from data1
  data2a$event_id = NULL # pull from data1

  # combine tables
  data2b = merge(data2a,data1,by='id')
  stopifnot(nrow(data2a) == nrow(data2b))

  # add annotation
  data2b$is_genic = ifelse(!is.na(data2b$gene1) | !is.na(data2b$gene2), 'Y','N')
  data2b$is_fusion = ifelse(!is.na(data2b$gene1) & !is.na(data2b$gene2), 'Y','N')
  data2b$is_cnv = ifelse(data2b$support == 'CN', 'Y','N')
  data2b$is_translocation = ifelse(data2b$chrom1 != data2b$chrom2, 'Y','N')
  data2b$is_select = 'Y'
  data2b$is_pc_gene = ifelse((!is.na(data2b$gene1) & data2b$gene1 %in% pc_genes) | (!is.na(data2b$gene2) & data2b$gene2 %in% pc_genes), 'Y','N')
  data2b$is_robinson_gene = ifelse((!is.na(data2b$gene1) & data2b$gene1 %in% robinson_genes) | (!is.na(data2b$gene2) & data2b$gene2 %in% robinson_genes), 'Y','N')
  data2b$is_ets_gene = ifelse((!is.na(data2b$gene1) & data2b$gene1 %in% ets_genes) | (!is.na(data2b$gene2) & data2b$gene2 %in% ets_genes), 'Y','N')

  # polish columns
  cols = colnames(data2b)
  prefix = c('sample','event_id','variant_id','mate_id','is_select','is_genic','is_fusion','is_cnv','is_pc_gene','is_robinson_gene','is_ets_gene','is_translocation','type','CN_overlap_type','gene1','gene2')
    cols1 = cols[!cols %in% prefix]
  cols2 = c(prefix,cols1)
  data2b = data2b[,cols2]

  out = rbind(out,data2b)
}

# write it yo
write.table(out,outFile1,quote=F,row.names=F,sep="\t")
