# 20241120 arb
# differential expression

options(stringsAsFactors=F,width=180)

library(Seurat)
library(SeuratWrappers)

inFile1 = 'intermediate/02.rds'
outFile1 = 'intermediate/04.txt'

read. = T

ccgenes = sort(unique(c(cc.genes$s.genes,cc.genes$g2m.genes)))

if(read.) { sos = readRDS(inFile1); }

out = out2 = NULL
labeling = 'Phase'

for(c1 in 1:length(sos)) {
so = sos[[c1]]
sample = names(sos)[c1]
print(sample)

clusters = sort(unique(so@meta.data[[labeling]]))
for(cluster in clusters) {
  print(cluster)
  cells = rownames(so@meta.data)[so@meta.data[[labeling]] == cluster]

  result = RunPresto(so,ident.1 = cells, assay = 'RNA')  # use much faster presto implementation 
  result1 = data.frame(sample=rep(sample,nrow(result)),cluster=rep(cluster,nrow(result)),gene = rownames(result),result)
  out = rbind(out,result1)
}
}

# order and write
out = as.data.frame(out)
out$order = as.numeric(out$pct.1)*as.numeric(out$avg_log2FC)*abs(as.numeric(out$pct.1)-as.numeric(out$pct.2))
out = out[order(out$sample,out$cluster,-out$order),]

idx = out$pct.1 > 0.5 & out$p_val_adj < 0.05 & out$avg_log2FC > 0
out$isSig = ifelse(idx,'Y','N')
out$isCCgene = ifelse(out$gene %in% ccgenes,'Y','N')

prefix = c('sample','cluster','gene','isSig','isCCgene','avg_log2FC','pct.1','pct.2','p_val')
cols = colnames(out)
cols1 = c(prefix,cols[!cols %in% prefix])
out = out[,cols1]

write.table(out,outFile1,quote=F,row.names=F,sep="\t")
