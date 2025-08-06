# 20230824 arb
# read each dataset and create a seurate object

library(Seurat)

inFile1 = 'intermediate/04.txt'
outDir = 'intermediate/05'

 # check/create dir
if(!dir.exists(outDir)) { dir.create(outDir); } 

data1 = read.delim(inFile1)

for(c1 in 1:nrow(data1)) {
  sample = data1$sample[c1]
  expression_matrix = Read10X(data1$path[c1])
  so = CreateSeuratObject(counts = expression_matrix, project = sample)
  outFile = paste0(outDir,'/',sample,'.rds')
  saveRDS(so, file=outFile)
  print(sample)
}









