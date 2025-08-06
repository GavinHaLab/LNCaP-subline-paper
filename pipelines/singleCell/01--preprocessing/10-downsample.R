# 20241025 arb
# downsample to consistent number of cells

library(Seurat)
library(scuttle)

inFile1 = 'intermediate/09b.txt'  # manifest
inDir = 'intermediate/06/'
outDir = 'intermediate/10/'
outFile1 = 'intermediate/10.txt'

data1 = read.delim(inFile1)
data1a = data1

if(!dir.exists(outDir)) { dir.create(outDir); }

# set goals for downsampling
targetNumCells = 1850
targetAvgNumReads = 19000

files = dir(inDir)
out = NULL
for(file in files) {

  sample = sub('.rds','',file)	 
  print(sample)
  inFile = paste0(inDir,file)
  so = readRDS(inFile)

  # downsample cells
  set.seed(123)
  so1 = subset(so,cells = sample(Cells(so),targetNumCells))

  # downsample counts
  so2 = so1
  counts1 = so2@assays$RNA@layers$counts

  prop = targetAvgNumReads/mean(apply(counts1,2,sum))

  counts2 = downsampleMatrix(counts1,prop = prop, bycol = T)
  so2@assays$RNA@layers$counts = counts2

  # update manifest and qc files
  so2$'nCount_RNA' = Matrix::colSums(counts2)
  so2$'nFeature_RNA' = Matrix::colSums(counts2 > 0)  

  so2 = PercentageFeatureSet(so2, "^MT-", col.name = 'percent_mito')
  so2 = PercentageFeatureSet(so2, "^RP[SL]", col.name = 'percent_ribo')
  so2 = PercentageFeatureSet(so2, "^HB[^(P)]", col.name = 'percent_hb')

  idx = data1a$sample == sample
  data1a$numCells[idx] = ncol(so2)
  data1a$avgNumReads[idx] = mean(so2@meta.data$nCount_RNA)
  data1a$avgNumGenes[idx] = mean(so2@meta.data$nFeature_RNA)  

  # write the new file
  outFile = paste0(outDir,file)
  saveRDS(so2,outFile)
}

# write it yo
write.table(data1a,outFile1,quote=F,row.names=F,sep="\t")


