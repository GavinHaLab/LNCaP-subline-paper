# 20240127 arb
# verison 2b - do not return only var genes 
# https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
# using ^^ "first scenario in which we normalize and then perform sct transform on all samples then run harmony

library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(dbscan)

read. = T
process. = T
umap. = T
cluster. = T
plot. = T
write. = T

inFile1 = 'input/manifest.txt'
inDir = 'input/sos/'
outFile1 = 'intermediate/02.rds'
outFilePrefix2 = 'figs/02-'
project = 'lncap'

data1 = read.delim(inFile1)

if(read.) { 

meta_cols = c('line','sample','sampleID')	
files = dir(inDir)
sos = list()
for(file in files) {
  sample = sub('.rds','',file)
  inFile = paste0(inDir,file)  
  so = readRDS(inFile)

  # add metadata
  for(meta_col in meta_cols) {
    so@meta.data[[meta_col]] = data1[[meta_col]][data1$sample == sample]
  }

  # update cell names for uniqueness across datasets
  so = RenameCells(so,add.cell.id = sample)

  # qc involves creating separate umaps... need to re-iniitalize **
  metaData = so@meta.data
  counts = so@assays$RNA$counts
  so = CreateSeuratObject(counts = counts,meta.data = metaData)

  sos[[sample]] = so
  print(sample)
}
}

if(process.) { 

  for(c1 in 1:length(sos)) {
    sample = names(sos)[c1]
    print(sample)
    so = sos[[c1]]

    # process
    so = NormalizeData(so)

    options(warn=-1)
    so = CellCycleScoring(so, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
    options(warn=0)

    # ** regress out cell cycle genes **
    so = FindVariableFeatures(so, nfeatures = 3000)
    so = ScaleData(so)
#    so = ScaleData(so,vars.to.regress = c('S.Score','G2M.Score'), features = rownames(so))
    # **

    so = RunPCA(so, npcs = 20)
    so = FindNeighbors(so, reduction = 'pca', dims = 1:20)
    so = FindClusters(so, resolution = .1)  # initial default cluster resolution  
    so = RunUMAP(so, reduction = 'pca', dims = 1:20, return.model = T) 
    
    # store it 
    sos[[c1]] = so
  }

}

if(umap.) { 


  for(c1 in 1:length(sos)) {
    sample = names(sos)[c1]
    print(sample)
    so = sos[[c1]]
    so = FindClusters(so, resolution = .5)  # initial default cluster resolution  
    so = RunUMAP(so, reduction = 'pca', dims = 1:20, return.model = T) #, min.dist = 0.1) # default min.dist = 0.3

    # store it 
    sos[[c1]] = so
  }

}

if(cluster.) { 

  print('clustering...')
  for(c1 in 1:length(sos)) {
    
    sample = names(sos)[c1]
    print(sample)
    so = sos[[c1]]
    set.seed(123)

    umap. = Embeddings(so, "umap")
    result = dbscan(umap.,eps=.425,minPts = 20,borderPoints = T) # eps = .4, 10 for c42b and fgcp17 .6/40 , .4/15, .5/25
    so@meta.data[['dbscanCluster']] = result$cluster

   sos[[sample]] = so
  }

}

if(plot.) { 
  for(c1 in 1:length(sos)) {
    sample = names(sos)[c1]
    print(sample)
    so = sos[[c1]]
    set.seed(123)

    # plot
    metas = c('dbscanCluster','seurat_clusters','kMeansCluster','Phase','kMeansCluster2')
    metas = c('dbscanCluster','seurat_clusters','Phase') 
    for(meta in metas) {
      plot = DimPlot(so, reduction = 'umap', group.by = meta, pt.size = 1.5, raster = F, label = T)  # pt.size = 0.1
      filename = paste0(outFilePrefix2,sample,'-dimplot-',meta,'.png')
      ggsave(filename, plot, width = 11, height = 10, units = 'in', dpi = 600)
      print(filename)
    }
  }
}

if(write.) { saveRDS(sos,outFile1); } 
