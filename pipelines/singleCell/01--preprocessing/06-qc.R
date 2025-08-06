# 20230831 arb
# generate qc updated versions of each dataset
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Calculate_QC
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(Seurat)
library(ggplot2)
library(scds)
library(SingleCellExperiment)

inDir = 'intermediate/05/'
outDir1 = 'intermediate/06/'
outDir2 = 'figs/06/'
outDir3 = 'figs/06/scds/'
outFile5 = 'intermediate/06.txt'

# check/create dir
if(!dir.exists(outDir1)) { dir.create(outDir1); }
if(!dir.exists(outDir2)) { dir.create(outDir2); }
if(!dir.exists(outDir3)) { dir.create(outDir3); } 

files = dir(inDir)
# files = files[1:3]
out = NULL
for(file in files) {
  inFile = paste0(inDir,file)
  so = readRDS(inFile)

  num_cells_before = ncol(so)

  # add qc metrics
  so = PercentageFeatureSet(so, "^MT-", col.name = 'percent_mito')
  so = PercentageFeatureSet(so, "^RP[SL]", col.name = 'percent_ribo')
  so = PercentageFeatureSet(so, "^HB[^(P)]", col.name = 'percent_hb')

  sample = sub('.rds','',file)

  # pre qc filtering view 
  outFile2 = paste0(outDir2,'vlnplot-',sample,'-before.png')
  feats = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
  VlnPlot(so, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3, assay = 'RNA', layer = 'counts') +
    NoLegend()
  ggsave(outFile2)
  print(outFile2)

  # some initial proessing
  so1 = NormalizeData(so)
  so1 = FindVariableFeatures(so1) 
  so1 = ScaleData(so1)
  so1 = RunPCA(so1, npcs = 20)
  so1 = RunUMAP(so1, reduction = 'pca', dims = 1:20, return.model = T)

  # *** run scds to identify doublets ***
  print('running scds...')
  sce1 = as.SingleCellExperiment(so1)
  sce1 = cxds(sce1, retRes = T)
  sce1 = bcds(sce1, retRes = T)
  sce1 = cxds_bcds_hybrid(sce1)
  so2 = as.Seurat(sce1)
  print('done')

  thresholds = c(0,.25,.5,.75,.9,.95,.96,.97,.98,.99,1)
  doubletThreshold = 0.90

  plotData = data.frame(readCount = so2@meta.data$nCount_RNA,hybridScore=so2@meta.data$hybrid_score)
  quartiles = quantile(plotData$hybridScore,thresholds)
  # *** draw doublet visuals ***
  # draw hist
  options(warn=-1)
  xlab = 'SCDS Hybrid Score'
  theme_set(theme_light(base_size=18))  
  p = ggplot(plotData,aes(hybridScore)) +
    geom_histogram(aes(y=..density..),position='identity',bins=50,fill='skyblue',color='grey') + 
    geom_density(alpha = 0.6) + 
    theme(plot.title = element_text(size=12,hjust=0.5)) + 
    theme(axis.text = element_text(color = 'black')) + 
    labs(x=xlab,y = 'Density',title=sample) + 
    geom_vline(xintercept=quartiles,color='red',linetype='dashed')
  outFile1 = paste0(outDir3,sample,'-hist.png')
  ggsave(outFile1,units='in',height=5,width=5,dpi=150)
  print(outFile1)
  options(warn=0)
  
  # draw scatter
  theme_set(theme_light(base_size=18))    
  p = ggplot(plotData,aes(readCount,hybridScore)) +
      theme_set(theme_gray(base_size=12)) +
      theme(axis.text = element_text(color = 'black')) +       
      theme(plot.title = element_text(size=12,hjust=0.5)) +       
      geom_point(shape = 1, color = 'skyblue') +
      geom_hline(yintercept=quartiles,color='red',linetype='dashed')      
      labs(title = sample)
    outFile1 = paste0(outDir3,sample,'-scatter.png')
    ggsave(outFile1,units='in',height=5,width=5,dpi=200)
    print(outFile1)

  # draw boxplots
  thresholdVals = quantile(plotData$hybridScore,thresholds)
  names(thresholdVals) = thresholds
  plotData2 = NULL
  for(c1 in 1:length(thresholdVals)) {
    thresholdVal = thresholdVals[c1]
    threshold = names(thresholdVals)[c1]
    tmp = plotData[plotData$hybridScore >= thresholdVal,]
    label = paste0(threshold,'__',round(thresholdVal,2))
    tmp$label = label

    if(is.null(plotData2)) { plotData2 = tmp; }
    else { plotData2 = rbind(plotData2,tmp); } 
  }

  plotData2$label = factor(plotData2$label)
  theme_set(theme_light(base_size=18))      
  p = ggplot(plotData2,aes(label,readCount)) + 
        geom_boxplot(fill='lightgray') + 
        geom_jitter(height=0,width=.25,shape=1,alpha=.25,color='blue') + 
        theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12)) + 
        theme(axis.title.y = element_text(size=12)) + 
        theme(plot.title = element_text(size=12,hjust=0.5)) +       	
         labs(title=sample, xlab = 'threshold (quantile__value)') 
  outFile1 = paste0(outDir3,sample,'-boxplot.png')
  ggsave(outFile1,units='in',height=5,width=5,dpi=200)
  print(outFile1)
  # ***

  iqr90 = quantile(so2@meta.data$nCount_RNA,.90)
  thresholdVal = quantile(so2@meta.data$hybrid_score,doubletThreshold)
  so2@meta.data$isDoublet = ifelse(so2@meta.data$hybrid_score > thresholdVal & so@meta.data$nCount_RNA > iqr90, 'Y','N')

  # ** draw umap with scds scores **
  plot = FeaturePlot(so2, features = 'hybrid_score', reduction = 'UMAP', pt.size = 1, raster = F)
  outFile1 = paste0(outDir3,sample,'-umap1.png')  
  ggsave(outFile1, plot, width = 13, height = 10, units = 'in', dpi = 600)  
  print(outFile1)
  # **

  # ** draw umap with doublet calls **
  plot = DimPlot(so2, group.by = 'isDoublet', reduction = 'UMAP', pt.size = 1, raster = F,label=F)
  outFile1 = paste0(outDir3,sample,'-umap2.png')    
  ggsave(outFile1, plot, width = 13, height = 10, units = 'in', dpi = 600)  
  print(outFile1)
  # **

  # filter
  # use outliers to get rid of doublets
  iqr75 = quantile(so2@meta.data$nCount_RNA,.75)
  iqr25 = quantile(so2@meta.data$nCount_RNA,.25)  
  iqr = iqr75 - iqr25
  iqr_3 = iqr*3
  outer_fence = iqr75 + iqr_3

  idx1 = so2@meta.data$nFeature_RNA > 200          # low number of detected genes means cells are dead
  idx2 = so2@meta.data$nCount_RNA < outer_fence    # high number of counts means potential doublet
  idx2a = so2@meta.data$nCount_RNA > 10000         # low read count cells may be low quality - in the process of dying - these often cluster by themselves
  idx3 = so2@meta.data$percent_hb < 25             # high proportion of expression in hb genes means red blood cells
  idx4 = so2@meta.data$percent_ribo < 50           # high proportion of ribosomal genes means our library prep may not be effectiv
  idx5 = so2@meta.data$percent_mito < 15           # high proportion of mitochondria expression indicates cells may undergoing apoptosis and are dying...
  idx6 = so2@meta.data$isDoublet == 'N'            # scds 90% quantile cutoffs

  idx = idx1 & idx2 & idx3 & idx4 & idx5 & idx6 & idx2a
  num_cells_after = sum(idx)
  num_features_filter_perc = round((1 - sum(idx1) / num_cells_before) * 100,1)
  RNA_counts_filter_perc = round((1 - sum(idx2) / num_cells_before) * 100,1)
  percent_hb_filter_perc = round((1 - sum(idx3) / num_cells_before) * 100,1)
  percent_ribo_filter_perc = round((1 - sum(idx4) / num_cells_before) * 100,1)
  percent_mito_filter_perc = round((1 - sum(idx5) / num_cells_before) * 100,1)
  percent_doublets = round((1 - sum(idx6) / num_cells_before) * 100,1)      

  so3 = so2[,idx]

  diff = num_cells_before - num_cells_after
  diff_perc = round(diff/num_cells_before*100,1)

  lineOut = c(sample,num_cells_before,num_cells_after,diff_perc,diff,num_features_filter_perc,RNA_counts_filter_perc,percent_hb_filter_perc,percent_ribo_filter_perc,percent_mito_filter_perc,percent_doublets)
  out = rbind(out,lineOut)

  # post filtering picture
  outFile2 = paste0(outDir2,'vlnplot-',sample,'-after.png')
  feats = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb",'hybrid_score')
  VlnPlot(so3, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3, assay = 'RNA', layer = 'counts') +
    NoLegend()
  ggsave(outFile2)
  print(outFile2)

  # remove reductions as we will likely use different methods laters
  so3[['UMAP']] = so3[['PCA']] = NULL

  # qc involves creating separate umaps... need to re-iniitalize **
  metaData = so3@meta.data
  counts = so3@assays$RNA$counts
  so4 = CreateSeuratObject(counts = counts,meta.data = metaData)

  # write the new file
  outFile = paste0(outDir1,file)
  saveRDS(so4,outFile)
}

# write it yo
out = data.frame(out)
colnames(out) = c('sample','num_cells_before','num_cells_after', 'percent_difference','difference','filter_percent_num_features','filter_percent_RNA_counts','filter_percent_hb','filter_percent_ribo','filter_percent_mito','percent_scds_doublets')
write.table(out,outFile5,quote=F,row.names=F,sep="\t")
