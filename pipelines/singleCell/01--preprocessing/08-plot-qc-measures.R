# 20240925 arb
# create per sample boxplots of qc metrics

options(stringsAsFactors=F)
library(ggplot2)

inFile1 = 'intermediate/07.txt'
outFile1 = 'figs/08-boxplot-numReads.png'
outFile2 = 'figs/08-boxplot-numGenes.png'

read.=T
if(read.) { data1 = read.delim(inFile1); }

plotData = data1

plotData$id = paste0(plotData$sample,'__',plotData$cellID)
numSamples=length(unique(plotData$sample))
numCells=prettyNum(length(unique(plotData$id)),big.mark=',')


avg = prettyNum(round(mean(plotData$numReads),1),big.mark=',')
title = paste0("Number of Reads Per Cell, Avg=",avg) 
subtitle = paste0(numSamples,' samples, ', numCells, ' cells')

p1 = ggplot(plotData,aes(x=sample,y=numReads)) +
  theme_light(base_size=18) +
  theme(axis.text = element_text(color = 'black')) +   
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=10)) +
  theme(plot.title = element_text(size=14,hjust=.5)) +
  theme(plot.subtitle = element_text(size=14,hjust=.5)) +     
  geom_boxplot(fill = 'skyblue') +
  labs(x = '', y = 'Number Reads',title=title,subtitle=subtitle) +
  scale_y_log10()

ggsave(outFile1,units='in',height=5,width=10,dpi=300)
print(outFile1)

avg = prettyNum(round(mean(plotData$numGenes),1),big.mark=',')
title = paste0("Number of Genes Per Cell, Avg=",avg) 
subtitle = paste0(numSamples,' samples, ', numCells, ' cells')

p1 = ggplot(plotData,aes(x=sample,y=numGenes)) +
  theme_light(base_size=18) +
  theme(axis.text = element_text(color = 'black')) +   
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=10)) +
  theme(plot.title = element_text(size=14,hjust=.5)) +
  theme(plot.subtitle = element_text(size=14,hjust=.5)) +     
  geom_boxplot(fill = 'skyblue') +
  labs(x = '', y = 'Number Genes',title=title,subtitle=subtitle) +
  scale_y_log10()

ggsave(outFile2,units='in',height=5,width=10,dpi=300)
print(outFile2)