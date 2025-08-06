# 20241125 arb
# generate barplots of the number of peaks per sample

options(stringsAsFactors=F)

library(ggplot2)

inFile1 = 'intermediate/04a.txt'
outPrefix = 'figs/04-barplot-'

read. = T
if(read.) { data1 = read.delim(inFile1); }

datasets = sort(unique(data1$dataset))

for(dataset in datasets) {
  data1a = data1[data1$dataset == dataset,]

  data1b = table(data1a$sample)
  data1c = data.frame(sample=names(data1b),numPeaks=as.vector(data1b))
  result = kruskal.test(numPeaks ~ sample, data = data1c)
  pval = signif(result$p.value,3)

  # perform single sample t test
  vals = log2(data1c$numPeaks+1)
  mean. = mean(vals)
  diffs = abs(mean. - vals)
  maxDiff = max(diffs)
  idxOutlier = which(diffs == maxDiff)
  vals1 = vals[-idxOutlier]
  outlierVal = vals[idxOutlier]

  result = t.test(vals1,mu = outlierVal)
  pval = signif(result$p.value,3)

  avg = prettyNum(round(mean(data1c$numPeaks)),big.mark=',')

  title = dataset
  subtitle = paste0('Avg = ', avg, ', T P = ', pval)
  plotData = data1c

  theme_set(theme_light(base_size=18))
  p = ggplot(plotData,aes(sample,numPeaks)) + 
      geom_bar(stat='identity',position = 'dodge', fill = 'skyblue') +
      theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
      theme(axis.text = element_text(color = 'black')) +       
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.subtitle = element_text(hjust = 0.5)) +
      labs(x = '', y = 'Number Peaks', title = title, subtitle = subtitle) 

  if(dataset == 'WCM') { p = p + theme(axis.text.x = element_text(size=8)); } 

  filename = paste0(outPrefix,dataset,'.png')
  ggsave(filename,units='in',height=5,width=5,dpi=200)
  print(filename)

}