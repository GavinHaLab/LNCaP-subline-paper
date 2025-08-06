# 20240706 arb

library(data.table)
setDTthreads(2)

inFile1 = 'intermediate/06.txt'
outFile1 = 'intermediate/07.txt'

data1 = fread(inFile1)
data1 = data.frame(data1)

data1$VAF1 = (data1$SR1+data1$DR1)/data1$DP1
data1$VAF2 = (data1$SR2+data1$DR2)/data1$DP2
data1$is_select = NULL

# ** baseline select **
minSpan = 1000
numCallers = 1

idxSpan = data1$SPAN == -1 | data1$SPAN >= minSpan
idxCallers = data1$num_sv_callers >= numCallers

data1$isSelect = ifelse(idxSpan & idxCallers, 'Y','N')
# **

# ** 10k span **
minSpan = 10000
idxSpan = data1$SPAN == -1 | data1$SPAN >= minSpan
data1$isFilterSpan = ifelse(idxSpan, 'Y','N')
# **

# ** 2 callers **
numCallers = 2
idxCallers = data1$num_sv_callers >= numCallers
data1$isFilterCallers = ifelse(idxCallers,'Y','N')
# **

# ** min vaf **
minVAF = 0.1
idxVAF = data1$VAF1 >= minVAF & data1$VAF2 >= minVAF
data1$isFilterVAF = ifelse(idxVAF, 'Y','N')
# **

# remove SRR
# data1a = data1[data1$sample != 'LNCaP_FGC_SRR7943697',]
data1a = data1

fwrite(data1a,outFile1,quote=F,row.names=F,sep="\t",na='NA')
