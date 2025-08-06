# 20241214 arb
# generate bedpe table filtered to variants for publication
# added rescues due to manual curation of robinson list

options(stringsAsFactors=F,width=180)

library(data.table)
setDTthreads(2)
library(dplyr)

inFile1 = 'intermediate/10.txt'
inFile2 = 'input/rescues20250118v1.txt'
outFile1 = 'intermediate/12.bedpe'

read. = T

if(read.) { 
  data1 = fread(inFile1)
  data1 = data.frame(data1)
  data2 = read.delim(inFile2)
}
	
# filter for sv calls by 2 or more callers
idx1 = data1$num_sv_callers >= 2

# allow for rescues
data2$id = paste0(data2$sample,'__',data2$event_id)
idxRescues = data1$id %in% data2$id

# filter for svs not in black list
# ** ignore black list region with 1mb of AR on chrx **
start = 67544032
start = start - 1000000
end = 67730619
end = end + 1000000
idx3 = data1$isENCODEBLv2 == 'N' | (data1$chrom1 == 'chrX' & data1$pos1 >= start & data1$pos1 <= end) | (data1$chrom2 == 'chrX' & data1$pos2 >= start & data1$pos2 <= end)
# rescue svs with 1MB of AR

# apply filter
idx = (idx1 & idx3) | idxRescues
data1a = data1[idx,]

# polish
data1a$strand1 = ifelse(data1a$orient_1 == 'fwd','+','-')
data1a$strand2 = ifelse(data1a$orient_2 == 'fwd','+','-')
data1a$start1 = data1a$pos1-1
data1a$start2 = data1a$pos2-1
data1a$end1 = data1a$pos1
data1a$end2 = data1a$pos2
data1a$span = data1a$SPAN
data1a$gene1 = data1a$transectGene1
data1a$gene2 = data1a$transectGene2
data1a$DP=data1a$DP1
data1a$DR=data1a$DR1
data1a$SR=data1a$SR1
data1a$VAF=data1a$VAF1
data1a$isGenic = data1a$is_genic
data1a$isTranslocation = data1a$is_translocation
data1a$svCallers = data1a$sv_callers
data1a$numSVCallers = data1a$num_sv_callers
data1a$eventID = data1a$event_id

# filter for svs that have span > 10k unless a fold back inversion
idx2 = data1a$SPAN == -1 | data1a$SPAN >= 10000 | (!is.na(data1a$CN_overlap_type) & data1a$CN_overlap_type == 'Inversion-FoldBack')
data1a$isLargeSpan = ifelse(idx2,'Y','N')

# **
prefix = c('chrom1','start1','end1','chrom2','start2','end2','strand1','strand2','sample')
cols = c('isLargeSpan','span','type','CN_overlap_type','support','ref1','alt1','ref2','alt2','gene1','gene2','DP','SR','DR','VAF','eventID','isTransectPCGene','isGenic','isTranslocation','svCallers','numSVCallers','isENCODEBLv2')
data1b = data1a[,c(prefix,cols)]

data1b = data1b[order(data1b$sample,data1b$chrom1,data1b$start1),]

# ** add rescue note **
data1b$note = '.'
data1b$id = paste0(data1b$sample,'__',data1b$eventID)
data1b$note[data1b$id %in% data2$id] = 'manual rescue'
data1b$id = NULL
# **

# update sample names to be consistent with manuscript
lookup = c('LNCaP_16D','LNCaP_42D','LNCaP_42F','LNCaP_95','LNCaP_ABL','LNCaP_APIPC','LNCaP_C4','LNCaP_C4-2','LNCaP_C4-2B','LNCaP_FGC_PRJNA361316','LNCaP_FGC','LNCaP_FGC_SRR7943697','LNCaP_shAR','LNCaP_AR907','LNCaP_AR909')
names(lookup) = c('LNCaP_16D','LNCaP_42D','LNCaP_42F','LNCaP_95','LNCaP_Abl','LNCaP_APIPC','LNCaP_C4','LNCaP_C42','LNCaP_C42B','LNCaP_clone_FGC','LNCaP_FGC','LNCaP_FGC_SRR7943697','LNCaP_shAR-pATK','LNCaP-AR-907','LNCaP-AR-909')
data1b$sample = lookup[data1b$sample]

# write yo bedbe yo
write.table(data1b,outFile1,quote=F,row.names=F,sep="\t")