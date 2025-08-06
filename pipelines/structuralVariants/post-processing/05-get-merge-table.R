# 20231125 arb
# get table with sv merge details - this is a super set of merged sv calls from our callers

options(stringsAsFactors=F)

library(dplyr)
library(data.table)
setDTthreads(2)

inDir='intermediate/02/results/'
outFile1 = 'intermediate/05.txt'

out = NULL
files=dir(inDir,pattern='sv-merge.txt')
# files = files[1:3]
for(file in files) {
  print(file)	   
  inFile = paste0(inDir,file)
  data1 = fread(inFile)
  data1 = data.frame(data1)

  data1a = select(data1,sample,event_id,variant_id,mate_id,sv_callers=callers,num_sv_callers=num_callers,chrom1=chrom,pos1=pos,DP1=tumor_dp,DR1=tumor_discordant_rs,SR1=tumor_spanning_rs,ref1=ref,alt1=alt)
  data1a$id = paste0(data1a$sample,'__',data1a$event_id)

  data1a = data1a[order(data1a$event_id,data1a$variant_id),]  

  # filter to first variant_id per event
  data1b = data1a[!duplicated(data1a$id),]  # first variant_id
  data1c = data1a[duplicated(data1a$id),]  # second variant_id
  #  stopifnot(nrow(data1b) == nrow(data1c))  # due to gridss it is no longer the case that we always have 2 - singletons should not have been annotated...

  # prep partner to be joined
  data1d = select(data1c,id,chrom2=chrom1,pos2=pos1,DP2=DP1,DR2=DR1,SR2=SR1,ref2=ref1,alt2=alt1)
  data1e = merge(data1b,data1d,by='id')
#  stopifnot(nrow(data1e) == nrow(data1b))
  stopifnot(nrow(data1e) == nrow(data1c))  

  # polish
  data1e = data1e[order(data1e$id),]
  data1e$id = NULL
  
  out = rbind(out,data1e)
}

fwrite(out,outFile1,quote=F,row.names=F,sep="\t",na="NA")