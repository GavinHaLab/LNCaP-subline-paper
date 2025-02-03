load("/path/to/image/projectName.RData") ## loading matrices from makeMatrixFromTITAN-ICHOR_segBased_hg38.R

samples = rownames(geneCNmat)
cf_threshold_set = c(0.5,0.6,0.7,0.8,0.9) #cellular prevalance threshold -- stick to clonal events
#cf_threshold_set = c(0) 

#load params to retreive ploidy information
params = read.table("/path/to/purity_ploidy.txt", head=T)
params = params[match(samples, params$Sample),]
sample_ploidy = params$Ploidy
#sample_ploidy = round(params$Ploidy)

gene_dat = read.table("/path/to/geneList/withCoords.bed", sep="\t") ## cols should be chr, start, end, gene_symbol, not expecting a header... only really need gene name

gene_list = as.character(gene_dat[!grepl("_", gene_dat$V4, fixed=TRUE),4]) # remove genes with _ in the name, 4th column is gene symbol
match_id = as.vector(match(gene_list, colnames(geneCNmat)))
match_id = match_id[!is.na(match_id)] 
names(match_id) = gene_list


for(cf_threshold in cf_threshold_set){
	print(cf_threshold)
	alteration_list = list()

	for(i in 1:length(match_id)){
		#print(i)
		id = match_id[i]
		#print(id)
		gene = names(id)
		#print(gene)
		chr_gene =  as.character(gene_dat[gene_dat$V4==gene,1])

		cn = geneCNmat[,id]
		loh = geneLOHmat[,id]
		cf = geneCFmat[,id]
		# sample_ploidy = round(params$Ploidy)
		sample_ploidy = params$Ploidy

		#remove segments where cf < threshold (0.5, 0.6,0,7, 0.8) however keep cf=NA (NEUT) and also Keep when cn >=4 #2021.03.11
		subclonal_id = which(cf < cf_threshold & cn <4)
		print(length(subclonal_id))
		
		if(length(subclonal_id)!=0){
			cn = cn[-subclonal_id]
			loh = loh[-subclonal_id]
			cf = cf[-subclonal_id]
			sample_ploidy = sample_ploidy[-subclonal_id]
		}
		##Adjust 'Corrected_Copy_Number' based on ploidy information
		if(chr_gene != "chrX") {cn_adj = cn/sample_ploidy
		} else {cn_adj = cn/(sample_ploidy/2)}


		cn_category = rep(NA,length(samples))
		cn_category[which((loh == 1))] = "LOH"

		cn_category[which((cn_adj >= 2))] = "Amplification"
		cn_category[which((cn_adj >= 2) & loh == 1)] = "Amplification LOH"
		cn_category[which((cn_adj > 0 & cn_adj < 1) & loh == 1)] = "Deletion_LOH"
		cn_category[which((cn_adj >= 0.95 & cn_adj <= 1.05) & loh ==1 )] = "Copy_Neutral_LOH"
		# cn_category[which((cn_adj >= 0 & cn_adj < 1))] = "Deletion"
		cn_category[which(cn_adj == 0)] == "CN_Homozygous_Del"
		# cn_category[which(cn_adj == 0)] = "CN_Homozygous_Del"
		# cn_category[which((cn_adj > 0.75 & cn_adj < 2.5))] = "Baseline"
		cn_category[which((cn_adj > 0 & cn_adj <=0.25))] = "Deep Deletion (3 copies)"
		# cn_category[which((cn_adj > 0.25 & cn_adj <= 0.75))] = "Shallow Deletion (1 copy)"
		cn_category[which((cn_adj > 0.5 & cn_adj <= 0.75))] = "Shallow Deletion (1 copy)"
		cn_category[which((cn_adj > 0.25 & cn_adj <= 0.5))] = "Deletion (2 copies)"

		idx_avail = which(!is.na(cn_category))
		#sample,category,value
		if( length(idx_avail)!=0 )	{
			alteration_list[[gene]] = cbind(names(cn_adj)[idx_avail], gene, cn_category[idx_avail])
		}else{
			message("not available for gene ", gene)
			print(cn_adj)
			print(loh)
			print(cn_category)
		}
	}

	out_alteration_all = do.call(rbind, alteration_list)
	colnames(out_alteration_all) = c("sample","category","value")

	write.csv(out_alteration_all, sprintf("/path/to/output/dir/comut_CN_cf%s.csv",cf_threshold), quote=F)
}
