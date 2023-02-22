library(dplyr)

tax=read.table("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/GreatApes.cluster_final_tax.tsv", head=T, stringsAsFactors=F, sep="\t") # nolint
colnames(tax) = c("bin","classification")

tax2 = sapply(paste0(tax$bin,";",tax$classification), function(x){
	x_split=strsplit(x, split=";")[[1]]
	y = x_split[2:8]
	y_low=y[max(grep("__$",y, invert=T))]
	if(max(grep("__$",y, invert=T)) == 6){
		y[7] = paste0(gsub("^g__","s__", y_low),"_sp", gsub("_cluster95_","",x_split[1]))
	}
	if(max(grep("__$",y, invert=T)) == 5){ # nolint
		y[6] = paste0("g__",gsub("_cluster95_","",x_split[1]))
		y[7] = paste0("s__",gsub("_cluster95_","",x_split[1]),"_sp", gsub("_cluster95_","",x_split[1]))
	}
		return(gsub(" ","_", c(x_split[1],y)))
	}) %>%
	t %>% data.frame(row.names=NULL, stringsAsFactors=F)

colnames(tax2) = c("bin","kingdom","phylum","class","order","family","genus","species")



for(this in list.files()){

#this=rev(strsplit(getwd(),split="/")[[1]])[1]
print(this)
if(file.exists(paste0(this,"/",this,".salmon.out"))==F) next
sf = read.table(paste0(this,"/",this,".salmon.out"),head=T, stringsAsFactors=F) %>% mutate(Coverage=NumReads*300/EffectiveLength)  %>% 
	mutate(EffectiveLengthCov=ifelse(Coverage>0.1, 1,0)*EffectiveLength)
readthresh = sum(sf$NumReads)*0.0001

sf_per_bin = sf %>% mutate(totreads=sum(NumReads), bin = gsub("_contig_[0-9]+", "", Name))%>%
	group_by(bin) %>% summarize(Length=sum(Length), EffectiveLength=sum(EffectiveLength), EffectiveLengthCov=sum(EffectiveLengthCov), Contigs=n(), ContigsCov=sum(Coverage>0.1), NumReads=sum(NumReads)) %>%
	mutate(rpk = NumReads/(EffectiveLength/readthresh)) %>% filter(NumReads>=1000, (EffectiveLengthCov/EffectiveLength)>=.2 ) %>%
	mutate(scaleFactor = sum(rpk)/1000000, tpm=rpk/scaleFactor , sample=this) %>% filter(tpm>=250) %>% mutate(scaleFactor = sum(rpk)/1000000, tpm=rpk/scaleFactor)


sf_per_bin = sf_per_bin %>% left_join(tax2)

all_lev = lapply(c("kingdom","phylum","class","order","family","genus","species","bin"), function(tl) return(sf_per_bin %>% select(tax=!!tl, tpm) %>% group_by(tax) %>% summarize(tpm=sum(tpm)) %>% mutate(level=tl))) %>%
	do.call("bind_rows",.) %>% mutate(sample=this)

write.table(all_lev, paste0(this,"/",this,".final_tax.tsv"), row.names=F, sep="\t")
}
