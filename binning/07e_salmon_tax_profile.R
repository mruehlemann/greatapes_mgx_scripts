library(dplyr)

this=rev(strsplit(getwd(),split="/")[[1]])[1]

sf = read.table(paste0(this,".salmon.out"),head=T, stringsAsFactors=F)
readthresh = sum(sf$NumReads)*0.0001

sf_per_bin = sf %>% mutate(totreads=sum(NumReads), bin = gsub("_contig_[0-9]+", "", Name)) %>%
	group_by(bin) %>% summarize(Length=sum(Length), EffectiveLength=sum(EffectiveLength), NumReads=sum(NumReads)) %>%
	mutate(rpk = NumReads/(EffectiveLength/readthresh)) %>% filter(NumReads>=1000) %>%
	mutate(scaleFactor = sum(rpk)/1000000, tpm=rpk/scaleFactor , sample=this) %>% filter(tpm>=250) %>% mutate(scaleFactor = sum(rpk)/1000000, tpm=rpk/scaleFactor)

tax=read.table("/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/GreatApes_final_tax.tsv", head=T, stringsAsFactors=F, sep="\t")
colnames(tax) = c("bin","classification")

tax2 = sapply(tax$classification, function(x) strsplit(x, split=";")[[1]]) %>%
	apply(., 2, function(x){x_low=x[max(grep("__$",x, invert=T))]; x[grepl("__$",x)] = paste0(x[grepl("__$",x)],"unclassified_", x_low); return(x)}) %>%
	t %>% data.frame(row.names=NULL, stringsAsFactors=F, bin=tax$bin, .)
colnames(tax2) = c("bin","kingdom","phylum","class","order","family","genus","species")

sf_per_bin = sf_per_bin %>% left_join(tax2)

all_lev = lapply(c("kingdom","phylum","class","order","family","genus","species","bin"), function(tl) return(sf_per_bin %>% select(tax=!!tl, tpm) %>% group_by(tax) %>% summarize(tpm=sum(tpm)) %>% mutate(level=tl))) %>%
	do.call("bind_rows",.) %>% mutate(sample=this)

write.table(all_lev, paste0(this,".final_tax.tsv"), row.names=F, sep="\t")
