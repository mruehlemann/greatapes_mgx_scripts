library(dplyr)
library(parallel)
args = commandArgs(trailingOnly=TRUE)
salmon_profile = args[1]

sample = strsplit(salmon_profile, split="[.]")[[1]][1]

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



sf = read.table(salmon_profile,head=T, stringsAsFactors=F) %>% mutate(Coverage=NumReads*300/EffectiveLength)  %>% 
   mutate(EffectiveLengthCov=ifelse(Coverage>0.1, 1,0)*EffectiveLength)
readthresh = sum(sf$NumReads)*0.0001
sf_simp = sf %>% filter(NumReads>0)
sf_name_count = rep(sf_simp$Name, times=sf_simp$NumReads)

this=sample


sf_rare_all = lapply(c(100000, 250000, 500000, 1000000, 2500000,5000000,10000000), function(depth){
    print(depth)
    if(depth > sum(sf$NumReads)) return()
    mclapply(1:5, function(i){
        print(i)
        readthresh = depth*0.0001
        sf_rare=sample(sf_name_count, size=depth) %>% data.frame(Name=.) %>% group_by(Name) %>% summarize(NumReadsRare=n()) %>%
            left_join(sf, .) %>% mutate(NumReadsRare = ifelse(is.na(NumReadsRare), 0, NumReadsRare), Coverage=NumReadsRare*300/EffectiveLength)  %>%  
            mutate(EffectiveLengthCov=ifelse(Coverage>0.1, 1,0)*EffectiveLength)
        sf_per_bin = sf_rare %>% mutate(totreads=sum(NumReadsRare), bin = gsub("_contig_[0-9]+", "", Name))%>%
            group_by(bin) %>% summarize(Length=sum(Length), EffectiveLength=sum(EffectiveLength), EffectiveLengthCov=sum(EffectiveLengthCov), Contigs=n(), ContigsCov=sum(Coverage>0.1), NumReadsRare=sum(NumReadsRare)) %>%
            filter(NumReadsRare>=1000, (EffectiveLengthCov/EffectiveLength)>=.2 ) %>% mutate(rpk = NumReadsRare/(EffectiveLength/1000)) %>% 
             mutate(scaleFactor = sum(rpk)/1000000, tpm=rpk/scaleFactor , sample=sample) %>% filter(tpm>=250) %>% mutate(scaleFactor = sum(rpk)/1000000, tpm=rpk/scaleFactor) %>%
            mutate(sample=this, depth=depth, i=i) %>% select(sample, depth, i, bin, tpm)
        return(sf_per_bin)
    }, mc.cores=5) %>% do.call("bind_rows",.)
}) %>% do.call("bind_rows",.)



write.table(sf_rare_all, paste0(this,".salmon_rarefied.tsv"), row.names=F, sep="\t")



# library(tidyverse)

# allsamples = list.files()
# #allsamples = allsamples[!allsamples %in% c("H07680-L1", "H07609-L1", "G02618-L1","H07610-L1","H07632-L1","H07664-L1","H07690-L1")]
# #allsamples = allsamples[grepl("^F",allsamples)==F]
# allrare = lapply(allsamples, function(x){print(x); tfile=paste0(x,"/",x,".salmon_rarefied.tsv"); if(file.exists(tfile)) read.table(tfile, head=T, stringsAsFactors=F)}) %>% do.call("bind_rows",.) %>% data.frame(stringsAsFactors=F)

# write.table(allrare, paste0("../GreatApes.salmon_rarefied.tsv"), row.names=F, sep="\t")


