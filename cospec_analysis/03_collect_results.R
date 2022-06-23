library(phangorn)
library(tidyverse)

setwd("/home/sukmb276/Isilon/Metagenomes/projects/ApesComplete/output/220616_analysis_final/cospec_analysis/subgroups")

allstats = readRDS("../groups_all_stats_candidates.Rds")

allgenomes = read.table("../../allgroups/GreatApes.all_genomes.tsv", head=T, stringsAsFactors=F) %>% filter(!is.na(GreatApes_95_final))

tax=read.table("../../allgroups/GreatApes.cleanbins_tax.tsv", head=T, stringsAsFactors=F, sep="\t")
colnames(tax) = c("bin","classification")

tax2 = sapply(tax$classification, function(x) strsplit(x, split=";")[[1]]) %>%
	apply(., 2, function(x){x_low=x[max(grep("__$",x, invert=T))]; x[grepl("__$",x)] = paste0(x[grepl("__$",x)],"unclassified_", x_low); return(x)}) %>%
	t %>% data.frame(row.names=NULL, stringsAsFactors=F, bin=tax$bin, .)
colnames(tax2) = c("bin","kingdom","phylum","class","order","family","genus","species")
tax2$species=gsub(" ","_",tax2$species)

allgenomes$host = sapply(allgenomes$cluster_97_final, function(x) strsplit(x, split="_")[[1]][1])
allgenomes = allgenomes %>% left_join(tax2, by=c("genome"="bin")) %>% select(genome, score, host, GreatApes_95_final, kingdom, phylum, class, order, family, genus, species)


allres=lapply(as.vector(allstats$group), function(s) list(res=read.table(paste0(s,"/",s,".cospec.out"), head=T, stringsAsFactors=F), bionj_tree=read.tree(paste0(s,"/trees/",s,".aln.bionj.bs.nwk")),
  bionj_tree_rooted=read.tree(paste0(s,"/trees/",s,".aln.bionj.rooted.bs.nwk"))))

names(allres) = as.vector(allstats$group)

alldata=list(stats=allstats, results=allres, genomes=allgenomes)

saveRDS(alldata, "../cosepec_alldata.Rds")



############
library(phangorn)
source("~/Isilon/Metagenomes/software/Rscripts/mad.R")

missin=read.table("../rooted_missing.txt", head=F, stringsAsFactors=F)
for(this in missin$V1){
  print(this)
  tree_bionj=read.tree(paste0(this,"/trees/",this, ".aln.bionj.bs.nwk"))
  bionj_tree_rooted<-mad(tree_bionj, output_mode="full")[[6]][[1]]
  write.tree(bionj_tree_rooted,  paste0(this,"/trees/",this, ".aln.bionj.rooted.bs.nwk"))
}
