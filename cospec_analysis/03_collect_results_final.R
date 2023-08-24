library(phangorn)
library(tidyverse)

setwd("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/cospec_analysis/subgroups")

allstats = readRDS("../groups_all_stats_candidates.Rds")

allgenomes = read.table("../../allgroups/GreatApes.all_genomes.tsv", head=T, stringsAsFactors=F) %>% filter(!is.na(GreatApes_95_final))
all_genomes=allgenomes

tax_cluster=read.table("../../allgroups/GreatApes.cluster_final_tax.tsv", head=T, stringsAsFactors=F, sep="\t")
colnames(tax_cluster) = c("bin","classification")
tax=allgenomes %>% left_join(tax_cluster, by=c("GreatApes_95_final" = "bin")) %>% filter(!is.na(classification)) %>% select(bin=genome, classification)

tax2 = sapply(paste0(tax_cluster$bin,";",tax_cluster$classification), function(x){
	x_split=strsplit(x, split=";")[[1]]
	y = x_split[2:8]
	y_low=y[max(grep("__$",y, invert=T))]
	if(max(grep("__$",y, invert=T)) == 6){
		y[7] = paste0(gsub("^g__","s__", y_low),"_sp", gsub("_cluster95_","",x_split[1]))
	}
	if(max(grep("__$",y, invert=T)) == 5){
		y[6] = paste0("g__",gsub("_cluster95_","",x_split[1]))
		y[7] = paste0("s__",gsub("_cluster95_","",x_split[1]),"_sp", gsub("_cluster95_","",x_split[1]))
	}
		return(gsub(" ","_", c(x_split[1],y)))
	}) %>%
	t %>% data.frame(row.names=NULL, stringsAsFactors=F)

colnames(tax2) = c("bin","kingdom","phylum","class","order","family","genus","species")

type_recode=data.frame(short=c("Ggorilla","Gberingei","HsA1","HsA2","HsG","Ptschw","Pttrog","Ptver","Ppan","HsUHGG"),
        short2=c("Ggorilla","Gberingei","Hs","Hs","Hs","Ptschw","Pttrog","Ptver","Ppan","Hs"),

long=c("Gorilla_gorilla_gorilla","Gorilla_beringei","Homo_sapiens","Homo_sapiens","Homo_sapiens","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pan_paniscus","Homo_sapiens"),stringsAsFactors=F)

host = read.table("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/groupings_nozoo_campbell_dk.txt", stringsAsFactors = F)
colnames(host) = c("sample","host")
host = bind_rows(host, data.frame(sample="MGYG00000",host="HsUHGG")) %>% mutate(host = case_when(host=="ZcambellGorilla" ~ "Ggorilla", host=="ZcambellPan"~ "Pttrog", TRUE ~ host),
    host_genus = ifelse(grepl("^P", host),"Pan",ifelse(grepl("^G", host),"Gorilla","Homo")))


all_genomes = all_genomes %>% left_join(tax2, by=c("GreatApes_95_final"="bin"))

all_genomes = all_genomes %>%
  filter(!is.na(GreatApes_95_final), !grepl("Manara", genome)) %>% mutate(sample=ifelse(grepl("^SRR", genome), substr(genome,1,10), substr(genome,1,9))) %>% left_join(host) %>%
  left_join(type_recode, by=c("host"="short")) %>% arrange(-score) %>% distinct_at(.vars=c("long","GreatApes_95_final"), .keep_all=T)


host_rep_genomes = all_genomes %>% arrange(-score) %>% distinct_at(.vars=c("long", "GreatApes_95_final"), .keep_all=T) %>% filter(!is.na(host))

#allstats = allstats[-234,]
#allstats = allstats %>% filter(grepl("Methano|UBA1407", group)==F)

allres=lapply(as.vector(allstats$group), function(s) list(res=read.table(paste0(s,"/",s,".out"), head=T, stringsAsFactors=F),
	bionj_tree=read.tree(paste0(s,"/trees/",s,".iqtree2.ml.nwk")),
  bionj_tree_rooted=read.tree(paste0(s,"/trees/",s,".iqtree2.ml.rooted.nwk")) ))

names(allres) = as.vector(allstats$group)

source("/work_beegfs/sukmb276/Metagenomes/software/Rscripts/mad.R")

host_tree = mad(read.nexus("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/cospec_analysis/host_trees_10kTrees/consensusTree_10kTrees_Primates_Version3_chrono.nex"),output_mode="full")[[6]][[1]]
#type_recode=data.frame(short=c("Ggorilla","Gberingei","HsA1","HsA2","HsG","Ptschw","Pttrog","Ptver","Ppan","HsUHGG"), short2=c("Ggorilla","Gberingei","Hs","Hs","Hs","Ptschw","Pttrog","Ptver","Ppan","Hs"), long=c("Gorilla_gorilla_gorilla","Gorilla_beringei","Homo_sapiens","Homo_sapiens","Homo_sapiens","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pan_paniscus","Homo_sapiens"),stringsAsFactors=F)

micro_tree=read.tree("../../allgroups/GreatApes.bac120+ar53.rooted.tre")

alldata=list(stats=allstats, results=allres, micro_tree=micro_tree, host_tree=host_tree)

saveRDS(alldata, "../cospec_final.Rds")



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
