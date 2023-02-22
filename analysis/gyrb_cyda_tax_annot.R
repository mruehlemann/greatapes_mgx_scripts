setwd("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups")
library(tidyverse)

rep_list=read.table("GreatApes.representatives_list", head=F, stringsAsFactors=F)
manara_hits = read.table("GreatApes.representatives.manara_mash.tsv",head=F, stringsAsFactors=F) %>%
  left_join(rep_list %>% mutate(V1=gsub("/home/sukmb276/Isilon","/work_beegfs/sukmb276", V1)), by=c("V2"="V1")) %>% filter(V3<=.065) %>%
  select(V2.y)

tax=read.table("GreatApes.cluster_final_tax.tsv", head=T, stringsAsFactors=F, sep="\t")
colnames(tax) = c("bin","classification")

tax2 = sapply(paste0(tax$bin,";",tax$classification), function(x){
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

colnames(tax2) = c("GreatApes_95_final","kingdom","phylum","class","order","family","genus","species")

type_recode=data.frame(short=c("Ggorilla","Gberingei","HsA1","HsA2","HsG","Ptschw","Pttrog","Ptver","Ppan","HsUHGG"), 
	short2=c("Ggorilla","Gberingei","Hs","Hs","Hs","Ptschw","Pttrog","Ptver","Ppan","Hs"), 
	long=c("Gorilla_gorilla_gorilla","Gorilla_beringei","Homo_sapiens","Homo_sapiens","Homo_sapiens","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pan_paniscus","Homo_sapiens"),stringsAsFactors=F)


host = read.table("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/groupings_nozoo_campbell.txt", stringsAsFactors = F) %>% mutate(V2=ifelse(V2=="ZcambellGorilla","Ggorilla", V2), V2=ifelse(V2=="ZcambellPan","Pttrog", V2))
colnames(host) = c("sample","host")
host = bind_rows(host, data.frame(sample="MGYG00000",host="HsUHGG"))

all_genomes=read.table("GreatApes.all_genomes.tsv",head=T, stringsAsFactors=F)
all_genomes = all_genomes %>% left_join(tax2)
all_genomes = all_genomes %>% mutate(sample=ifelse(grepl("^SRR", genome), substr(genome,1,10), substr(genome,1,9))) %>% left_join(host) %>% left_join(type_recode, by=c("host"="short")) %>% filter(sample!="ManaraSGB")


#### protein catalogue
library(data.table)
prot_to_genome = read_delim('GreatApes.prot_to_genome.tsv', delim='\t', quote='', num_threads=10,col_names=c('mag_id','prot_id'))
prot_to_cluster = read_delim('protein_catalogs/GreatApes.prot95_cluster.tsv', delim='\t', col_names=c('cluster_id','prot_id'))
cluster_annot = readRDS('protein_catalogs/GreatApes.emapper.Rds')

cluster_to_sgb = ag_clean %>% select(mag_id = genome, GreatApes_95_final, host) %>% left_join(prot_to_genome) %>% left_join(prot_to_cluster) %>% mutate(host_group=substr(host,1,1)) %>% select(GreatApes_95_final, mag_id, host_group, prot_id, cluster_id) %>% distinct()

cyda_cluster = cluster_annot %>% filter(Preferred_name == "cydA") %>% pull(cluster_id)
cydb_cluster = cluster_annot %>% filter(Preferred_name == "cydB") %>% pull(cluster_id)

cyda_prot = cluster_to_sgb %>% filter(cluster_id %in% cyda_cluster) %>% left_join(tax2) %>% filter(genus %in% c("g__Prevotella", "g__Paraprevotella"))
cydb_prot = cluster_to_sgb %>% filter(cluster_id %in% cydb_cluster) %>% left_join(tax2) %>% filter(genus %in% c("g__Prevotella", "g__Paraprevotella"))


dir.create("../analyses/cydab", recursive = T)
write.table(cyda_prot, "../analyses/cydab/cyda_prot_to_genome.tsv", sep="\t", quote=F)
write.table(cydb_prot, "../analyses/cydab/cydb_prot_to_genome.tsv", sep="\t", quote=F)



### get all gyrB annotated protein clusters; n=1329
gyrb_clusters = cluster_annot %>% filter(Preferred_name == "gyrB") %>% pull(cluster_id)
### get all gyrB proteins; n=7916
gyrb_prot = prot_to_cluster %>% filter(cluster_id %in% gyrb_clusters) %>% pull(prot_id)
### get all gyrB to genome mappings
gyrb_prot_to_genome = prot_to_genome %>% filter(prot_id %in% gyrb_prot) %>% left_join(all_genomes %>% select(mag_id=genome, host, GreatApes_95_final)) %>% left_join(tax2)

dir.create("../analyses/gyrb", recursive = T)
write.table(gyrb_prot_to_genome, "../analyses/gyrb/gyrb_prot_to_genome.tsv", sep="\t", quote=F)

#####