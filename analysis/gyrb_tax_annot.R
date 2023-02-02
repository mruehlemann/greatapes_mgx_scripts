setwd("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups")
library(tidyverse)

rep_list=read.table("GreatApes.representatives_list", head=F, stringsAsFactors=F)
manara_hits = read.table("GreatApes.representatives.manara_mash.tsv",head=F, stringsAsFactors=F) %>%
  left_join(rep_list %>% mutate(V1=gsub("/home/sukmb276/Isilon","/work_ifs/sukmb276", V1)), by=c("V2"="V1")) %>% filter(V3<=.065) %>%
  select(V2.y)

ag=read.table("GreatApes.all_genomes.tsv",head=T, stringsAsFactors=F)
ag_clean = ag %>% filter(!is.na(GreatApes_95_final)) %>% mutate(host=substr(cluster_97_final,1,3))

sgb_stats = ag_clean %>% select(GreatApes_95_final, host) %>% reshape2::dcast(GreatApes_95_final ~ host, data=.) %>%
  mutate(Manara=ifelse(GreatApes_95_final %in% manara_hits$V2.y,1,0)) %>% (function(df) data.frame(df, ngroup=rowSums(df[,-1]>0))) %>%
  filter(!(MGY>0 & ngroup==1)) %>% mutate(novel=ifelse(Manara==0 & MGY==0, 1,0)) %>% mutate(H=HsA+HsG+MGY, G=Gbe+Ggo, P=Ppa+Pts+Ptt+Ptv, ngen=(H>0)+(P>0)+(G>0)) %>%
  column_to_rownames("GreatApes_95_final")

group_novelty = lapply(colnames(sgb_stats), function(x) table(this=sgb_stats[,x]>0,novel=sgb_stats$novel==1) %>% reshape2::melt(.) %>% mutate(group=x)) %>% do.call("bind_rows", .) %>% filter(this, !group %in% c("Manara","ngroup","novel","MGY") ) %>% select(group,novel, value)


tax=read.table("GreatApes.cluster_final_tax.tsv", head=T, stringsAsFactors=F, sep="\t")
colnames(tax) = c("GreatApes_95_final","classification")

tax2 = sapply(paste0(tax$GreatApes_95_final,";",tax$classification), function(x){
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

#### protein catalogue
library(data.table)
prot_to_genome = read_delim('prot_to_genome.tsv', delim='\t', quote='', num_threads=10,col_names=c('mag_id','prot_id'))
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
gyrb_clusters = cluster_annot %>% filter(V9 == "gyrB") %>% pull(V1)
### get all gyrB proteins; n=7916
gyrb_prot = prot_to_cluster %>% filter(cluster %in% gyrb_clusters) %>% pull(protid)
### get all gyrB to genome mappings
gyrb_prot_to_genome = prot_to_genome %>% filter(protid %in% gyrb_prot) %>% left_join(ag_clean %>% select(genome, host, GreatApes_95_final)) %>% left_join(tax2)

dir.create("../analyses/gyrb", recursive = T)
write.table(gyrb_prot_to_genome, "../analyses/cyda/cyda_prot_to_genome.tsv", sep="\t", quote=F)

#####