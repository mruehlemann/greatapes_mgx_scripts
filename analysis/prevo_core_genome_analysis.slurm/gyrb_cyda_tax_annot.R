setwd("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups")
library(tidyverse)

rep_list=read.table("GreatApes.representatives_list", head=F, stringsAsFactors=F)

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

colnames(tax2) = c("final","kingdom","phylum","class","order","family","genus","species")

type_recode=data.frame(short=c("Ggorilla","Gberingei","HsA1","HsA2","HsG","Ptschw","Pttrog","Ptver","Ppan","HsUHGG","HsD"), 
	short2=c("Ggorilla","Gberingei","Hs","Hs","Hs","Ptschw","Pttrog","Ptver","Ppan","Hs","Hs"), 
	long=c("Gorilla_gorilla_gorilla","Gorilla_beringei","Homo_sapiens","Homo_sapiens","Homo_sapiens","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pan_paniscus","Homo_sapiens","Homo_sapiens"),stringsAsFactors=F)


host = read_tsv("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/groupings_nozoo_campbell_dk.txt", col_names=c('sample','host')) %>% 
    mutate(host=ifelse(host=="ZcambellGorilla","Ggorilla", host), host=ifelse(host=="ZcambellPan","Pttrog", host)) %>%
    mutate(host_group=substr(host,1,1)) %>% 
    bind_rows(., data.frame(sample='MGYG00000', host='HsUHGG',host_group='H')) %>% 
    bind_rows(., data.frame(sample='Manara', host='Other',host_group='O')) 

all_genomes = read_delim('GreatApes.all_genomes.tsv', delim='\t',col_names=T) %>% 
    filter(!is.na(final)) %>% mutate(GreatApes_95_final=final) %>% 
    filter(!(grepl('MGYG', genome) & genome != cluster_99_final)) %>%
    mutate(sample = case_when(grepl('MGYG', genome) ~ 'MGYG00000', grepl('Manara', genome) ~ 'Manara', T ~ gsub('_cleanbin_[0-9]+','', genome))) %>%
    left_join(host) %>% filter(!duplicated(genome))


cluster_annot = read_delim('/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/protein_catalogs/GreatApes.prot95_rep_seq.emapper.annotations',col_names=c('cluster_id','seed_ortholog','evalue','score','eggNOG_OGs','max_annot_lvl','COG_category','Description','Preferred_name','GOs','EC','KEGG_ko','KEGG_Pathway','KEGG_Module','KEGG_Reaction','KEGG_rclass','BRITE','KEGG_TC','CAZy','BiGG_Reaction','PFAMs'), delim='\t', quote='', num_threads=10,comment='')

prot_to_cluster = read_tsv('protein_catalogs/GreatApes.prot95_cluster.tsv', col_names=c('cluster_id','prot_id')) %>% mutate(contig_id = gsub('_[0-9]+$','',prot_id))
contig_to_sgb = read_tsv('GreatApes.contig_to_bin.tsv', col_names=c('mag_id','contig_id'))
mag_to_sgb = all_genomes %>% select(mag_id=genome, sample, GreatApes_95_final, host, host_group) 

cluster_to_mag  = prot_to_cluster %>% left_join(contig_to_sgb) %>% select(cluster_id, mag_id)
sgb_host_to_cluster = mag_to_sgb  %>% left_join(cluster_to_mag) %>% select(GreatApes_95_final, host_group, cluster_id) %>% distinct()

#sgb_host_w_annot = sgb_host_to_cluster  %>% left_join(cluster_annot)

library(data.table)

cluster_to_sgb = all_genomes %>% select(mag_id = genome, GreatApes_95_final, host, host_group) %>%
	left_join(contig_to_sgb) %>%
	left_join(prot_to_cluster) %>% 
	select(GreatApes_95_final, mag_id, host_group, prot_id, cluster_id) %>% distinct() 


### get all gyrB annotated protein clusters; n=1329
gyrb_clusters = cluster_annot %>% filter(Preferred_name == "gyrB") %>% pull(cluster_id)

### get all gyrB proteins; n=7916
gyrb_prot = prot_to_cluster %>% filter(cluster_id %in% gyrb_clusters) %>% pull(prot_id)
### get all gyrB to genome mappings
gyrb_prot_to_genome = prot_to_genome %>% filter(prot_id %in% gyrb_prot) %>% left_join(all_genomes %>% select(mag_id=genome, host, GreatApes_95_final)) %>% left_join(tax2)

gyrb_prot_to_genome = cluster_to_sgb %>% filter(cluster_id %in% gyrb_clusters) %>% left_join(all_genomes %>% select(mag_id=genome, host, GreatApes_95_final)) %>% left_join(tax2, by=c("GreatApes_95_final"="final"))

dir.create("../analyses/gyrb", recursive = T)
write.table(gyrb_prot_to_genome, "../analyses/gyrb/gyrb_prot_to_genome.tsv", sep="\t", quote=F)

#####