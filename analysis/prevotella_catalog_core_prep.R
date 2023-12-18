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


#### protein catalogue
rep_list_prevo = rep_list %>% filter(V2 %in% (tax2 %>% filter(genus=="g__Prevotella") %>% pull(final)))
rep_list_prevo$genome = lapply(rep_list_prevo$V1, function(x) gsub(".fa$|.fasta$|.fna$","", rev(strsplit(x, split="/")[[1]])[1])) %>% unlist
rep_list_prevo = rep_list_prevo %>% filter(!grepl("Manara", genome))

rep_list_outgroup = rep_list %>% filter(V2 == "GreatApes_cluster95_001501")
rep_list_outgroup$genome = lapply(rep_list_outgroup$V1, function(x) gsub(".fa$|.fasta$|.fna$","", rev(strsplit(x, split="/")[[1]])[1])) %>% unlist

rename_list = bind_rows(rep_list_prevo, rep_list_outgroup) %>% select(cluster=V2, genome) %>% mutate(rename = ifelse(cluster=="GreatApes_cluster95_001501", "Outgroup", gsub("_cluster95_","",cluster)))
write.table(rename_list, "/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/rename_list.tsv", sep="\t", row.names=F, quote=F)


cluster_annot = read_delim('/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/protein_catalogs/GreatApes.prot95_rep_seq.emapper.annotations',col_names=c('cluster_id','seed_ortholog','evalue','score','eggNOG_OGs','max_annot_lvl','COG_category','Description','Preferred_name','GOs','EC','KEGG_ko','KEGG_Pathway','KEGG_Module','KEGG_Reaction','KEGG_rclass','BRITE','KEGG_TC','CAZy','BiGG_Reaction','PFAMs'), delim='\t', quote='', num_threads=10,comment='')

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


prot_to_cluster = read_tsv('protein_catalogs/GreatApes.prot95_cluster.tsv', col_names=c('cluster_id','prot_id')) %>% mutate(contig_id = gsub('_[0-9]+$','',prot_id))
contig_to_sgb = read_tsv('GreatApes.contig_to_bin.tsv', col_names=c('mag_id','contig_id'))
mag_to_sgb = all_genomes %>% select(mag_id=genome, sample, GreatApes_95_final, host, host_group) 

cluster_to_mag  = prot_to_cluster %>% left_join(contig_to_sgb) %>% select(cluster_id, mag_id)
sgb_host_to_cluster = mag_to_sgb  %>% left_join(cluster_to_mag) %>% select(GreatApes_95_final, host_group, cluster_id) %>% distinct()

sgb_host_w_annot = sgb_host_to_cluster  %>% left_join(cluster_annot)

library(data.table)
#prot_to_genome = read_delim('protein_catalogs/GreatApes.prot_to_genome.tsv', delim='\t', quote='', num_threads=8,col_names=c('mag_id','prot_id'))
prot_to_cluster = read_delim('protein_catalogs/GreatApes.prot95_cluster.tsv', delim='\t', col_names=c('cluster_id','prot_id')) %>% 
	mutate(contig_id = gsub("_[0-9]+$","",prot_id))

cluster_to_sgb = all_genomes %>% select(mag_id = genome, GreatApes_95_final, host, host_group) %>%
	left_join(contig_to_sgb) %>%
	left_join(prot_to_cluster) %>% 
	select(GreatApes_95_final, mag_id, host_group, prot_id, cluster_id) %>% distinct()

cluster_to_sgb_tax = cluster_to_sgb %>% left_join(tax2, by=c("GreatApes_95_final"="final")) %>% 
	filter(genus %in% c("g__Prevotella"), host_group!="O")

### keep only annotations of prot clusters in prevo / paraprevo
annot_prevo = cluster_annot %>% filter(cluster_id %in% cluster_to_sgb_tax$cluster_id)

cluster_to_sgb_tax_annot = cluster_to_sgb_tax %>% left_join(annot_prevo %>% select(cluster_id, Preferred_name)) %>% filter(Preferred_name != "-")
n_per_host = cluster_to_sgb_tax_annot %>% select(GreatApes_95_final, host_group) %>% distinct %>% group_by(host_group) %>% summarise(n_total=n())
n_per_prot = cluster_to_sgb_tax_annot %>% select(GreatApes_95_final, host_group, Preferred_name) %>% distinct %>% group_by(host_group, Preferred_name) %>% summarise(n=n()) %>% left_join(n_per_host)  %>% mutate(rel_host = n/n_total)
n_per_all = cluster_to_sgb_tax_annot %>% select(GreatApes_95_final, Preferred_name) %>% distinct %>% group_by(Preferred_name) %>% summarise(n_across=n()) %>% left_join(n_per_prot)  %>% mutate(rel_across = n_across/(cluster_to_sgb_tax$GreatApes_95_final %>% unique %>% length))



##### Test set: genes found in 80% of Prevotella genomes (n=878) and in max 10% of genomes in > 2 copies (n = 626)

host_core_genome = n_per_prot %>% filter(host_group!="G", host_group!="O", rel_host>=.8)
host_core_genome = host_core_genome %>% left_join(cluster_to_sgb_tax_annot %>% filter(Preferred_name %in% host_core_genome$Preferred_name) %>% group_by(mag_id, Preferred_name) %>% summarise(n=n()) %>% group_by(Preferred_name) %>% summarise(mean_occ = mean(n)))


host_core_genome_testset = host_core_genome %>% filter(rel_host >= .8, mean_occ<1.05)
write.table(host_core_genome_testset, "/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/host_core_genome_testset.tsv", sep="\t", row.names=F, quote=F)

analysis_gene_set = cluster_to_sgb_tax_annot %>% filter(Preferred_name %in% host_core_genome_testset$Preferred_name)  %>% 
	left_join(rep_list_prevo %>% select(GreatApes_95_final=V2, rep_genome = genome))  %>% mutate(is_rep = mag_id == rep_genome) %>% 
	arrange(Preferred_name, GreatApes_95_final, -is_rep) 

analysis_gene_set_rep = analysis_gene_set %>% filter(!duplicated(paste(GreatApes_95_final, Preferred_name)))

write.table(analysis_gene_set_rep, "/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/analysis_gene_set_complete.tsv", sep="\t", row.names=F, quote=F)

#### Test set: genes with similar prevalence as cydA across all SGBs (~ 55.7 +/- 10%)

host_core_genome2 = n_per_all %>% filter(host_group=="H") %>% filter(between(rel_across, 0.557-0.1, 0.557+0.1))
host_core_genome2 = host_core_genome2 %>% left_join(cluster_to_sgb_tax_annot %>% filter(Preferred_name %in% host_core_genome2$Preferred_name) %>% group_by(mag_id, Preferred_name) %>% summarise(n=n()) %>% group_by(Preferred_name) %>% summarise(mean_occ = mean(n)))
write.table(host_core_genome2, "/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/host_core_genome_freqacross_testset.tsv", sep="\t", row.names=F, quote=F)



analysis_gene_set2 = cluster_to_sgb_tax_annot %>% filter(Preferred_name %in% host_core_genome2$Preferred_name)  %>% 
	left_join(rep_list_prevo %>% select(GreatApes_95_final=V2, rep_genome = genome))  %>% mutate(is_rep = mag_id == rep_genome) %>% 
	arrange(Preferred_name, GreatApes_95_final, -is_rep) 

analysis_gene_set_rep2 = analysis_gene_set2 %>% filter(!duplicated(paste(GreatApes_95_final, Preferred_name)))

write.table(analysis_gene_set_rep2, "/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/analysis_gene_set_complete_freqacross.tsv", sep="\t", row.names=F, quote=F)



#### Test set: 1000 random genes with at least 35% prevalence across the dataset

host_core_genome3 = n_per_all %>% filter(host_group=="H") %>% filter(rel_across >= 0.2)
host_core_genome3 = host_core_genome3 %>% left_join(cluster_to_sgb_tax_annot %>% filter(Preferred_name %in% host_core_genome3$Preferred_name) %>% group_by(mag_id, Preferred_name) %>% summarise(n=n()) %>% group_by(Preferred_name) %>% summarise(mean_occ = mean(n)))
write.table(host_core_genome3, "/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/host_core_genome_minprev.tsv", sep="\t", row.names=F, quote=F)

host_core_genome3_rand = host_core_genome3 %>% filter(mean_occ<1.05)
write.table(host_core_genome3_rand, "/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/host_core_genome_minprev_rand.tsv", sep="\t", row.names=F, quote=F)


analysis_gene_set3 = cluster_to_sgb_tax_annot %>% filter(Preferred_name %in% host_core_genome3$Preferred_name)  %>% 
	left_join(rep_list_prevo %>% select(GreatApes_95_final=V2, rep_genome = genome))  %>% mutate(is_rep = mag_id == rep_genome) %>% 
	arrange(Preferred_name, GreatApes_95_final, -is_rep) 

analysis_gene_set_rep3 = analysis_gene_set3 %>% filter(!duplicated(paste(GreatApes_95_final, Preferred_name)))

write.table(analysis_gene_set_rep3, "/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/analysis_gene_set_complete_minprev20.tsv", sep="\t", row.names=F, quote=F)








annot_og = cluster_annot %>% filter(cluster_id %in% cluster_to_sgb_tax_og$cluster_id)%>% filter(Preferred_name %in% host_core_genome_testset$Preferred_name)
cluster_to_sgb_tax_annot_og = cluster_to_sgb_tax_og %>% left_join(annot_og %>% select(cluster_id, Preferred_name)) %>% filter(Preferred_name != "-")

analysis_gene_set_og_part1 = cluster_to_sgb_tax_annot_og %>% filter(genus=="g__Paraprevotella", !duplicated(paste(Preferred_name)))
analysis_gene_set_og_part2 = cluster_to_sgb_tax_annot_og %>% filter(!Preferred_name %in% analysis_gene_set_og_part1$Preferred_name, genus=="g__Alloprevotella", !duplicated(paste(Preferred_name)))
analysis_gene_set_og_part3 = cluster_to_sgb_tax_annot_og %>% filter(!Preferred_name %in% c(analysis_gene_set_og_part1$Preferred_name,analysis_gene_set_og_part2$Preferred_name), genus=="g__Bacteroides", !duplicated(paste(Preferred_name)))

analysis_gene_set_og = bind_rows(analysis_gene_set_og_part1, analysis_gene_set_og_part2, analysis_gene_set_og_part3)

analysis_gene_set_complete = bind_rows(analysis_gene_set, analysis_gene_set_og) %>% arrange(Preferred_name, genus)

#####









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