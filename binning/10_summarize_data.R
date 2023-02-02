setwd("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups")
.libPaths("/oldhome/sukmb276/R/x86_64-pc-linux-gnu-library/3.5/")
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


### Cospec results
cosp_cand=read.table("../cospec_analysis/groups_all_candidates.tsv",head=T, stringsAsFactors=F)
cosp_cand_group=read.table("../cospec_analysis/groups_all_candidates.tsv",head=T, stringsAsFactors=F)
cospec_res_smry = readRDS("../cospec_analysis/cospec_final.Rds")
cospec_res = lapply(cospec_res_smry$results, function(x) x$res) %>% do.call("bind_rows",.)

#### protein catalogue
library(data.table)
prot_to_genome=fread("protein_catalogs/GreatApes.protid_to_genomes.tsv", sep="\t", col.names=c("genome","protid"))
prot_to_cluster=fread("protein_catalogs/GreatApes.prot95_cluster.tsv", sep="\t", col.names=c("cluster","protid"))
cluster_annot = fread('grep -v "^#" protein_catalogs/GreatApes.prot95_rep_seq.emapper.annotations',sep="\t")
cluster_to_sgb = ag_clean %>% select(genome, GreatApes_95_final, host) %>% left_join(prot_to_genome) %>% left_join(prot_to_cluster) %>% mutate(host_group=substr(host,1,1)) %>% select(GreatApes_95_final, host_group, cluster) %>% distinct()

### parse eggnogs
cluster_eggnogs=apply(cluster_annot[,c(1,5)], 1, function(x){ogs=gsub("@.+$","",strsplit(x[2],split=",")[[1]]); return(data.frame(ogs=unique(ogs), cluster=x[1], stringsAsFactors=F))}) %>% do.call("bind_rows", .)
saveRDS(cluster_eggnogs, "protein_catalogs/GreatApes.cluster_eggnogs.Rds")

eggnog_to_sgb = cluster_to_sgb %>% left_join(cluster_eggnogs) %>% select(-cluster) %>% filter(!is.na(ogs)) %>% distinct() %>% left_join(tax2)
saveRDS(eggnog_to_sgb, "protein_catalogs/GreatApes.eggnog_to_sgb.Rds")
#eggnog_to_sgb=readRDS("protein_catalogs/GreatApes.eggnog_to_sgb.Rds")

### parse cog categories
cluster_cogcat=apply(cluster_annot[,c(1,7)] %>% filter(V7!="-"), 1, function(x){ogs=strsplit(x[2],split="")[[1]]; return(data.frame(ogs=unique(ogs), cluster=x[1], stringsAsFactors=F))}) %>% do.call("bind_rows", .)
saveRDS(cluster_cogcat, "protein_catalogs/GreatApes.cluster_cogcat.Rds")

cogcat_to_sgb = cluster_to_sgb %>% left_join(cluster_cogcat) %>% select(-cluster) %>% filter(!is.na(ogs)) %>% distinct() %>% left_join(tax2)
saveRDS(cogcat_to_sgb, "protein_catalogs/GreatApes.cogcat_to_sgb.Rds")

### parse EC
cluster_ecs = apply(cluster_annot[,c(1,11)] %>% filter(V11!="-"), 1, function(x){ogs=strsplit(x[2],split=",")[[1]]; return(data.frame(ogs=unique(ogs), cluster=x[1], stringsAsFactors=F))}) %>% do.call("bind_rows", .)
saveRDS(cluster_ecs, "protein_catalogs/GreatApes.cluster_ecs.Rds")

ecs_to_sgb = cluster_to_sgb %>% left_join(cluster_ecs) %>% select(-cluster) %>% filter(!is.na(ogs)) %>% distinct() %>% left_join(tax2)
saveRDS(ecs_to_sgb, "protein_catalogs/GreatApes.ecs_to_sgb.Rds")

### parse keggko
cluster_keggko = apply(cluster_annot[,c(1,12)] %>% filter(V12!="-"), 1, function(x){ogs=strsplit(x[2],split=",")[[1]]; return(data.frame(ogs=unique(ogs), cluster=x[1], stringsAsFactors=F))}) %>% do.call("bind_rows", .)
saveRDS(cluster_keggko, "protein_catalogs/GreatApes.cluster_keggko.Rds")

keggko_to_sgb = cluster_to_sgb %>% left_join(cluster_keggko) %>% select(-cluster) %>% filter(!is.na(ogs)) %>% distinct() %>% left_join(tax2)
saveRDS(keggko_to_sgb, "protein_catalogs/GreatApes.keggko_to_sgb.Rds")

### parse keggpwy
cluster_keggpwy = apply(cluster_annot[,c(1,13)] %>% filter(V13!="-"), 1, function(x){ogs=strsplit(x[2],split=",")[[1]]; return(data.frame(ogs=unique(ogs), cluster=x[1], stringsAsFactors=F))}) %>% do.call("bind_rows", .)
saveRDS(cluster_keggpwy, "protein_catalogs/GreatApes.cluster_keggpwy.Rds")

keggpwy_to_sgb = cluster_to_sgb %>% left_join(cluster_keggpwy) %>% select(-cluster) %>% filter(!is.na(ogs)) %>% distinct() %>% left_join(tax2)
saveRDS(keggpwy_to_sgb, "protein_catalogs/GreatApes.keggpwy_to_sgb.Rds")

### parse keggmod
cluster_keggmod = apply(cluster_annot[,c(1,14)] %>% filter(V14!="-"), 1, function(x){ogs=strsplit(x[2],split=",")[[1]]; return(data.frame(ogs=unique(ogs), cluster=x[1], stringsAsFactors=F))}) %>% do.call("bind_rows", .)
saveRDS(cluster_keggmod, "protein_catalogs/GreatApes.cluster_keggmod.Rds")

keggmod_to_sgb = cluster_to_sgb %>% left_join(cluster_keggmod) %>% select(-cluster) %>% filter(!is.na(ogs)) %>% distinct() %>% left_join(tax2)
saveRDS(keggmod_to_sgb, "protein_catalogs/GreatApes.keggmod_to_sgb.Rds")

### parse keggrc
cluster_keggrc= apply(cluster_annot[,c(1,16)] %>% filter(V16!="-"), 1, function(x){ogs=strsplit(x[2],split=",")[[1]]; return(data.frame(ogs=unique(ogs), cluster=x[1], stringsAsFactors=F))}) %>% do.call("bind_rows", .)
saveRDS(cluster_keggrc, "protein_catalogs/GreatApes.cluster_keggrc.Rds")

keggrc_to_sgb = cluster_to_sgb %>% left_join(cluster_keggrc) %>% select(-cluster) %>% filter(!is.na(ogs)) %>% distinct() %>% left_join(tax2)
saveRDS(keggrc_to_sgb, "protein_catalogs/GreatApes.keggrc_to_sgb.Rds")

### parse cazy
cluster_cazy= apply(cluster_annot[,c(1,19)] %>% filter(V19!="-"), 1, function(x){ogs=strsplit(x[2],split=",")[[1]]; return(data.frame(ogs=unique(ogs), cluster=x[1], stringsAsFactors=F))}) %>% do.call("bind_rows", .)
saveRDS(cluster_cazy, "protein_catalogs/GreatApes.cluster_cazy.Rds")

cazy_to_sgb = cluster_to_sgb %>% left_join(cluster_cazy) %>% select(-cluster) %>% filter(!is.na(ogs)) %>% distinct() %>% left_join(tax2)
saveRDS(cazy_to_sgb, "protein_catalogs/GreatApes.cazy_to_sgb.Rds")

### parse pfam
cluster_pfam= apply(cluster_annot[,c(1,21)] %>% filter(V21!="-"), 1, function(x){ogs=strsplit(x[2],split=",")[[1]]; return(data.frame(ogs=unique(ogs), cluster=x[1], stringsAsFactors=F))}) %>% do.call("bind_rows", .)
saveRDS(cluster_pfam, "protein_catalogs/GreatApes.cluster_pfam.Rds")

pfam_to_sgb = cluster_to_sgb %>% left_join(cluster_pfam) %>% select(-cluster) %>% filter(!is.na(ogs)) %>% distinct() %>% left_join(tax2)
saveRDS(pfam_to_sgb, "protein_catalogs/GreatApes.pfam_to_sgb.Rds")


###
cosp_clusters = cosp_cand %>% filter(group %in% c(cospec_res %>% filter(hommola_p_b01<0.05) %>% pull(clade))) %>% pull(1) %>% unique()
cosp_candidates = cosp_cand %>% pull(1) %>% unique()


this_to_sgb = readRDS("protein_catalogs/GreatApes.keggko_to_sgb.Rds")
cands= this_to_sgb %>% filter(host_group!="G", !is.na(genus)) %>% select(GreatApes_95_final,host_group,genus) %>% distinct %>% group_by(host_group, genus) %>% summarize(n=n()) %>% reshape2::dcast(genus ~ host_group, value.var="n", fill=0) %>% filter(H>=3, P>=3)


cands_res = lapply(cands$genus, function(tax){print(tax); tax_sub = this_to_sgb %>% filter(genus==tax, host_group!="G") %>% mutate(presence=1) %>% reshape2::dcast(GreatApes_95_final+host_group ~ ogs, value.var="presence", fill=0)
tax_res=lapply(colnames(tax_sub %>% select(-GreatApes_95_final, -host_group))[between(colMeans(tax_sub[,-1][,-1]),.1,.9)], function(x) return(list(ogs=x, H_prev=mean(tax_sub[tax_sub$host_group=="H",x]), P_prev=mean(tax_sub[tax_sub$host_group=="P",x]), pval=fisher.test(table(tax_sub$host_group, tax_sub[,x]))$p.value))) %>% do.call("bind_rows",.) %>% mutate(tax=tax)
tax_res %>% arrange(pval) %>% return()}) %>% do.call("bind_rows",.)


KO_TO_MODULE=read.table("kegg_ko_to_module/ALL_KEGG_MODULES.tsv",head=F,stringsAsFactors=F)
MODULE_NAMES=read.table("kegg_ko_to_module/ALL_KEGG_MODULES_NAMES.tsv",head=F,stringsAsFactors=F,sep="\t")

cands_res %>% left_join(tax2 %>% select(tax=genus, family)  %>% group_by(tax, family) %>% summarize(n=n()))  %>% mutate(effect=log2((H_prev+0.01)/(P_prev+0.01))) %>% mutate(Zscore=-qnorm(pval/2)*sign(effect)) %>% group_by(ogs, family) %>% summarize(Z_mean=mean(Zscore), p_meta=2*pnorm(-abs(Z_mean))) %>% arrange(p_meta)

cands_res %>% left_join(tax2 %>% select(tax=genus, family)  %>% group_by(tax, family) %>% summarize(n=n()))  %>%
  mutate(effect=log2((H_prev+0.01)/(P_prev+0.01))) %>%
  mutate(Zscore=-qnorm(pval/2)*sign(effect)) %>%
  (function(df) df %>% left_join(df %>% group_by(ogs, family) %>% summarize(n_tot=sum(n)))) %>%
  mutate(w=sqrt(n), Z_w=Zscore*w) %>% group_by(ogs, family) %>% summarize(Z_sum = sum(Z_w), w_sum=sum(w^2)) %>%
  mutate(Z_meta=Z_sum/sqrt(w_sum), p_meta=2*pnorm(-abs(Z_meta))) %>% arrange(p_meta) %>% filter(p_meta < 0.05) %>% data.frame


cands_res %>% left_join(KO_TO_MODULE, by=c("ogs"="V2")) %>% filter(!is.na(V1)) %>%
  mutate(effect=log2((H_prev+0.01)/(P_prev+0.01))) %>%
  mutate(Zscore=-qnorm(pval/2)*sign(effect)) %>%
  (function(df) df %>% left_join(df %>% group_by(ogs, V1))) %>%
  mutate(w=1, Z_w=Zscore*w) %>% group_by(tax, V1) %>% summarize(Z_sum = sum(Z_w), w_sum=sum(w^2)) %>%
  mutate(Z_meta=Z_sum/sqrt(w_sum), p_meta=2*pnorm(-abs(Z_meta))) %>% arrange(p_meta) %>% left_join(MODULE_NAMES)

#####
this_to_sgb_sub = this_to_sgb %>% filter(GreatApes_95_final %in% cosp_candidates)
cands2 = this_to_sgb_sub %>% mutate(cosp=ifelse(GreatApes_95_final %in% cosp_clusters,"yes","no")) %>% select(GreatApes_95_final,cosp,family) %>% distinct %>% group_by(cosp, family) %>% summarize(n=n()) %>% reshape2::dcast(family ~ cosp, value.var="n", fill=0) %>% filter(yes>=3, no>=3)

cands2_res = lapply(cands2$family, function(tax){print(tax); tax_sub = this_to_sgb_sub %>% mutate(cosp=ifelse(GreatApes_95_final %in% cosp_clusters,"yes","no")) %>% filter(family==tax) %>% select(GreatApes_95_final, cosp, ogs) %>% distinct %>% mutate(presence=1) %>% reshape2::dcast(GreatApes_95_final+cosp ~ ogs, value.var="presence", fill=0)
tax_res=lapply(colnames(tax_sub %>% select(-GreatApes_95_final, -cosp))[between(colMeans(tax_sub[,-1][,-1]),.1,.9)], function(x) return(list(ogs=x, cosp_prev=mean(tax_sub[tax_sub$cosp=="yes",x]), noncosp_prev=mean(tax_sub[tax_sub$cosp=="no",x]), pval=fisher.test(table(tax_sub$cosp, tax_sub[,x]))$p.value))) %>% do.call("bind_rows",.) %>% mutate(tax=tax)
tax_res %>% arrange(pval) %>% return()}) %>% do.call("bind_rows",.)


cands2_res %>% left_join(KO_TO_MODULE, by=c("ogs"="V2")) %>% filter(!is.na(V1)) %>%
mutate(effect=log2((cosp_prev+0.01)/(noncosp_prev+0.01))) %>%
mutate(Zscore=-qnorm(pval/2)*sign(effect)) %>%
(function(df) df %>% left_join(df %>% group_by(ogs, V1))) %>% #data.frame %>% arrange(pval)%>% head(20)
mutate(w=2, Z_w=Zscore*w) %>% group_by(tax, V1) %>% summarize(Z_sum = sum(Z_w), w_sum=sum(w^2)) %>%
mutate(Z_meta=Z_sum/(w_sum), p_meta=2*pnorm(-abs(Z_meta))) %>% arrange(p_meta) %>% left_join(MODULE_NAMES) %>% data.frame() %>% head(20)

cands2_res %>% left_join(KO_TO_MODULE, by=c("ogs"="V2"))  %>%
mutate(effect=log2((cosp_prev+0.01)/(noncosp_prev+0.01))) %>%
mutate(Zscore=-qnorm(pval/2)*sign(effect)) %>%
#(function(df) df %>% left_join(df %>% group_by(ogs, V1))) %>% #data.frame %>% arrange(pval)%>% head(20)
mutate(w=1, Z_w=Zscore*w) %>% group_by(ogs) %>% summarize(Z_sum = sum(Z_w), w_sum=sum(w^2)) %>%
mutate(Z_meta=Z_sum/(w_sum), p_meta=2*pnorm(-abs(Z_meta))) %>% arrange(p_meta) %>% filter(w_sum>3)

 left_join(MODULE_NAMES) %>% data.frame() %>% head(20)

####

cands3= this_to_sgb %>% filter( GreatApes_95_final %in% cosp_clusters, !is.na(family)) %>% select(GreatApes_95_final,host_group,family) %>% distinct %>% group_by(host_group, family) %>% summarize(n=n()) %>% reshape2::dcast(family ~ host_group, value.var="n", fill=0) %>% filter(H>=3, P>=3)


cands_res = lapply(cands$family, function(tax){print(tax); tax_sub = this_to_sgb %>% filter(family==tax, host_group!="G") %>% mutate(presence=1) %>% reshape2::dcast(GreatApes_95_final+host_group ~ ogs, value.var="presence", fill=0)
tax_res=lapply(colnames(tax_sub %>% select(-GreatApes_95_final, -host_group))[between(colMeans(tax_sub[,-1][,-1]),.1,.9)], function(x) return(list(ogs=x, H_prev=mean(tax_sub[tax_sub$host_group=="H",x]), P_prev=mean(tax_sub[tax_sub$host_group=="P",x]), pval=fisher.test(table(tax_sub$host_group, tax_sub[,x]))$p.value))) %>% do.call("bind_rows",.) %>% mutate(tax=tax)
tax_res %>% arrange(pval) %>% return()}) %>% do.call("bind_rows",.)
