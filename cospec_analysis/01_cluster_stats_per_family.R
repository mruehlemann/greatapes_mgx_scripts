library(tidyverse)

all_genomes = read.table("GreatApes.all_genomes.tsv", head=T, stringsAsFactors=F) %>% filter(!is.na(GreatApes_95_final))
#all_genomes_uhgg = read.table("../subgroups/UHGG_v2/genomes-all_metadata.tsv", head=T, stringsAsFactors=F, sep="\t")
#all_genomes_uhgg = all_genomes_uhgg %>% filter(Contamination < 10, Completeness>50)


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

colnames(tax2) = c("bin","kingdom","phylum","class","order","family","genus","species")

type_recode=data.frame(short=c("Ggorilla","Gberingei","HsA1","HsA2","HsG","Ptschw","Pttrog","Ptver","Ppan","HsUHGG"), 
	short2=c("Ggorilla","Gberingei","Hs","Hs","Hs","Ptschw","Pttrog","Ptver","Ppan","Hs"), 
	long=c("Gorilla_gorilla_gorilla","Gorilla_beringei","Homo_sapiens","Homo_sapiens","Homo_sapiens","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pan_paniscus","Homo_sapiens"),stringsAsFactors=F)

host = read.table("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/groupings_nozoo_campbell.txt", stringsAsFactors = F) %>% mutate(V2=ifelse(V2=="ZcambellGorilla","Ggorilla", V2), V2=ifelse(V2=="ZcambellPan","Pttrog", V2))
colnames(host) = c("sample","host")
host = bind_rows(host, data.frame(sample="MGYG00000",host="HsUHGG"))

all_genomes = all_genomes %>% left_join(tax2, by=c("GreatApes_95_final"="bin"))
all_genomes = all_genomes %>% mutate(sample=ifelse(grepl("^SRR", genome), substr(genome,1,10), substr(genome,1,9))) %>% left_join(host) %>% left_join(type_recode, by=c("host"="short")) %>% filter(sample!="ManaraSGB")

host_rep_genomes = all_genomes %>% arrange(-score) %>% distinct_at(.vars=c("long", "GreatApes_95_final"), .keep_all=T)

#### TREE BASED CLUSTERING

library(ape)
library(digest)
#library(phangorn)

allcand=list()

fams=list.files("fam_trees")
for(f in fams){
	tree_genomes_bac = read.tree(paste0("fam_trees/",f,"/",f,".faa.treefile"))
	tree_genomes_bac_rooted = root(tree_genomes_bac,"outgroup", resolve.root = TRUE)
	tree_genomes_bac_subtrees = subtrees(tree_genomes_bac_rooted)
	tree_genomes_bac_subtrees_matrix = lapply(seq_along(tree_genomes_bac_subtrees), function(x) data.frame(tips=tree_genomes_bac_subtrees[[x]]$tip.label, subtree=x)) %>% 
		do.call("rbind",.) %>% data.frame(stringsAsFactors=F) # %>% mutate(in_tree=1) %>% reshape2::dcast(subtree ~ tips, value.var="in_tree", fill=0) 

	tip_labs = (tree_genomes_bac %>% drop.tip(., "outgroup"))$tip.label
	ntips=length(tip_labs)

	
	subtree_stats = lapply(seq_along(tree_genomes_bac_subtrees), function(x){
	tree_sub = host_rep_genomes %>% filter(GreatApes_95_final %in% tree_genomes_bac_subtrees[[x]]$tip.label) %>% mutate(long=substr(host,1,2))
	this_stats=list(
			name=paste0(f,"_subtree_",x), id=x, n_genomes=nrow(tree_sub), n_clusters = length(unique(tree_sub$GreatApes_95_final)), tot_host = length(unique(tree_sub$short2)),
			tot_hostn2 = sum(table(tree_sub$short2) > 1), tot_hostgroups=length(unique(substr(tree_sub$host,1,1))), n_spe=length(unique(tree_sub$species)),
			n_gen=length(unique(tree_sub$genus)), n_fam=length(unique(tree_sub$family)))
	this_stats$tax = case_when(this_stats$n_spe == 1 ~ tree_sub$species[1], this_stats$n_gen == 1 ~ tree_sub$genus[1], this_stats$n_fam == 1 ~ tree_sub$family[1], TRUE ~ "none")
	return(this_stats)
	}) %>% do.call("bind_rows", .) %>% mutate(level="BAC120SUB") %>% mutate(edge_id=id+ntips)

	sst_bac=lapply(tip_labs, function(tl){
		print(tl)
		st_np = subtree_stats %>% filter(id %in% (tree_genomes_bac_subtrees_matrix %>% filter(tips==tl) %>% pull(subtree))) %>% filter(tot_host>=4)
		min_4 = st_np %>% arrange(n_clusters) %>% pull(id) %>% head(1)
		min_5 = st_np %>% filter(tot_host>=5) %>% arrange(n_clusters) %>% pull(id) %>% head(1)
		min_6 = st_np %>% filter(tot_host>=6)  %>% arrange(n_clusters) %>% pull(id) %>% head(1)
		min_7 = st_np %>% filter(tot_host>=7)  %>% arrange(n_clusters) %>% pull(id) %>% head(1)
		return(list(tl=tl, min_4=min_4, min_5=min_5, min_6=min_6,min_7=min_7))
	})%>% do.call("rbind",.) %>% data.frame

	subtree_stats$edge_id = subtree_stats$id + ntips
	min_4_excl = tree_genomes_bac_subtrees_matrix %>% left_join(subtree_stats, by=c("subtree"="id")) %>% filter(subtree %in% sst_bac$min_4) %>% arrange(n_genomes) %>% mutate(dup=duplicated(tips)) %>% filter(dup) %>% pull(subtree) %>% unique
	min_5_excl = tree_genomes_bac_subtrees_matrix %>% left_join(subtree_stats, by=c("subtree"="id")) %>% filter(subtree %in% sst_bac$min_5) %>% arrange(n_genomes) %>% mutate(dup=duplicated(tips)) %>% filter(dup) %>% pull(subtree) %>% unique
	min_6_excl = tree_genomes_bac_subtrees_matrix %>% left_join(subtree_stats, by=c("subtree"="id")) %>% filter(subtree %in% sst_bac$min_6) %>% arrange(n_genomes) %>% mutate(dup=duplicated(tips)) %>% filter(dup) %>% pull(subtree) %>% unique
	min_7_excl = tree_genomes_bac_subtrees_matrix %>% left_join(subtree_stats, by=c("subtree"="id")) %>% filter(subtree %in% sst_bac$min_7) %>% arrange(n_genomes) %>% mutate(dup=duplicated(tips)) %>% filter(dup) %>% pull(subtree) %>% unique

	min_4_incl = subtree_stats %>% filter(id %in% (sst_bac %>% pull(min_4) %>% unique ), !id %in% min_4_excl, tax!="none")
	min_5_incl = subtree_stats %>% filter(id %in% (sst_bac %>% pull(min_5) %>% unique ), !id %in% min_5_excl, tax!="none")
	min_6_incl = subtree_stats %>% filter(id %in% (sst_bac %>% pull(min_6) %>% unique ), !id %in% min_6_excl, tax!="none")
	min_7_incl = subtree_stats %>% filter(id %in% (sst_bac %>% pull(min_7) %>% unique ), !id %in% min_7_excl, tax!="none")

	candidates = bind_rows(min_4_incl, min_5_incl, min_6_incl, min_7_incl) %>% distinct(.) %>% arrange(-n_genomes)

	if(nrow(candidates)==0){next}

	candidates_long = apply(candidates, 1, function(x){
		data.frame(GreatApes_95_final=tree_genomes_bac_subtrees[[as.numeric(x[2])]]$tip.label) %>% mutate(group=x[1])
	}) %>% do.call("bind_rows",.)


	candidates_stats=lapply(candidates_long %>% pull(group) %>% unique, function(x){
		print(x)
		df=candidates_long %>% filter(group==x) %>% left_join(host_rep_genomes) %>% mutate(hostgroup = substr(host,1,1), host2=ifelse(grepl("Hs",host),"Hs",host)) %>% filter(!is.na(host))
		df2 = df %>% group_by(group) %>% summarize(n_genomes=n()) %>%
		left_join(df %>% group_by(group) %>% summarize(n_clusters=n())) %>%
		left_join(df %>% group_by(host2, group) %>% summarize(n=n()) %>% group_by(group) %>% summarize(tot_host=n())) %>%
		left_join(df %>% group_by(host2, group) %>% summarize(n=n()) %>% filter(n>=2) %>% group_by(group) %>% summarize(tot_hostn2=n())) %>%
		left_join(df %>% group_by(hostgroup, group) %>% summarize(n=n()) %>% group_by(group) %>% summarize(tot_hostgroups=n()))
		df3 =  table(factor(df$host, levels=unique(all_genomes$host))) %>% data.frame() %>% column_to_rownames("Var1") %>% t %>% data.frame(row.names=x) %>% rownames_to_column("group") %>% left_join(df2, .) %>% data.frame
		df3$n_phy = length(unique(df$phylum))
		df3$n_cla = length(unique(df$class))
		df3$n_ord = length(unique(df$order))
		df3$n_fam = length(unique(df$family))
		df3$n_gen = length(unique(df$genus))
		df3$n_spe = length(unique(df$species))
		df3$tax_cons = with(df3, case_when(n_spe == 1 ~ df$species[1], n_gen == 1 ~ df$genus[1],n_fam == 1 ~ df$family[1],n_ord == 1 ~ df$order[1],n_cla == 1 ~ df$class[1],n_phy == 1 ~ df$phylum[1], TRUE ~ df$kingdom[1]))
		df3$grouping = df$level[1]
		df3$hash = df$genome %>% sort %>% paste(., collapse=",") %>% digest(algo="md5")
		return(df3)
		}) %>% do.call("bind_rows",.)

		allcand[[f]] = list(candidates_stats=candidates_stats, candidates_long=candidates_long)
}


#####

groups_all = lapply(allcand, function(ac) return(ac$candidates_long)) %>% do.call("bind_rows",.) %>% filter(GreatApes_95_final!="outgroup")



library(digest)

groups_all_stats=lapply(groups_all %>% pull(group) %>% unique, function(x){
	print(x)
	df=groups_all %>% filter(group==x) %>% left_join(host_rep_genomes) %>% mutate(hostgroup = substr(host,1,1), host2=ifelse(grepl("Hs",host),"Hs",host)) %>% filter(!is.na(host))
	df2 = df %>% group_by(group) %>% summarize(n_genomes=n()) %>%
	left_join(df %>% group_by(group) %>% summarize(n_clusters=n())) %>%
	left_join(df %>% group_by(host2, group) %>% summarize(n=n()) %>% group_by(group) %>% summarize(tot_host=n())) %>%
	left_join(df %>% group_by(host2, group) %>% summarize(n=n()) %>% filter(n>=2) %>% group_by(group) %>% summarize(tot_hostn2=n())) %>%
	left_join(df %>% group_by(hostgroup, group) %>% summarize(n=n()) %>% group_by(group) %>% summarize(tot_hostgroups=n()))
	df3 =  table(factor(df$host, levels=unique(all_genomes$host))) %>% data.frame() %>% column_to_rownames("Var1") %>% t %>% data.frame(row.names=x) %>% rownames_to_column("group") %>% left_join(df2, .) %>% data.frame
	df3$n_phy = length(unique(df$phylum))
	df3$n_cla = length(unique(df$class))
	df3$n_ord = length(unique(df$order))
	df3$n_fam = length(unique(df$family))
	df3$n_gen = length(unique(df$genus))
	df3$n_spe = length(unique(df$species))
	df3$tax_cons = with(df3, case_when(n_spe == 1 ~ df$species[1], n_gen == 1 ~ df$genus[1],n_fam == 1 ~ df$family[1],n_ord == 1 ~ df$order[1],n_cla == 1 ~ df$class[1],n_phy == 1 ~ df$phylum[1], TRUE ~ df$kingdom[1]))
	df3$grouping = df$level[1]
	df3$hash = df$genome %>% sort %>% paste(., collapse=",") %>% digest(algo="md5")
	return(df3)
	}) %>% do.call("bind_rows",.)


groups_all_stats_candidates = groups_all_stats %>% filter(!duplicated(hash), tot_host >= 4, tot_hostgroups >= 2, n_fam==1) %>% arrange(-n_genomes)
groups_all_candidates = groups_all %>% filter(group %in% groups_all_stats_candidates$group) %>% left_join(host_rep_genomes) %>% mutate(hostgroup = substr(host,1,1), host2=ifelse(grepl("Hs",host),"Hs",host))

dir.create("../cospec_analysis", recursive=T)

write.table(host_rep_genomes, "../cospec_analysis/genomes_included.tsv", sep="\t", row.names=F, quote=F)

saveRDS(groups_all,"../cospec_analysis/groups_all.Rds")
saveRDS(groups_all_candidates,"../cospec_analysis/groups_all_candidates.Rds")
saveRDS(groups_all_stats, "../cospec_analysis/groups_all_stats.Rds")
saveRDS(groups_all_stats_candidates, "../cospec_analysis/groups_all_stats_candidates.Rds")
write.table(groups_all_stats_candidates, "../cospec_analysis/groups_all_stats_candidates.tsv", sep="\t", row.names=F, quote=F)
write.table(groups_all_candidates, "../cospec_analysis/groups_all_candidates.tsv", sep="\t", row.names=F, quote=F)
