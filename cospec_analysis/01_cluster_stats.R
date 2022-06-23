library(tidyverse)

all_genomes = read.table("GreatApes.all_genomes.tsv", head=T, stringsAsFactors=F) %>% filter(!is.na(GreatApes_95_final))
tax=read.table("GreatApes.cleanbins_tax.tsv", head=T, stringsAsFactors=F, sep="\t")
colnames(tax) = c("bin","classification")

tax2 = sapply(tax$classification, function(x) strsplit(x, split=";")[[1]]) %>%
	apply(., 2, function(x){x_low=x[max(grep("__$",x, invert=T))]; x[grepl("__$",x)] = paste0(x[grepl("__$",x)],"unclassified_", x_low); return(x)}) %>%
	t %>% data.frame(row.names=NULL, stringsAsFactors=F, bin=tax$bin, .)
colnames(tax2) = c("bin","kingdom","phylum","class","order","family","genus","species")
tax2$species=gsub(" ","_",tax2$species)

all_genomes = all_genomes %>% left_join(tax2, by=c("genome"="bin"))
all_genomes$host = sapply(all_genomes$cluster_97_final, function(x) strsplit(x, split="_")[[1]][1])

group_stats=lapply(c("family","genus","species","GreatApes_95_final"), function(x){
	all_genomes %>% select(host, group = one_of(x)) %>% filter(!grepl("unclassified", group)) %>% mutate(hostgroup = substr(host,1,1)) %>%
	(function(df){ df %>% group_by(group) %>% summarize(n_genomes=n()) %>%
	left_join(df %>% group_by(host, group) %>% summarize(n=n()) %>% group_by(group) %>% summarize(tot_host=n())) %>%
	left_join(df %>% group_by(host, group) %>% summarize(n=n()) %>% filter(n>=2) %>% group_by(group) %>% summarize(tot_hostn2=n())) %>%
	left_join(df %>% group_by(hostgroup, group) %>% summarize(n=n()) %>% group_by(group) %>% summarize(tot_hostgroups=n()))}) %>% mutate(level=x)
}
) %>% do.call("bind_rows",.)


taxgroups_long = apply(group_stats,1, function(x) all_genomes %>% filter(get(x[6]) == x[1]) %>% select(genome) %>% mutate(level=x[6], group=x[1])) %>% do.call("bind_rows",.) %>% select(group, level, genome)

#### TREE BASED CLUSTERING

library(ape)
#library(phangorn)

tree_genomes_bac = read.tree("gtdbtk_r207v2_out/align/bac120_infer_tree/gtdbtk.unrooted.tree")
tree_genomes_bac_subtrees = subtrees(tree_genomes_bac)

subtree_stats = lapply(seq_along(tree_genomes_bac_subtrees), function(x){
tree_sub = all_genomes %>% filter(GreatApes_95_final %in% tree_genomes_bac_subtrees[[x]]$tip.label)
this_stats=list(
		name=paste0("BAC120SUB_",x), id=x, n_genomes=nrow(tree_sub), n_clusters = length(unique(tree_sub$GreatApes_95_final)), tot_host = length(unique(tree_sub$host)),
		tot_hostn2 = sum(table(tree_sub$host) > 1), tot_hostgroups=length(unique(substr(tree_sub$host,1,1))), n_spe=length(unique(tree_sub$species)),
		n_gen=length(unique(tree_sub$genus)), n_fam=length(unique(tree_sub$family)))
this_stats$tax = case_when(this_stats$n_spe == 1 ~ tree_sub$species[1], this_stats$n_gen == 1 ~ tree_sub$genus[1], this_stats$n_fam == 1 ~ tree_sub$family[1], TRUE ~ "none")
return(this_stats)
}) %>% do.call("bind_rows", .) %>% mutate(level="BAC120SUB")

treegroupsbac_long = lapply(subtree_stats$id, function(x) all_genomes %>% filter(GreatApes_95_final %in% tree_genomes_bac_subtrees[[x]]$tip.label) %>% select(genome) %>% mutate(level="BAC120SUB", group=paste0("BAC120SUB_",x))) %>% do.call("bind_rows",.) %>% select(group, level, genome)


tree_genomes_ar = read.tree("gtdbtk_r207v2_out/align/ar53_infer_tree/gtdbtk.unrooted.tree")
tree_genomes_ar_subtrees = subtrees(tree_genomes_ar)

subtree_stats_ar = lapply(seq_along(tree_genomes_ar_subtrees), function(x){
tree_sub = all_genomes %>% filter(GreatApes_95_final %in% tree_genomes_bac_subtrees[[x]]$tip.label)
this_stats=list(
		name=paste0("AR53SUB_",x), id=x, n=nrow(tree_sub), n_clusters = length(unique(tree_sub$GreatApes_95_final)), tot_host = length(unique(tree_sub$host)),
		tot_hostn2 = sum(table(tree_sub$host) > 1), tot_hostgroups=length(unique(substr(tree_sub$host,1,1))), n_spe=length(unique(tree_sub$species)), n_gen=length(unique(tree_sub$genus)), n_fam=length(unique(tree_sub$family)))
this_stats$tax = case_when(this_stats$n_spe == 1 ~ tree_sub$species[1], this_stats$n_gen == 1 ~ tree_sub$genus[1], this_stats$n_fam == 1 ~ tree_sub$family[1], TRUE ~ "none")
return(this_stats)
}) %>% do.call("bind_rows", .) %>% mutate(level="AR53SUB")

treegroupsar_long = lapply(subtree_stats_ar$id, function(x) all_genomes %>% filter(GreatApes_95_final %in% tree_genomes_ar_subtrees[[x]]$tip.label) %>% select(genome) %>% mutate(level="AR53SUB", group=paste0("AR53SUB_",x))) %>% do.call("bind_rows",.) %>% select(group, level, genome)

#####
groups_all = bind_rows(taxgroups_long, treegroupsbac_long,treegroupsar_long)

### L1_cleanbin_000073 is an MX02 archaean wrongly clustered with bacteria by dRep which causes trouble in the marker alignment
groups_all[groups_all$level=="BAC120SUB" & groups_all$genome=="H07606-L1_cleanbin_000073","group"]="NoCluster"


library(digest)

groups_all_stats=lapply(unique(groups_all$group), function(x){
	print(x)
	df=groups_all %>% filter(group==x) %>% left_join(all_genomes) %>% mutate(hostgroup = substr(host,1,1), host=ifelse(grepl("Hs",host),"Hs",host))
	df2 = df %>% group_by(group) %>% summarize(n_genomes=n()) %>%
	left_join(df %>% group_by(group) %>% summarize(n_clusters=n())) %>%
	left_join(df %>% group_by(host, group) %>% summarize(n=n()) %>% group_by(group) %>% summarize(tot_host=n())) %>%
	left_join(df %>% group_by(host, group) %>% summarize(n=n()) %>% filter(n>=2) %>% group_by(group) %>% summarize(tot_hostn2=n())) %>%
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


groups_all_stats_candidates = groups_all_stats %>% arrange(grouping) %>% map_df(rev) %>% filter(!duplicated(hash), n_genomes >= 10, tot_host >= 4, tot_hostn2 >=2, tot_hostgroups >= 2, n_genomes <= 500, n_fam==1) %>% arrange(-n_genomes)
groups_all_candidates = groups_all %>% filter(group %in% groups_all_stats_candidates$group)

dir.create("../cospec_analysis", recursive=T)

saveRDS(groups_all,"../cospec_analysis/groups_all.Rds")
saveRDS(groups_all_candidates,"../cospec_analysis/groups_all_candidates.Rds")
saveRDS(groups_all_stats, "../cospec_analysis/groups_all_stats.Rds")
saveRDS(groups_all_stats_candidates, "../cospec_analysis/groups_all_stats_candidates.Rds")
write.table(groups_all_stats_candidates, "../cospec_analysis/groups_all_stats_candidates.tsv", sep="\t", row.names=F, quote=F)
write.table(groups_all_candidates, "../cospec_analysis/groups_all_candidates.tsv", sep="\t", row.names=F, quote=F)
