library(dplyr)

clusters = read.table("GreatApes.rep97.dRep_95.csv", stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_95=paste0(secondary_cluster)) %>% select(genome, cluster_95) %>% mutate(genome=gsub(".fa$","",genome))


scores = read.table("GreatApes.rep97.refined.out", head=T, stringsAsFactors=F) %>% tibble::rownames_to_column("genome") %>% mutate(score.gtdb_rel207_ar53 = ifelse(is.na(score.gtdb_rel207_ar53), 0, score.gtdb_rel207_ar53), score.gtdb_rel207_bac120 = ifelse(is.na(score.gtdb_rel207_bac120), 0, score.gtdb_rel207_ar53)) %>%
	mutate(domain=ifelse(max==score.gtdb_rel207_ar53,"Archaea","Bacteria"), genome=gsub("fasta","fa",genome)) %>% select(genome, domain, score=max)

clusters = left_join(clusters, scores)

cl95_reps = clusters %>% arrange(cluster_95, -score) %>% distinct_at(.vars="cluster_95", .keep_all=T) %>% arrange(-score) %>% mutate(cluster_95_final=paste0("GreatApes","_cluster95_",formatC(seq_along(cluster_95), width = 6, format = "d", flag = "0")))

clusters = clusters %>% left_join(cl95_reps %>% select(cluster_95, cluster_95_final))

write.table(clusters, paste0("GreatApes", ".dRep_cluster.tsv"), sep="\t", row.names=F,quote=F)
write.table(cl95_reps, paste0("GreatApes", ".dRep_cluster95_representatives.tsv"), sep="\t", row.names=F,quote=F)

all_sub = list.files("../subgroups/")
all_sub = all_sub[all_sub!="alltax"]
all_clusters = sapply(all_sub, function(x) read.table(paste0("../subgroups/",x,"/",x,".dRep_cluster.tsv"), head=T, stringsAsFactors=F), simplify=F) %>% do.call("rbind",.)

all_reps97 = sapply(all_sub, function(x) read.table(paste0("../subgroups/",x,"/",x,".dRep_cluster97_representatives.tsv"), head=T, stringsAsFactors=F), simplify=F) %>% do.call("rbind",.)

all_clusters_w_rep = all_clusters %>% left_join(all_reps97 %>% select(cluster_97_final, genome_rep97=genome)) %>% left_join(clusters %>% select(genome_rep97=genome, GreatApes_95_final = cluster_95_final))

write.table(all_clusters_w_rep, paste0("GreatApes", ".all_genomes.tsv"), sep="\t", row.names=F,quote=F)
