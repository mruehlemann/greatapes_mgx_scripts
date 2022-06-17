library(dplyr)

local=rev(strsplit(getwd(),split="/")[[1]])[1]

cluster_95 = read.table(paste0(local,".dRep_fastANI_95.csv"), stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_95=paste0(secondary_cluster)) %>% select(genome, cluster_95)
cluster_97 = read.table(paste0(local,".dRep_fastANI_97.csv"), stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_97=paste0(secondary_cluster)) %>% select(genome, cluster_97)
cluster_99 = read.table(paste0(local,".dRep_fastANI_99.csv"), stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_99=paste0(secondary_cluster)) %>% select(genome, cluster_99)

clusters= cluster_95 %>% left_join(cluster_97)  %>% left_join(cluster_99)


scores = read.table(paste0(local,".refined.out"), head=T, stringsAsFactors=F) %>% tibble::rownames_to_column("genome") %>% mutate(score.gtdb_rel207_ar53 = ifelse(is.na(score.gtdb_rel207_ar53), 0, score.gtdb_rel207_ar53), score.gtdb_rel207_bac120 = ifelse(is.na(score.gtdb_rel207_bac120), 0, score.gtdb_rel207_ar53)) %>%
	mutate(domain=ifelse(max==score.gtdb_rel207_ar53,"Archaea","Bacteria"), genome=gsub("fasta","fa",genome)) %>% select(genome, domain, score=max)

clusters = left_join(clusters, scores)

cl95_reps = clusters %>% arrange(cluster_95, -score) %>% distinct_at(.vars="cluster_95", .keep_all=T) %>% select(-cluster_99, -cluster_97) %>% arrange(-score) %>% mutate(cluster_95_final=paste0(local,"_cluster95_",formatC(seq_along(cluster_95), width = 6, format = "d", flag = "0")))
cl97_reps = clusters %>% arrange(cluster_97, -score) %>% distinct_at(.vars="cluster_97", .keep_all=T) %>% select(-cluster_99, -cluster_95) %>% arrange(-score) %>% mutate(cluster_97_final=paste0(local,"_cluster97_",formatC(seq_along(cluster_97), width = 6, format = "d", flag = "0")))
cl99_reps = clusters %>% arrange(cluster_95, -score) %>% distinct_at(.vars="cluster_99", .keep_all=T) %>% select(-cluster_95, -cluster_97) %>% arrange(-score) %>% mutate(cluster_99_final=paste0(local,"_cluster99_",formatC(seq_along(cluster_99), width = 6, format = "d", flag = "0")))

clusters = clusters %>% left_join(cl95_reps %>% select(cluster_95, cluster_95_final)) %>% left_join(cl97_reps %>% select(cluster_97, cluster_97_final)) %>% left_join(cl99_reps %>% select(cluster_99, cluster_99_final))

write.table(clusters, paste0(local, ".dRep_cluster.tsv"), sep="\t", row.names=F,quote=F)
write.table(cl95_reps, paste0(local, ".dRep_cluster95_representatives.tsv"), sep="\t", row.names=F,quote=F)
write.table(cl97_reps, paste0(local, ".dRep_cluster97_representatives.tsv"), sep="\t", row.names=F,quote=F)
write.table(cl99_reps, paste0(local, ".dRep_cluster99_representatives.tsv"), sep="\t", row.names=F,quote=F)
