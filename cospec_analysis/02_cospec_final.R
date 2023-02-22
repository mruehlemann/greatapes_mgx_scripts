
library(ape)
library(phangorn)
library(tidyverse)
library(TreeDist)
library(Quartet)
library(combinat)
###
this = rev(strsplit(getwd(), split="/")[[1]])[1]

type_recode=data.frame(short=c("Ggorilla","Gberingei","HsA1","HsA2","HsG","Ptschw","Pttrog","Ptver","Ppan","HsUHGG"),
        short2=c("Ggorilla","Gberingei","Hs","Hs","Hs","Ptschw","Pttrog","Ptver","Ppan","Hs"),

long=c("Gorilla_gorilla_gorilla","Gorilla_beringei","Homo_sapiens","Homo_sapiens","Homo_sapiens","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pan_paniscus","Homo_sapiens"),stringsAsFactors=F)

host = read.table("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/groupings_nozoo_campbell.txt", stringsAsFactors = F)
colnames(host) = c("sample","host")
host = bind_rows(host, data.frame(sample="MGYG00000",host="HsUHGG")) %>% mutate(host = case_when(host=="ZcambellGorilla" ~ "Ggorilla", host=="ZcambellPan"~ "Pttrog", TRUE ~ host),
    host_genus = ifelse(grepl("^P", host),"Pan",ifelse(grepl("^G", host),"Gorilla","Homo")))

all_groups = readRDS("../../groups_all.Rds")
this_groups = all_groups %>% filter(group == this)

para.D=readDist(paste0(this,".faa.mldist")) %>% as.matrix()

all_genomes = read.table("../../../allgroups/GreatApes.all_genomes.tsv", head=T, stringsAsFactors=F) %>%
  filter(!is.na(GreatApes_95_final), !grepl("Manara", genome), GreatApes_95_final %in% this_groups$GreatApes_95_final, genome %in% rownames(para.D)) %>% mutate(sample=ifelse(grepl("^SRR", genome), substr(genome,1,10), substr(genome,1,9))) %>% left_join(host) %>%
  left_join(type_recode, by=c("host"="short")) %>% arrange(-score) %>% distinct_at(.vars=c("long","GreatApes_95_final"), .keep_all=T)

source("/work_beegfs/sukmb276/Metagenomes/software/Rscripts/mad.R")
bac_tree = read.tree(paste0(this,".faa.treefile"))
bac_tree_rooted = bac_tree %>% root(outgroup="outgroup", resolve.root=T)


dir.create("trees",recursive=T, showWarnings=F)
write.tree(bac_tree,  paste0("trees/",this, ".iqtree2.ml.nwk"))
write.tree(bac_tree_rooted,  paste0("trees/",this, ".iqtree2.ml.rooted.nwk"))


all_genomes = read.table("../../../allgroups/GreatApes.all_genomes.tsv", head=T, stringsAsFactors=F) %>%
  filter(!is.na(GreatApes_95_final), genome %in% bac_tree$tip.label) %>%
  mutate(sample=ifelse(grepl("^SRR", genome), substr(genome,1,10), substr(genome,1,9))) %>%
  left_join(host) %>% left_join(type_recode, by=c("host"="short")) %>% arrange(-score) %>%
  distinct_at(.vars=c("long","GreatApes_95_final"), .keep_all=T)


this_tree = bac_tree
this_tree_labels = this_tree$tip.label

##############

host_tree = mad(read.nexus("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/cospec_analysis/host_trees_10kTrees/consensusTree_10kTrees_Primates_Version3_chrono.nex"),output_mode="full")[[6]][[1]]
host.D = cophenetic(host_tree)

n_host = length(unique(all_genomes$long))
n_genomes = nrow(all_genomes)

nperm_initial = 999

####

all_genome_comb = do.call(expand.grid, split(all_genomes$genome, all_genomes$long))

all_genome_comb_sub = all_genome_comb %>% sample_n(nperm_initial, replace=T)

####

p_rf_rnd=apply(all_genome_comb_sub,1,function(co){
  genomes_sub = all_genomes %>% filter(genome %in% as.character(unlist(co)))

  bac_tree_sub = keep.tip(this_tree, genomes_sub$genome)
  host_tree_sub = keep.tip(host_tree, genomes_sub$long)

  ####

  assoc= genomes_sub %>% select(long, genome)
  HP = assoc %>% mutate(count=1) %>% reshape2::dcast(long ~ genome, value.var="count", fill=0, data=.) %>% column_to_rownames("long")

  bac_tree_sub_renam = bac_tree_sub
  bac_tree_sub_renam$tip.label = assoc[match(bac_tree_sub_renam$tip.label,assoc$genome),"long"]

  ### RF-distance based
  cosp=cospeciation(host_tree_sub, bac_tree_sub, nsim=2, assoc=as.matrix(assoc))$d

  allperm = 1:999
  cosp_nperm = length(allperm)
  cosp_permut=lapply(allperm, function(z){
    assoc_shuf = assoc
    assoc_shuf$long = sample(assoc_shuf$long)
    cospeciation(host_tree_sub, bac_tree_sub, nsim=2, assoc=as.matrix(assoc_shuf))$d %>% return()})

  cosp_p = mean(c(cosp, unlist(cosp_permut)) <= cosp)
  return(cosp_p)

})

#####

p_homm_rnd = apply(all_genome_comb_sub,1,function(co){
  genomes_sub = all_genomes %>% filter(genome %in% as.character(unlist(co)))
  bac_tree_sub = keep.tip(this_tree, genomes_sub$genome)
  para.D.sub=para.D[ genomes_sub$genome, genomes_sub$genome]

  co1 = para.D.sub %>% reshape2::melt() %>%
    left_join(genomes_sub %>% select(genome1=genome, long1=long), by=c("Var1"="genome1")) %>%
    left_join(genomes_sub %>% select(genome2=genome, long2=long), by=c("Var2"="genome2")) %>%
    filter(Var1!=Var2, Var1>Var2) %>%
    left_join(reshape2::melt(host.D) %>% select(long1=Var1, long2=Var2, hdist=value), by=c("long1","long2")) %>%
    (function(df) cor(x=df$hdist, y=df$value))

  allperm = 1:999

  hommola_nperm=length(allperm)

  permut=lapply(allperm, function(z){
    genomes_shuf=genomes_sub
    genomes_shuf$long = sample(genomes_sub$long)
    para.D.sub %>% reshape2::melt() %>%
      left_join(genomes_shuf %>% select(genome1=genome, long1=long), by=c("Var1"="genome1")) %>%
      left_join(genomes_shuf %>% select(genome2=genome, long2=long), by=c("Var2"="genome2")) %>%
      filter(Var1!=Var2, Var1>Var2) %>%
      left_join(reshape2::melt(host.D) %>% select(long1=Var1, long2=Var2, hdist=value), by=c("long1","long2")) %>%
      (function(df) cor(x=df$hdist, y=df$value))})

  hommola_p = mean(c(co1, unlist(permut)) >= co1)
  return(hommola_p)
  })


### hommola all tips

genomes_sub = all_genomes

bac_tree_sub = keep.tip(this_tree, genomes_sub$genome)

para.D.sub=para.D

co1 = para.D.sub %>% reshape2::melt() %>%
  left_join(genomes_sub %>% select(genome1=genome, long1=long), by=c("Var1"="genome1")) %>%
  left_join(genomes_sub %>% select(genome2=genome, long2=long), by=c("Var2"="genome2")) %>%
  filter(Var1!=Var2, Var1>Var2) %>%
  left_join(reshape2::melt(host.D) %>% select(long1=Var1, long2=Var2, hdist=value), by=c("long1","long2")) %>%
  (function(df) cor(x=df$hdist, y=df$value))

allperm=1:999
hommola_nperm=length(allperm)

permut=lapply(allperm, function(z){
  genomes_shuf=genomes_sub
  genomes_shuf$long = sample(genomes_sub$long)
  para.D.sub %>% reshape2::melt() %>%
    left_join(genomes_shuf %>% select(genome1=genome, long1=long), by=c("Var1"="genome1")) %>%
    left_join(genomes_shuf %>% select(genome2=genome, long2=long), by=c("Var2"="genome2")) %>%
    filter(Var1!=Var2, Var1>Var2) %>%
    left_join(reshape2::melt(host.D) %>% select(long1=Var1, long2=Var2, hdist=value), by=c("long1","long2")) %>%
    (function(df) cor(x=df$hdist, y=df$value))})

hommola_p_alltip = mean(c(co1, unlist(permut)) >= co1)




output = list(clade=this, n_host = n_host, n_genomes = n_genomes, tree_permut=length(p_homm_rnd),
  rf_p_mean=mean(p_rf_rnd), rf_p_b01=1-mean(p_rf_rnd<.1), rf_p_b005=1-mean(p_rf_rnd<.05), hommola_p_mean=mean(p_homm_rnd), hommola_p_b01=1-mean(p_homm_rnd<.1), hommola_p_b005=1-mean(p_homm_rnd<.05), hommola_p_alltip=hommola_p_alltip)

write.table(output, paste0(this,".out"), row.names=F, sep="\t")

saveRDS(list(clade=this, p_homm_rnd=p_homm_rnd, p_rf_rnd=p_rf_rnd), paste0(this,".permut_p.Rds"))