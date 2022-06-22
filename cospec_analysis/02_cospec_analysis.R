library(ape)
library(phangorn)
library(tidyverse)
library(TreeDist)
library(Quartet)

###
this = rev(strsplit(getwd(), split="/")[[1]])[1]

type_recode=data.frame(short=c("Ggorilla","Gberingei","HsA1","HsA2","HsG","Ptschw","Pttrog","Ptver","Ppan"), short2=c("Ggorilla","Gberingei","Hs","Hs","Hs","Ptschw","Pttrog","Ptver","Ppan"), long=c("Gorilla_gorilla_gorilla","Gorilla_beringei_beringei","Homo_sapiens_sapiens","Homo_sapiens_sapiens","Homo_sapiens_sapiens","Pan_troglodytes_schweinfurthii","Pan_troglodytes_troglodytes","Pan_troglodytes_verus","Pan_paniscus"),stringsAsFactors=F)

host = read.table("/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/groupings.txt", stringsAsFactors = F)
colnames(host) = c("sample","host")
all_groups = readRDS("../../groups_all.Rds")
this_groups = all_groups %>% filter(group == this)
all_genomes = read.table("../../../allgroups/GreatApes.all_genomes.tsv", head=T, stringsAsFactors=F) %>%
  filter(!is.na(GreatApes_95_final), genome %in% this_groups$genome) %>% mutate(sample=substr(genome,1,9)) %>% left_join(host) %>% left_join(type_recode, by=c("host"="short"))


### alignemnt to distance matrix
aln <- read.phyDat(paste0(this, ".user_msa.fasta"), type="AA", format="fasta")
aln2=as.character(aln)

para.D  <- as.matrix(dist.ml(aln, "WAG"))
saveRDS(para.D, paste0(this, ".user_msa.matrix.Rds"))

source("~/Isilon/Metagenomes/software/Rscripts/mad.R")

dir.create("trees",recursive=T, showWarnings=F)

# tree_fastme = fastme.bal(para.D, nni=T, spr=T, tbr=T)
# tree_fastme_bs=boot.phylo(tree_fastme, aln2, FUN=function(xx) fastme.bal(dist.ml(as.phyDat(xx, type="AA"), "WAG"), nni=T, spr=T, tbr=T), B=1000)
# tree_fastme$node.label = tree_fastme_bs/1000
# write.tree(tree_fastme, paste0("trees/",this, ".aln.fastme.bs.nwk"))
# fastme_tree_rooted<-mad(tree_fastme, output_mode="full")[[6]][[1]]
# write.tree(fastme_tree_rooted, paste0("trees/",this, ".aln.fastme.bs.nwk"))

tree_bionj = bionj(para.D)
tree_bionj_bs=boot.phylo(tree_bionj, aln2, FUN=function(xx) bionj(dist.ml(as.phyDat(xx, type="AA"), "WAG")), B=1000, mc.cores = 1)
tree_bionj$node.label = tree_bionj_bs/1000
write.tree(tree_bionj,  paste0("trees/",this, ".aln.bionj.bs.nwk"))
bionj_tree_rooted<-mad(tree_bionj, output_mode="full")[[6]][[1]]
write.tree(bionj_tree_rooted,  paste0("trees/",this, ".aln.bionj.bs.nwk"))

this_tree = tree_bionj
this_tree_labels = this_tree$tip.label

### load functions

### Test Prep

host.D <- as.matrix(readDist("/home/sukmb276/Isilon/references/mtDNA/subset/pwdistmat"))
rownames(host.D)<-colnames(host.D)<-sapply(rownames(host.D), function(x){gsub("[.]","_",strsplit(x,split="-")[[1]][2])})

host_tree = mad(bionj(host.D), output_mode="full")[[6]][[1]]
this_host_tree = drop.tip(host_tree, type_recode %>% filter(!short %in% all_genomes$host) %>% pull(long))

n_host = length(this_host_tree$tip.label)
n_genomes = nrow(all_genomes)

nperm_initial = 100

source("/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts/cospec_analysis/02_functions.R")

### treedist based

nperm=nperm_initial
all_treedist_p=cospec_treebased(all_genomes, this_tree, this_host_tree, nperm)
p_final = 1-colSums((all_treedist_p %>% t) < 0.05)/(nperm+1)

if(any(p_final < 0.1)){
all_treedist_p2=cospec_treebased(all_genomes, this_tree, this_host_tree, 900)
all_treedist_p = cbind(all_treedist_p, all_treedist_p2)
nperm=1000
p_final = 1-colSums((all_treedist_p %>% t) < 0.05)/(nperm+1)
}

output_treebased = data.frame(p_final) %>% t %>% data.frame(subtree=this, n_host=n_host, n_genomes=n_genomes, nperm_treebased = nperm, .)

write.table(output_treebased, paste0(this,".cospec.out"), row.names=F, sep="\t")


### Hommola

nperm_hommola=nperm_initial
all_cospec_hommola=cospec_hommola(all_genomes, para.D, host.D, nperm_hommola)
p_final_hommola = 1-sum(all_cospec_hommola < 0.05)/(nperm_hommola+1)

if(p_final_hommola < 0.1){
all_cospec_hommola2=cospec_hommola(all_genomes, para.D, host.D, 900)
all_hommola_p = c(all_cospec_hommola, all_cospec_hommola2)
nperm_hommola=1000
p_final_hommola = 1-colSums((all_hommola_p %>% t) < 0.05)/(nperm_hommola+1)
}

output_hommola = data.frame(output_treebased, nperm_hommola=nperm_hommola, p_hommola=p_final_hommola)

write.table(output_hommola, paste0(this,".cospec.out"), row.names=F, sep="\t")


### ParaFit

nperm_parafit=nperm_initial
all_cospec_parafit=cospec_parafit(all_genomes, para.D, host.D, nperm_hommola)
p_final_parafit = 1-sum(all_cospec_parafit < 0.05)/(nperm_parafit+1)

if(p_final_parafit < 0.1){
all_cospec_parafit2=cospec_parafit(all_genomes, para.D, host.D, 900)
all_parafit_p = c(all_cospec_parafit, all_cospec_parafit2)
nperm_parafit=1000
p_final_parafit = 1-colSums((all_hommola_p %>% t) < 0.05)/(nperm_parafit+1)
}

output_parafit = data.frame(output_hommola, nperm_parafit=nperm_parafit, p_parafit=p_final_parafit)

write.table(output_parafit, paste0(this,".cospec.out"), row.names=F, sep="\t")
