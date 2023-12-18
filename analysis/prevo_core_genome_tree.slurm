#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80Gb
#SBATCH --time=0-08:00
#SBATCH --output=/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="ranger"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
workfolder="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final"
source $scriptdir/00_sources.txt
#####
module load miniconda3
conda activate /work_beegfs/sukmb276/Cerascreen_Data/cerascreen_new_env/

all_genes=($(cut -f 1 $workfolder/analyses/host_core_genome_minprev_rand.tsv))


mkdir -p $workfolder/analyses/1krandom_analysis


R --vanilla <<< '''
library(Biostrings)
library(tidyverse)

genes=read_tsv("../host_core_genome_minprev_rand.tsv")
genes_prev = genes %>% filter(rel_across >= .95)

allgene =  read_tsv("../analysis_gene_set_complete.tsv")


all_aln = lapply(genes_prev$Preferred_name, function(g){
    thisgene = allgene %>% filter(Preferred_name == g)
    thisaln = readAAMultipleAlignment(paste0(g,"/",g,".aln.faa"))
    rownames(thisaln) = rownames(thisaln) %>% gsub(" .+","",.)
    rownames(thisaln) = thisgene[match(rownames(thisaln), thisgene$prot_id),"GreatApes_95_final"] %>% unlist
    return(thisaln)}
)

dss2df <- function(dss) data.frame(names=rownames(dss), seq=as.character(dss), stringsAsFactors=F)

rep = which(lapply(all_aln, function(x) dim(x)[1]) == 185)[1]

all_aln_start = dss2df(all_aln[[rep]])
all_aln_next = all_aln_start

for(i in (1:length(all_aln))[-rep]){
    all_aln_next = all_aln_next %>% 
        left_join(dss2df(all_aln[[i]]) %>% select(names, seq2=seq)) %>% 
        mutate(seq2=ifelse(is.na(seq2), paste0(rep("-",  dim(all_aln[[i]])[2]), collapse=""), seq2)) %>% 
        mutate(seq = paste0(seq, seq2)) %>% 
        select(names, seq)
}

all_aln_merged = AAMultipleAlignment(x=all_aln_next$seq, use.names=F)
rownames(all_aln_merged) = all_aln_next$names
AAStr=as(all_aln_merged,"AAStringSet")
writeXStringSet(AAStr,file="all_aln_merged.faa")

'''

/work_beegfs/sukmb276/software/bin/iqtree2 -s all_aln_merged.faa -T $SLURM_CPUS_PER_TASK --alrt 1000 -B 1000 -M -mset LG,WAG



R --vanilla <<< '''
library(ape)
library(phangorn)
library(dplyr)

spetree=read.tree("../../markertree_ranger_w_bl.tre")

gsc=read.table("../../analysis_gene_set_complete.tsv", head=T, stringsAsFactors=F)
rena= gsc %>% select(cluster=GreatApes_95_final, genus) %>%
    mutate(rename=gsub("_cluster95_","",cluster)) %>% distinct %>%
    filter(rename %in% spetree$tip.label)
genetree = read.tree(list.files(pattern=".aln.faa.contree"))
#rena=read.table("../../rename_list.tsv", head=T, stringsAsFactors=F)
genes=read.table("geneset",head=F,stringsAsFactors=F)
genes = genes %>% left_join(rena, by=c("V1"="cluster")) %>% 
    mutate(genename=paste0(rename,"_gene1")) %>% 
    select(geneid=V4, genename, genome=V2, genome_rename=rename) %>% 
    filter(!is.na(genome_rename))

genetree = genetree %>% keep.tip(genes$geneid)
genetree$tip.label = genes[match(genetree$tip.label, genes$geneid),"genename"]
genetree$edge.length = NULL
genetree$node.label = (genetree$node.label %>% as.numeric)
genetree$node.label = ifelse(is.na(genetree$node.label),"",genetree$node.label)
write.tree(genetree, "genetree_ranger.tre")

genetree2 = genetree
genetree2$tip.label = gsub("_gene[0-9]","",genetree2$tip.label)

spetree2= keep.tip(spetree, genetree2$tip.label)


'''


cat ../../markertree_ranger_w_bl.tre genetree_ranger.tre > ranger_input.tre
../../../analyses_old2//cydab/Linux/SupplementaryPrograms/OptResolutions.linux -i ranger_input.tre -B 50 -N 100 | tee ranger_optresolutions.log
grep Optimal ranger_optresolutions.log -A 1 | tail -n 1 > $gene.ranger_optim.tre


cat ../../markertree_ranger_w_bl.tre genetree_ranger.tre > ranger_input_optim.tre
mkdir -p recon

/work_beegfs/sukmb276/software/bin/parallel -j $SLURM_CPUS_PER_TASK ../../../analyses_old2//cydab/Linux/CorePrograms/Ranger-DTL.linux --seed {} -i ranger_input_optim.tre -o recon/rangerOutput{} ::: $(seq 1 100) | tee ranger_core.log
../../../analyses_old2//cydab/Linux/CorePrograms/AggregateRanger.linux recon/rangerOutput > ${gene}.ranger_AggregateOutput.txt
grep -E "is: [0-9]+" recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Cost",i, c, c/i}' > ${gene}.ranger_final_stats.tsv
grep -E "Duplications: [0-9]+" recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Duplications",i, c, c/i}' >> ${gene}.ranger_final_stats.tsv
grep -E "Transfers: [0-9]+" recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Transfers",i, c, c/i}' >> ${gene}.ranger_final_stats.tsv
grep -E "Losses: [0-9]+" recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Losses",i, c, c/i}' >> ${gene}.ranger_final_stats.tsv
grep -E -o "Alignment has [0-9]+ sequences" ${gene}.aln.faa.log | head -n 1 |  cut -d " " -f 3 | awk '{print "N_seq", $1, $1, $1}' >> ${gene}.ranger_final_stats.tsv

grep Transfer, recon/* | tr ":" "," | cut -d "," -f 1,5-6 | sed 's/ Mapping --> //' | sed 's/ Recipient --> //' | sed 's/recon\///' | awk '{print "transfer,"$0}' > $gene.events.csv
grep Speciation, recon/* | tr ":" "," | cut -d "," -f 1,5 | sed 's/ Mapping --> //'  | sed 's/recon\///' | awk '{print "speciation,"$0","}' >> $gene.events.csv

grep "Species Tree" recon/rangerOutput1 -A 1 | tail -n 1 > species_tree.nodenames.tre

exit




