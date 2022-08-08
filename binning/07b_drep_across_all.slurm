#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_04b_vamb"


###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME
scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

####################################
####    SAMPLE SELECTION     ########
####################################

mkdir -p $workfolder/allgroups

###########################
####    DREP all        #######
###########################

module load miniconda3
source activate binning_env

head -n 1 $workfolder/samples/$(ls $workfolder/samples | head -n 1)/*.refined.out > $workfolder/allgroups/${PROJECTID}.rep97.refined.out
head -n 1 $workfolder/samples/$(ls $workfolder/samples | head -n 1)/*.refined.out > $workfolder/allgroups/${PROJECTID}.rep97_mq.refined.out

for this_local in Gberingei Ggorilla HsA1 HsA2 HsG Ppan Ptschw Pttrog Ptver;do
cd $workfolder/subgroups/$this_local/
tail -n+2 ${this_local}.dRep_cluster97_representatives.tsv | awk '{if($4>=.7) print $1}'  | xargs -I {} grep {} ${this_local}.cleanbins_list >>  $workfolder/allgroups/${PROJECTID}.rep97.cleanbins_list
tail -n+2 ${this_local}.dRep_cluster97_representatives.tsv | awk '{if($4>=.7) print $1}'  | xargs -I {} grep {} ${this_local}.refined.out  >>  $workfolder/allgroups/${PROJECTID}.rep97.refined.out
tail -n+2 ${this_local}.dRep_cluster97_representatives.tsv | awk '{if($4>=.5 && $4<.7) print $1}'  | xargs -I {} grep {} ${this_local}.cleanbins_list >>  $workfolder/allgroups/${PROJECTID}.rep97_mq.cleanbins_list
tail -n+2 ${this_local}.dRep_cluster97_representatives.tsv | awk '{if($4>=.5 && $4<.7) print $1}'  | xargs -I {} grep {} ${this_local}.refined.out  >>  $workfolder/allgroups/${PROJECTID}.rep97_mq.refined.out
done

cd $TMPDIR
cut -f 1 $workfolder/allgroups/${PROJECTID}.rep97.cleanbins_list | sed '/^$/d' > genomes_list
dRep compare -g genomes_list -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -sa 0.95 -d .
#cp data_tables/Cdb.csv $workfolder/allgroups/${PROJECTID}.rep97.dRep_95.csv
$PYTHON $scriptdir/binning/07a_drep_fastANI99.py
#cp data_tables/Cdb_97.csv $workfolder/allgroups/${PROJECTID}.rep97.dRep_97.csv
#cp data_tables/Cdb_99.csv $workfolder/allgroups/${PROJECTID}.rep97.dRep_99.csv

R EXC CMD "
library(tidyverse)
clusters = read.table("data_tables/Cdb.csv", stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_95=paste0(secondary_cluster)) %>% select(genome, cluster_95) %>% mutate(genome=gsub(".fa$","",genome))

scores = read.table("/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/GreatApes.rep97.refined.out", head=T, stringsAsFactors=F) %>% tibble::rownames_to_column("genome")  %>% mutate(score.gtdb_rel207_ar53 = ifelse(is.na(score.gtdb_rel207_ar53), 0, score.gtdb_rel207_ar53), score.gtdb_rel207_bac120 = ifelse(is.na(score.gtdb_rel207_bac120), 0, score.gtdb_rel207_ar53)) %>%
	mutate(domain=ifelse(max==score.gtdb_rel207_ar53,"Archaea","Bacteria"), genome=gsub(".fasta","",genome)) %>% select(genome, domain, score=max)

genome_paths=read.table("genomes_list", head=F, stringsAsFactors=F)
genome_paths$genome = sapply(genome_paths$V1, function(x) gsub(".fa","",rev(strsplit(x,split="/")[[1]])[1]))

clusters = left_join(clusters, scores) %>% arrange(-score) %>% distinct_at(.vars=c("cluster_95"), .keep_all=T) %>% left_join(genome_paths)

write.table(clusters$V1, "fastani_ref_list", col.names=F, row.names=F, quote=F)
"

cut -f 1 $workfolder/allgroups/${PROJECTID}.rep97_mq.cleanbins_list > fastani_query_list
fastANI --refList fastani_ref_list --ql fastani_query_list -t $SLURM_CPUS_PER_TASK --fragLen 10000 --minFraction 0.8 -o fastani_mq.out

module load R
R EXC CMD "library(tidyverse); fa=read.table("fastani_mq.out",head=F, stringsAsFactors=F);
fa$V1=sapply(fa$V1, function(x) rev(strsplit(x,split="/")[[1]])[1]);
fa$V2=sapply(fa$V2, function(x) rev(strsplit(x,split="/")[[1]])[1]);
fa_filt = fa %>% filter(fa$V3 >= 95) ;
mq=read.table("fastani_query_list", head=F, stringsAsFactors=F);
mq$genome = sapply(mq$V1, function(x) rev(strsplit(x,split="/")[[1]])[1]) ;
mq_nohit = mq %>% filter(!genome %in% fa_filt$V1);
write.table(mq_nohit$V1, "drep_l2_input", col.names=F, row.names=F, quote=F)
"

dRep compare -g drep_l2_input -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -sa 0.95 -d drep_l2_out

module load R
R EXC CMD "
library(tidyverse);
clusters = read.table("data_tables/Cdb.csv", stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_95=paste0(secondary_cluster)) %>% select(genome, cluster_95) %>% mutate(genome=gsub(".fa$","",genome))

scores = read.table("/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/GreatApes.rep97.refined.out", head=T, stringsAsFactors=F) %>% tibble::rownames_to_column("genome") %>% mutate(score.gtdb_rel207_ar53 = ifelse(is.na(score.gtdb_rel207_ar53), 0, score.gtdb_rel207_ar53), score.gtdb_rel207_bac120 = ifelse(is.na(score.gtdb_rel207_bac120), 0, score.gtdb_rel207_ar53)) %>%
	mutate(domain=ifelse(max==score.gtdb_rel207_ar53,"Archaea","Bacteria"), genome=gsub("fasta","fa",genome)) %>% select(genome, domain, score=max)

clusters = left_join(clusters, scores)

scores2 = read.table("/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/GreatApes.rep97_mq.refined.out", head=T, stringsAsFactors=F) %>% tibble::rownames_to_column("genome") %>% mutate(score.gtdb_rel207_ar53 = ifelse(is.na(score.gtdb_rel207_ar53), 0, score.gtdb_rel207_ar53), score.gtdb_rel207_bac120 = ifelse(is.na(score.gtdb_rel207_bac120), 0, score.gtdb_rel207_ar53)) %>%
	mutate(domain=ifelse(max==score.gtdb_rel207_ar53,"Archaea","Bacteria"), genome=gsub("fasta","fa",genome)) %>% select(genome, domain, score=max)

fa=read.table("fastani_mq.out",head=F, stringsAsFactors=F);
fa$V1=sapply(fa$V1, function(x) gsub(".fa","",rev(strsplit(x,split="/")[[1]])[1]));
fa$V2=sapply(fa$V2, function(x) gsub(".fa","",rev(strsplit(x,split="/")[[1]])[1]));
fa_filt = fa %>% filter(fa$V3 >= 95) %>% left_join(scores2, by=c("V1"="genome")) %>% left_join(clusters %>% select(-score, -domain), by=c("V2"="genome")) %>% arrange(-V3) %>% distinct_at(.vars="V1", .keep_all=T) %>% select(genome=V1, cluster_95, domain, score)

clusters = bind_rows(clusters, fa_filt) %>% arrange(-score)

clusters2 = read.table("drep_l2_out/data_tables/Cdb.csv", stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_95=paste0(secondary_cluster)) %>% select(genome, cluster_95) %>% mutate(genome=gsub(".fa$","",genome))

clusters2 = left_join(clusters2, scores2) %>% left_join(clusters2 %>% group_by(cluster_95) %>% summarize(n=n())) %>% filter(n>1) %>% mutate(cluster_95=paste0(cluster_95,"_MQ")) %>% arrange(-score)

clusters = bind_rows(clusters, clusters2 %>% select(-n))

cl95_reps = clusters %>% arrange(-score, cluster_95) %>% distinct_at(.vars="cluster_95", .keep_all=T) %>% arrange(-score) %>% mutate(cluster_95_final=paste0("GreatApes","_cluster95_",formatC(seq_along(cluster_95), width = 6, format = "d", flag = "0")))

clusters = clusters %>% left_join(cl95_reps %>% select(cluster_95, cluster_95_final))

write.table(clusters, paste0("GreatApes", ".dRep_cluster.tsv"), sep="\t", row.names=F,quote=F)
write.table(cl95_reps, paste0("GreatApes", ".dRep_cluster95_representatives.tsv"), sep="\t", row.names=F,quote=F)

all_sub = list.files("/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final//subgroups/")
all_sub = all_sub[all_sub!="alltax"]
all_clusters = sapply(all_sub, function(x) read.table(paste0("/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final//subgroups/",x,"/",x,".dRep_cluster.tsv"), head=T, stringsAsFactors=F), simplify=F) %>% do.call("rbind",.)

all_reps97 = sapply(all_sub, function(x) read.table(paste0("/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final///subgroups/",x,"/",x,".dRep_cluster97_representatives.tsv"), head=T, stringsAsFactors=F), simplify=F) %>% do.call("rbind",.)

all_clusters_w_rep = all_clusters %>% left_join(all_reps97 %>% select(cluster_97_final, genome_rep97=genome)) %>% left_join(clusters %>% select(genome_rep97=genome, GreatApes_95_final = cluster_95_final))

write.table(all_clusters_w_rep, paste0("GreatApes", ".all_genomes.tsv"), sep="\t", row.names=F,quote=F)
"

conda deactivate

cp ${PROJECTID}.* $workfolder/allgroups

cd $workfolder/allgroups

while read line; do
read genome cluster_95 domain score cluster_95_final <<< $(echo $line)
grep -w $genome ${PROJECTID}.rep97.cleanbins_list | awk -v cl=$cluster_95_final '{print $1"\t"cl}'
done < ${PROJECTID}.dRep_cluster95_representatives.tsv | tee ${PROJECTID}.representatives_list
