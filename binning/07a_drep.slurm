#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_04b_vamb"


###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt


####################################
####    SAMPLE SELECTION     ########
####################################

## selects all samples that belong to the respective subgroup (locality = column 3 of the metadata sheet)
this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples_indices=$(awk -veth=${this_local} '{if($2==eth) print (NR-1)}' $groupingfile | tr '\n' ',' | sed 's/[,]$//')
local_samples=$(awk -veth=${this_local} '{if($2==eth) print $1}' $groupingfile)

###########################
####    drep        #######
###########################

module load miniconda3

cd $workfolder/subgroups/$this_local/
#mkdir -p genome_links

for i in $local_samples; do ls $workfolder/samples/$i/binning/cleanbins/*.fa; done | awk -F '/' '{print $0"\t"$NF}' | sed 's/[.]fa$//' >  $workfolder/subgroups/$this_local/$this_local.cleanbins_list

head -n 1 $workfolder/samples/$(ls $workfolder/samples  | head -n 1 | cut -d ' ' -f 1)/binning/*.refined.out > ${this_local}.refined.out
for i in $local_samples; do cat $workfolder/samples/$i/binning/$i.refined.out | grep -v "uniqueSCGs.gtdb_rel207_bac120"; done >> ${this_local}.refined.out

source activate r_microbiome_env

R --vanilla <<< "
library(tidyverse)
this=rev(strsplit(getwd(),split='/')[[1]])[1]
scores=read.table(paste0(this,'.refined.out'), head=T, stringsAsFactors=F)
if(has_rownames(scores)==F){scores = scores  %>% column_to_rownames('bin')}
if('max_rat' %in% colnames(scores)){scores = scores %>% mutate(max_complete = max_rat)}
scores = scores %>% rownames_to_column('MAG.ID') %>% mutate(MAG.ID=gsub('.fa(sta)*$','', MAG.ID)) %>% 
    mutate(Domain = ifelse(max == score.gtdb_rel207_bac120, 'Bacteria', 'Archaea')) %>%
    mutate(
        Score = max, Completeness = max_complete, 
        Contamination=ifelse(Domain=='Bacteria', (sumSCGs.gtdb_rel207_bac120 - uniqueSCGs.gtdb_rel207_bac120)/120, (sumSCGs.gtdb_rel207_ar53 - uniqueSCGs.gtdb_rel207_ar53)/ 53)
        ) 

scores_reduced = scores %>% select(MAG.ID, Domain, Score, Completeness, Contamination)
write.table(scores_reduced, paste0(this, '.scores_final.out'), col.names=T, row.names=F, quote=F, sep='\t')


cb=read.table(paste0(this,'.cleanbins_list'), head=F, stringsAsFactors=F)

cb_filt = cb %>% filter(V2 %in% (scores_reduced %>% filter(Contamination < .1) %>% pull(MAG.ID)))

write.table(cb_filt, paste0(this, '.cleanbins_filtered'), col.names=F, row.names=F, quote=F, sep='\t')
"

source activate binning_env

cd $TMPDIR
cut -f 1 $workfolder/subgroups/$this_local/$this_local.cleanbins_filtered > genomes_list
dRep compare -g genomes_list -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -sa 0.95 -d .
cp data_tables/Cdb.csv ${this_local}.dRep_fastANI_95.csv
python $scriptdir/binning/07a_drep_fastANI99.py
cp data_tables/Cdb_97.csv ${this_local}.dRep_fastANI_97.csv
cp data_tables/Cdb_99.csv ${this_local}.dRep_fastANI_99.csv
cp ${this_local}.dRep_fastANI_95.csv ${this_local}.dRep_fastANI_97.csv ${this_local}.dRep_fastANI_99.csv $workfolder/subgroups/$this_local
conda deactivate

cd $workfolder/subgroups/$this_local
conda activate r_microbiome_env

R --vanilla <<< "
library(dplyr)

this=rev(strsplit(getwd(),split='/')[[1]])[1]

cluster_95 = read.table(paste0(this,'.dRep_fastANI_95.csv'), stringsAsFactors=F, head=T, sep=',') %>% mutate(cluster_95=paste0(secondary_cluster)) %>% select(genome, cluster_95)
cluster_97 = read.table(paste0(this,'.dRep_fastANI_97.csv'), stringsAsFactors=F, head=T, sep=',') %>% mutate(cluster_97=paste0(secondary_cluster)) %>% select(genome, cluster_97)
cluster_99 = read.table(paste0(this,'.dRep_fastANI_99.csv'), stringsAsFactors=F, head=T, sep=',') %>% mutate(cluster_99=paste0(secondary_cluster)) %>% select(genome, cluster_99)

clusters= cluster_95 %>% left_join(cluster_97)  %>% left_join(cluster_99) %>% mutate(genome=gsub('.fa(sta)*$|.fna$','',genome))

scores_reduced = read.table(paste0(this, '.scores_final.out'), stringsAsFactors=F, head=T)

clusters = left_join(clusters, scores_reduced, by=c('genome'='MAG.ID'))

cl95_reps = clusters %>% arrange(cluster_95, -Score) %>% distinct_at(.vars='cluster_95', .keep_all=T) %>% select(-cluster_99, -cluster_97) %>% arrange(-Score) %>% mutate(cluster_95_final=paste0(this,'_cluster95_',formatC(seq_along(cluster_95), width = 6, format = 'd', flag = '0')))
cl97_reps = clusters %>% arrange(cluster_97, -Score) %>% distinct_at(.vars='cluster_97', .keep_all=T) %>% select(-cluster_99, -cluster_95) %>% arrange(-Score) %>% mutate(cluster_97_final=paste0(this,'_cluster97_',formatC(seq_along(cluster_97), width = 6, format = 'd', flag = '0')))
cl99_reps = clusters %>% arrange(cluster_95, -Score) %>% distinct_at(.vars='cluster_99', .keep_all=T) %>% select(-cluster_95, -cluster_97) %>% arrange(-Score) %>% mutate(cluster_99_final=paste0(this,'_cluster99_',formatC(seq_along(cluster_99), width = 6, format = 'd', flag = '0')))

clusters = clusters %>% left_join(cl95_reps %>% select(cluster_95, cluster_95_final)) %>% 
    left_join(cl97_reps %>% select(cluster_97, cluster_97_final)) %>% 
    left_join(cl99_reps %>% select(cluster_99, cluster_99_final)) %>% select(-cluster_95, -cluster_97, -cluster_99)

write.table(clusters, paste0(this, '.dRep_cluster.tsv'), sep='\t', row.names=F,quote=F)
write.table(cl95_reps, paste0(this, '.dRep_cluster95_representatives.tsv'), sep='\t', row.names=F,quote=F)
write.table(cl97_reps, paste0(this, '.dRep_cluster97_representatives.tsv'), sep='\t', row.names=F,quote=F)
write.table(cl99_reps, paste0(this, '.dRep_cluster99_representatives.tsv'), sep='\t', row.names=F,quote=F)
"

### store primary clustring information

R --vanilla <<< "
library(dplyr)

this=rev(strsplit(getwd(),split='/')[[1]])[1]

cluster_primary = read.table(paste0(this,'.dRep_fastANI_95.csv'), stringsAsFactors=F, head=T, sep=',') %>% 
    mutate(genome=gsub('.fa(sta)*$|.fna$','',genome), 
    primary_final=paste0(this,'_primary_',formatC(primary_cluster, width = 6, format = 'd', flag = '0'))) %>% 
    select(genome, primary_final)

scores_reduced = read.table(paste0(this, '.scores_final.out'), stringsAsFactors=F, head=T)
cluster_primary = left_join(cluster_primary, scores_reduced, by=c('genome'='MAG.ID'))

cluster_primary_reps = cluster_primary %>% arrange(primary_final, -Score) %>% distinct_at(.vars='primary_final', .keep_all=T) %>% arrange(-Score, -Completeness) 

write.table(cluster_primary, paste0(this, '.dRep_primary.tsv'), sep='\t', row.names=F,quote=F)
write.table(cluster_primary_reps, paste0(this, '.dRep_primary_reps.tsv'), sep='\t', row.names=F,quote=F)
"

mkdir -p primary_clustering/representatives

while read line; do
read genome primary domain score comp cont <<< $(echo $line)
if [ "$genome" == "genome" ]; then continue; fi
st=`echo "$score < 0.7" | bc`
if [ $st -eq 1 ]; then echo "$genome is medium quality. skipping" ; continue; fi
path=$(grep -w $genome ${this_local}.cleanbins_filtered | cut -f 1)
echo $path
ln -s $path primary_clustering/representatives/${primary}.fna
done < ${this_local}.dRep_primary_reps.tsv

source activate binning_env

cd primary_clustering/representatives
mash sketch -s 1000 -p 24 -o ../$this_local.msh *