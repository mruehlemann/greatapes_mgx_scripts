#!/bin/bash
#SBATCH --reservation=sukmb412_13
#SBATCH --gres=gpu:1
#SBACTH --nodelist=medg01
#SBATCH -c 16
#SBATCH --mem=500gb
#SBATCH --time=3-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_04b_vamb"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME
scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

###########################
####    VAMB        #######
###########################
module load miniconda3

source activate mgx_tools_env
cd $TMPDIR

while read line; do
echo $line
read path genome <<< $(echo $line)
awk -v genome=$genome '/^>/ {print ">"genome"_contig_" ++i; next}{print}' $path > ${genome}.fna
done < $workfolder/allgroups/${PROJECTID}.representatives_list

cat ${PROJECTID}_cluster95_*.fna >  ${PROJECTID}.catalogue.fna

#minimap2 -I100G -d ${PROJECTID}.catalogue.mmi ${PROJECTID}.catalogue.fna # make index

salmon index -t ${PROJECTID}.catalogue.fna -i ${PROJECTID}_salmon_index -p 24

cp -r ${PROJECTID}_salmon_index $workfolder/allgroups/salmon_index

#gzip ${PROJECTID}.catalogue.fna
cp ${PROJECTID}.catalogue.fna.gz  $workfolder/allgroups

gzip ${PROJECTID}_cluster95_00*.fna
mkdir $workfolder/allgroups/reference_collection
mv ${PROJECTID}_cluster95_00*.fna.gz $workfolder/allgroups/reference_collection
#mash dist -d 0.2 -p 24 /work_beegfs/sukmb276/Metagenomes/reference_data/manara_sgb.msh $(cut -f 1 GreatApes.representatives_list | sed 's/\/home\/sukmb276\/Isilon/\/work_beegfs\/sukmb276/') | tee GreatApes.representatives.manara_mash.tsv


conda activate ccmetagen_env
cp $workfolder/allgroups/${PROJECTID}.catalogue.fna.gz .

mkdir -p cleanbins_faa
cat $workfolder/subgroups/*/*.cleanbins_list | sed 's/work_ifs/work_beegfs/' > $workfolder/allgroups/${PROJECTID}.cleanbins_list
cat $workfolder/subgroups/*/*.cleanbins_filtered | sed 's/work_ifs/work_beegfs/' > $workfolder/allgroups/${PROJECTID}.cleanbins_filtered

/work_beegfs/sukmb276/software/bin/parallel -j 24 --colsep '\t' prodigal -i {1} -a cleanbins_faa/{2}.faa :::: $workfolder/allgroups/${PROJECTID}.cleanbins_filtered

cat cleanbins_faa/*.faa > $TMPDIR/${PROJECTID}.allprot.faa

cd $TMPDIR
mmseqs easy-linclust --min-seq-id 0.95 --threads 24 --cov-mode 1 -c 0.8 --kmer-per-seq 80 ${PROJECTID}.allprot.faa ${PROJECTID}.prot95 $TMPDIR

mkdir $workfolder/allgroups/protein_catalogs/

### rename representative genomes 
awk -v collection=$PROJECTID  '/^>/ {s = "0000000000"++i; j=substr(s, 1 + length(s) - 9); print ">"collection"_PROT95_" j; print $1"\t"collection"_PROT95_"j > "cluster_renamed.tsv"; next}{print}' ${PROJECTID}.prot95_rep_seq.fasta > ${PROJECTID}.prot95_rep_seq.faa
sed -i "s/>//" cluster_renamed.tsv
cp cluster_renamed.tsv ${PROJECTID}.PROT95_renamed.tsv


source activate r_microbiome_env
ln -s ${PROJECTID}.prot95_cluster.tsv clustering95.tsv

R --vanilla <<< """
library(tidyverse)
cl = read_tsv('clustering95.tsv', col_names=c('clusterold', 'protid'))
rename=read_tsv('cluster_renamed.tsv', col_names=c('clusterold', 'cluster'))
cl2 = cl %>% left_join(rename) %>% select(cluster, protid)
write.table(cl2, 'cl_renamed.tsv', sep='\t', row.names=F, col.names=F, quote=F)
"""


mv cl_renamed.tsv ${PROJECTID}.prot95_cluster.tsv

cp ${PROJECTID}.PROT95_renamed.tsv ${PROJECTID}.allprot.faa ${PROJECTID}.prot95_rep_seq.faa  ${PROJECTID}.prot95_cluster.tsv $workfolder/allgroups/protein_catalogs/

cp $workfolder/allgroups/protein_catalogs/${PROJECTID}.prot95_rep_seq.faa .
mkdir $TMPDIR/emapper_tmp
emapper.py --cpu ${SLURM_CPUS_PER_TASK} -i ${PROJECTID}.prot95_rep_seq.faa --tax_scope 1,2,2157,10239 --itype proteins --data_dir /work_beegfs/ikmb_repository/databases/eggnog_mapper/220425/ --output ${PROJECTID}.prot95_rep_seq --temp_dir $TMPDIR/emapper_tmp --dbmem -m mmseqs


conda activate r_microbiome_env

R --vanilla < """
library(tidyverse)
cluster_annot = read_delim('protein_catalogs/GreatApes.prot95_rep_seq.emapper.annotations',col_names=c('cluster_id','seed_ortholog','evalue','score','eggNOG_OGs','max_annot_lvl','COG_category','Description','Preferred_name','GOs','EC','KEGG_ko','KEGG_Pathway','KEGG_Module','KEGG_Reaction','KEGG_rclass','BRITE','KEGG_TC','CAZy','BiGG_Reaction','PFAMs'), delim='\t', quote='', num_threads=10,comment='')

host = read_tsv("/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/groupings_nozoo_campbell_dk.txt", col_names=c('sample','host')) %>% 
    mutate(host=ifelse(host=="ZcambellGorilla","Ggorilla", host), host=ifelse(host=="ZcambellPan","Pttrog", host)) %>%
    mutate(host_group=substr(host,1,1)) %>% 
    bind_rows(., data.frame(sample='MGYG00000', host='HsUHGG',host_group='H')) %>% 
    bind_rows(., data.frame(sample='Manara', host='Other',host_group='O')) %>% 
    select(-host) 

all_genomes = read_delim('GreatApes.all_genomes.tsv', delim='\t',col_names=T) %>% 
    filter(!is.na(final)) %>% mutate(GreatApes_95_final=final) %>% 
    filter(!(grepl('MGYG', genome) & genome != cluster_99_final)) %>%
    mutate(sample = case_when(grepl('MGYG', genome) ~ 'MGYG00000', grepl('Manara', genome) ~ 'Manara', T ~ gsub('_cleanbin_[0-9]+','', genome))) %>%
    left_join(host) %>% filter(!duplicated(genome))


prot_to_cluster = read_tsv('protein_catalogs/GreatApes.prot95_cluster.tsv', col_names=c('cluster_id','prot_id')) %>% mutate(contig_id = gsub('_[0-9]+$','',prot_id))
contig_to_sgb = read_tsv('GreatApes.contig_to_bin.tsv', col_names=c('mag_id','contig_id'))
mag_to_sgb = all_genomes %>% select(mag_id=genome, sample, GreatApes_95_final, host_group) 

cluster_to_mag  = prot_to_cluster %>% left_join(contig_to_sgb) %>% select(cluster_id, mag_id)
sgb_host_to_cluster = mag_to_sgb  %>% left_join(cluster_to_mag) %>% select(GreatApes_95_final, host_group, cluster_id) %>% distinct()

sgb_host_w_annot = sgb_host_to_cluster  %>% left_join(cluster_annot)

for(annot in c('KEGG_ko', 'Preferred_name', 'GOs', 'EC', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs')){
    print(annot)
    this_annot = sgb_host_w_annot %>% select(GreatApes_95_final, host_group, ogs=!!annot) %>% filter(ogs!='-') %>% separate_rows(ogs) %>% distinct 
    saveRDS(this_annot, paste0('protein_catalogs/GreatApes.',annot,'_to_sgb.Rds'))
}

saveRDS(cluster_annot, 'GreatApes.emapper.Rds')

"""




#
#
# source activate mashtee_env
# cut -f 1 GMbC.representatives_list > GMbC.representatives_genomes
# mashtree --numcpus 24 --file-of-files GMbC.representatives_genomes --outtree GMbC.mashtree.tre
