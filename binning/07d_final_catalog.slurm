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

salmon index -t ${PROJECTID}.catalogue.fna -i ${PROJECTID}_salmon_index -p 20

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
/work_beegfs/sukmb276/software/bin/parallel -j 24 --colsep '\t' prodigal -i {1} -a cleanbins_faa/{2}.faa :::: $workfolder/allgroups/${PROJECTID}.cleanbins_list
#/work_beegfs/sukmb276/software/bin/parallel -j 24 --colsep '\t' prodigal -i {1} -a cleanbins_faa/{2}.faa :::: /work_beegfs/sukmb276/Metagenomes/reference_data/UHGG_v2/UHGG.cleanbins_list
cp -r cleanbins_faa $workfolder/allgroups

cat cleanbins_faa/*.faa | gzip -c > protein_catalogs/${PROJECTID}.allprot.faa.gz
cat cleanbins_faa/*.faa > $TMPDIR/${PROJECTID}.allprot.faa
gzip $TMPDIR/${PROJECTID}.allprot.faa
cd cleanbins_faa 
grep '>' * | awk '{print $1}' | sed -r 's/[.]faa[:][>]/\t/' > $workfolder/allgroups/${PROJECTID}.prot_to_genome.tsv

cd $TMPDIR
mmseqs easy-linclust --min-seq-id 0.95 --threads 24 --cov-mode 1 -c 0.8 --kmer-per-seq 80 ${PROJECTID}.allprot.faa ${PROJECTID}.prot95 $TMPDIR
#mmseqs easy-linclust --min-seq-id 0.80 --threads 24 --cov-mode 1 -c 0.8 --kmer-per-seq 80 allprot.faa GreatApes.prot80 $TMPDIR
#mmseqs easy-linclust --min-seq-id 0.50 --threads 24 --cov-mode 1 -c 0.8 --kmer-per-seq 80 allprot.faa GreatApes.prot50 $TMPDIR
mkdir $workfolder/allgroups/protein_catalogs/
cp ${PROJECTID}.prot*_cluster.tsv $workfolder/allgroups/protein_catalogs/
cp ${PROJECTID}.prot*_rep_seq.fasta $workfolder/allgroups/protein_catalogs/


mkdir $TMPDIR/emapper_tmp
cp -r /work_beegfs/ikmb_repository/databases/eggnog_mapper/220425/ eggnog_db
emapper.py --cpu  24 -i ${PROJECTID}.prot95_rep_seq.fasta --tax_scope 1,2,2157,10239 --itype proteins --data_dir /work_beegfs/ikmb_repository/databases/eggnog_mapper/220425/ --output ${PROJECTID}.prot95_rep_seq --temp_dir $TMPDIR/emapper_tmp --dbmem -m mmseqs

# mkdir $workfolder/allgroups/protein_catalogs/emapper_batches
# cut -f 1 ${PROJECTID}.prot95_cluster.tsv | uniq > ${PROJECTID}.prot95_cluster_repseq.list
# split ${PROJECTID}.prot95_cluster_repseq.list -l 100000 -a 3 -d $workfolder/allgroups/protein_catalogs/emapper_batches/batch_

####
module load miniconda3
conda activate mgx_tools_env
btchs=$(ls $workfolder/allgroups/protein_catalogs/emapper_batches/ | head -n 1 | tail -n 1)
cp $workfolder/allgroups/protein_catalogs/emapper_batches/$btchs $btchs.emapper.hits
mkdir $TMPDIR/${btchs}_tmp
emapper.py --cpu  ${SLURM_CPUS_PER_TASK} -i ${PROJECTID}.prot95_rep_seq.fasta --tax_scope 1,2,2157,10239 --itype proteins --data_dir /work_beegfs/ikmb_repository/databases/eggnog_mapper/220425/ --output $btchs --temp_dir $TMPDIR/${btchs}_tmp --dbmem -m mmseqs --resume


conda activate r_microbiome_env

R --vanilla < """
library(tidyverse)
prot_to_cluster = read_delim('GreatApes.prot95_cluster.tsv', delim='\t', col_names=c('cluster_id','prot_id'))

cluster_annot = read_delim('GreatApes.prot95_rep_seq.emapper.annotations',col_names=T, delim='\t', quote='', num_threads=10,comment='', skip=4)

cluster_annot_parts = lapply(c('test/GreatApes.prot95_rep_seq_part_05.emapper.annotations', 'test/GreatApes.prot95_rep_seq_part_06.emapper.annotations', 
    'test/GreatApes.prot95_rep_seq_part_07.emapper.annotations', 'test/GreatApes.prot95_rep_seq_part_08.emapper.annotations', 'test/GreatApes.prot95_rep_seq_part_09.emapper.annotations'), function(x)
 read_delim(x, delim='\t', quote='', num_threads=10,comment='', skip=4, col_names=colnames(cluster_annot))) %>% do.call('bind_rows', .)

cluster_annot = bind_rows(cluster_annot, cluster_annot_parts)
colnames(cluster_annot)[1] = 'cluster_id'

cluster_annot = cluster_annot %>% filter(!duplicated(cluster_id))

prot_to_genome = read_delim('../GreatApes.prot_to_genome.tsv', delim='\t', quote='', num_threads=10,col_names=c('mag_id','prot_id'))
genome_to_sgb = read_delim('../GreatApes.all_genomes.tsv', delim='\t',col_names=T) %>% select(mag_id=genome, GreatApes_95_final) %>% filter(!is.na(GreatApes_95_final)) %>% mutate(sample=substr(mag_id,1,9))
host = read_delim('../../../../groupings_nozoo.txt', delim='\t', col_names=c('sample', 'host')) %>% mutate(host_group=substr(host,1,1)) %>% bind_rows(., data.frame(sample='MGYG00000', host='HsUHGG',host_group='H')) %>% select(-host)

sgb_host_to_cluster = genome_to_sgb %>% left_join(host) %>% left_join(prot_to_genome) %>% left_join(prot_to_cluster) %>% select(GreatApes_95_final, host_group, cluster_id) %>% distinct()
sgb_host_w_annot = sgb_host_to_cluster  %>% left_join(cluster_annot)

for(annot in c('Preferred_name', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs')){
    print(annot)
    this_annot = sgb_host_w_annot %>% select(GreatApes_95_final, host_group, ogs=!!annot) %>% filter(ogs!='-') %>% separate_rows(ogs) %>% distinct 
    saveRDS(this_annot, paste0('GreatApes.',annot,'_to_sgb.Rds'))
}

saveRDS(cluster_annot, 'GreatApes.emapper.Rds')

"""




#
#
# source activate mashtee_env
# cut -f 1 GMbC.representatives_list > GMbC.representatives_genomes
# mashtree --numcpus 24 --file-of-files GMbC.representatives_genomes --outtree GMbC.mashtree.tre
