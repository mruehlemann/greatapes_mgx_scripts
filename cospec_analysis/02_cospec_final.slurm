#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=20gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/cospec_analysis/log/%A_%a.out
#SBATCH --job-name="ca_02"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

###########################
####    SUbGROUP selection    ########
####################################
cd $workfolder/cospec_analysis

#zcat ../allgroups/GreatApes*.user_msa.fasta.gz > tmp_all_msa.fasta

allsubgroups=($(cut -f 1 groups_all_stats_candidates.tsv | tr -d '\"'))
thissubgroup=${allsubgroups[$SLURM_ARRAY_TASK_ID]}

genomes=$(awk -v tg=$thissubgroup '{if($2 == tg){print $3}}' groups_all_candidates.tsv)
this_fam=$(echo $thissubgroup | cut -d '_' -f 1-3)
this_ord=$(grep -w $this_fam ../allgroups/GreatApes.cluster_final_tax.tsv | head -n 1 | cut -f 2 | cut -d ';' -f 4)
this_dom=$(grep -w $this_fam ../allgroups/GreatApes.cluster_final_tax.tsv | head -n 1 | cut -f 2 | cut -d ';' -f 1)

mkdir -p subgroups/${thissubgroup}
echo ${genomes[@]} | tr ' ' "\n" | /work_beegfs/sukmb276/software/bin/faSomeRecords tmp_all_msa.fasta /dev/stdin subgroups/${thissubgroup}/${thissubgroup}.user_msa.fasta
og=$(grep -w $this_ord ../allgroups/GreatApes.cluster_final_tax.tsv | grep -v $this_fam | cut -f 1 | xargs -I {} grep -w {} ../allgroups/GreatApes.dRep_cluster95_representatives.tsv | awk '{if($4>.8) print}' | shuf -n 1 | cut -f 1)
if [ "$og" == "" ]; then
this_cla=$(grep -w $this_fam ../allgroups/GreatApes.cluster_final_tax.tsv | head -n 1 | cut -f 2 | cut -d ';' -f 3)
og=$(grep -w $this_cla ../allgroups/GreatApes.cluster_final_tax.tsv | grep -v $this_fam | cut -f 1 | xargs -I {} grep -w {} ../allgroups/GreatApes.dRep_cluster95_representatives.tsv | awk '{if($4>.8) print}' | shuf -n 1 | cut -f 1)
fi

if [ "$og" == "" ]; then
this_phy=$(grep -w $this_fam ../allgroups/GreatApes.cluster_final_tax.tsv | head -n 1 | cut -f 2 | cut -d ';' -f 2)
og=$(grep -w $this_phy ../allgroups/GreatApes.cluster_final_tax.tsv | grep -v $this_fam | cut -f 1 | xargs -I {} grep -w {} ../allgroups/GreatApes.dRep_cluster95_representatives.tsv | awk '{if($4>.8) print}' | shuf -n 1 | cut -f 1)
fi

if [ "$og" == "" ]; then
og=$(grep -w $this_dom ../allgroups/GreatApes.cluster_final_tax.tsv | grep -v $this_fam | cut -f 1 | xargs -I {} grep -w {} ../allgroups/GreatApes.dRep_cluster95_representatives.tsv | awk '{if($4>.8) print}' | shuf -n 1 | cut -f 1)
fi


echo $og | /work_beegfs/sukmb276/software/bin/faSomeRecords tmp_all_msa.fasta /dev/stdin subgroups/${thissubgroup}/${thissubgroup}.og.fasta

cd subgroups/${thissubgroup}
cat ${thissubgroup}.og.fasta | grep -v '>' | awk 'BEGIN{print ">outgroup"}{print}' > ${thissubgroup}.faa
cat ${thissubgroup}.user_msa.fasta >> ${thissubgroup}.faa

#module load miniconda3
#source activate gtdbtk_2.1.0

/work_beegfs/sukmb276/software/bin/iqtree2 -s ${thissubgroup}.faa -T $SLURM_CPUS_PER_TASK --alrt 1000 -B 1000 -m WAG


source activate r_microbiome_env
Rscript $workfolder/cospec_analysis/cospec_final.R