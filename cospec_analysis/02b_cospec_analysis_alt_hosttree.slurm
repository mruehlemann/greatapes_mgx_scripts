#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="ca_02"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

###########################
####    SUbGROUP selection    ########
####################################
cd $workfolder/cospec_analysis

#zcat ../allgroups/GreatApes*.user_msa.fasta.gz > tmp_all_msa.fasta

allsubgroups=($(cut -f 1 groups_all_stats_candidates.tsv | tr -d '\"'))
thissubgroup=${allsubgroups[$SLURM_ARRAY_TASK_ID]}

genomes=$(awk -v tg=$thissubgroup '{if($1 == tg){print $3}}' groups_all_candidates.tsv)

mkdir -p subgroups/${thissubgroup}
echo ${genomes[@]} | tr ' ' "\n" | ~/Isilon/software/bin/faSomeRecords tmp_all_msa.fasta /dev/stdin subgroups/${thissubgroup}/${thissubgroup}.user_msa.fasta
cd subgroups/${thissubgroup}

source activate r_microbiome_env
Rscript $scriptdir/cospec_analysis/02_cospec_analysis.R
