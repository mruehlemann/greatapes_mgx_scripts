#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_08a"

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

allsubgroups=($(cut -f 1 groups_all_stats_candidates.tsv | tr -d '\"'))
thissubgroup=${allsubgroups[$SLURM_ARRAY_TASK_ID]}

genomes=$(awk -v tg=$thissubgroup '{if($1 == tg){print $3}}' groups_all_candidates.tsv)

mkdir -p subgroups/${thissubgroup}
cd subgroups/${thissubgroup}


echo ${genomes[@]} | tr ' ' "\n" | ~/Isilon/software/bin/faSomeRecords /work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/GreatApes_cleanbins/align/GreatApes_cleanbins.bac120.user_msa.fasta /dev/stdin ${thissubgroup}.user_msa.fasta


ls fasta/ | cut -d '_' -f 2-4 | ~/Isilon/software/bin/faSomeRecords /home/sukmb276/Isilon/Metagenomes/projects/ApesComplete/output/201125_analysis_new/02c_Taxonomy_GTDBr202/gtdbtk_output/all.bac120.user_msa.fasta /dev/stdin ${thisbin}.bac120.user_msa.fasta
if [ "$(wc -l ${thisbin}.bac120.user_msa.fasta | cut -f 1 -d ' ')" == "0" ]; then
ls fasta/ | cut -d '_' -f 2-4 | ~/Isilon/software/bin/faSomeRecords /home/sukmb276/Isilon/Metagenomes/projects/ApesComplete/output/201125_analysis_new/02c_Taxonomy_GTDBr202/gtdbtk_output/all.ar122.user_msa.fasta /dev/stdin ${thisbin}.bac120.user_msa.fasta
fi

source activate r_microbiome_env
