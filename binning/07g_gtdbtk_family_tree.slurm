#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_03a_cds"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME
scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

###########################
####    VAMB        #######
###########################

module load miniconda3
source activate gtdbtk_2.1.0
GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/databases/GTDB-TK/release207_v2

cd $workfolder/allgroups
allfam=($(cat ${PROJECTID}.families.txt))
this_fam=${allfam[$SLURM_ARRAY_TASK_ID]}

mkdir -p fam_trees/$this_fam
grep -w $this_fam ${PROJECTID}.cluster_final_tax.tsv | cut -f 1 | faSomeRecords gtdbtk_r207v2_out/align/${PROJECTID}.bac120.user_msa.fasta.gz /dev/stdin fam_trees/$this_fam/${this_fam}.msa.faa
grep -w $this_fam ${PROJECTID}.cluster_final_tax.tsv -v | grep -v Archaea | head -n 2 | tail -n 1 | cut -f 1 | faSomeRecords gtdbtk_r207v2_out/align/${PROJECTID}.bac120.user_msa.fasta.gz /dev/stdin fam_trees/$this_fam/${this_fam}.outgroup.faa

#grep -w $this_fam ${PROJECTID}.cluster_final_tax.tsv | cut -f 1 | faSomeRecords gtdbtk_r207v2_out/align/${PROJECTID}.ar53.user_msa.fasta.gz /dev/stdin fam_trees/$this_fam/${this_fam}.msa.faa
#grep -w $this_fam ${PROJECTID}.cluster_final_tax.tsv -v | grep Archaea | head -n 2 | tail -n 1 | cut -f 1 | faSomeRecords gtdbtk_r207v2_out/align/${PROJECTID}.ar53.user_msa.fasta.gz /dev/stdin fam_trees/$this_fam/${this_fam}.outgroup.faa


cd fam_trees/$this_fam
cat ${this_fam}.outgroup.faa | grep -v '>' | awk 'BEGIN{print ">outgroup"}{print}' > $this_fam.og.faa
cat ${this_fam}.msa.faa ${this_fam}.og.faa > $this_fam.faa

iqtree2 -s ${this_fam}.faa -T $SLURM_CPUS_PER_TASK --alrt 1000 -B 1000 -m WAG
