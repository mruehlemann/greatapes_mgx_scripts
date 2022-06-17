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
scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

###########################
####    VAMB        #######
###########################

module load miniconda3
source activate gtdbtk_2.1.0
GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/databases/GTDB-TK/release207_v2

cd $workfolder/allgroups

gtdbtk classify_wf --batchfile $workfolder/allgroups/${PROJECTID}.representatives_list --cpus $SLURM_CPUS_PER_TASK -x fa --pplacer_cpus $SLURM_CPUS_PER_TASK --full_tree --out_dir gtdbtk_r207v2_out --prefix ${PROJECTID}_fulltree

cd gtdbtk_r207v2_out/align
gtdbtk infer --msa_file ${PROJECTID}_fulltree.ar53.user_msa.fasta.gz --cpus 24 --out_dir ar53_infer_tree
gtdbtk infer --msa_file ${PROJECTID}_fulltree.bac120.user_msa.fasta.gz --cpus 24 --out_dir bac120_infer_tree

conda deactivate
#
#
# source activate mashtee_env
# cut -f 1 GMbC.representatives_list > GMbC.representatives_genomes
# mashtree --numcpus 24 --file-of-files GMbC.representatives_genomes --outtree GMbC.mashtree.tre
