#!/bin/bash
#SBATCH --reservation=sukmb412_13
#SBATCH --gres=gpu:1
#SBACTH --nodelist=medg01
#SBATCH -c 4
#SBATCH --mem=100gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_04b_vamb"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt
###########################
####    SAMPLE SELECTION     ########
####################################

cd $TMPDIR

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}

####    VAMB        #######
###########################

module load miniconda3
source activate vamb_gpu_env

for i in $(seq 0 12); do
this_local=${all_local[$i]}
vamb --cuda -p ${SLURM_CPUS_PER_TASK} --outdir $workfolder/subgroups/${this_local}/VAMB_out --fasta $workfolder/subgroups/${this_local}/${this_local}.catalogue.fna.gz --jgi $workfolder/subgroups/${this_local}/${this_local}.all_depth.tsv -o _contig_
done
