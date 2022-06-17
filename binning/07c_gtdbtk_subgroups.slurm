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

## selects all samples that belong to the respective subgroup (locality = column 3 of the metadata sheet)
this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples_indices=$(awk -veth=${this_local} '{if($2==eth) print (NR-1)}' $groupingfile | tr '\n' ',' | sed 's/[,]$//')
local_samples=$(awk -veth=${this_local} '{if($2==eth) print $1}' $groupingfile)

###########################
####    GTDBTK        #######
###########################

module load miniconda3
source activate gtdbtk_2.1.0
GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/databases/GTDB-TK/release207_v2

cd $TMPDIR
#cd $workfolder/allgroups

cut -f 1 $workfolder/subgroups/${this_local}/${this_local}.dRep_cluster97_representatives.tsv | xargs -I {} grep -w {} $workfolder/subgroups/${this_local}/${this_local}.cleanbins_list > batch

gtdbtk classify_wf --batchfile batch --cpus $SLURM_CPUS_PER_TASK -x fa --pplacer_cpus $SLURM_CPUS_PER_TASK --out_dir gtdbtk_r207v2_out --prefix ${PROJECTID}_${this_local}

conda deactivate

cp gtdbtk_r207v2_out/${PROJECTID}_${this_local}.*.summary.tsv $workfolder/subgroups/${this_local}
