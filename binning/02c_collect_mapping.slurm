#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10Gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_02c"

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

##################################
####    start collection   ########
##################################

source activate metagenome_env

module load R

cd $workfolder/subgroups/${this_local}
Rscript $scriptdir/binning/02c_collect_mapping.R
sed -i "s/[.]var/-var/g" ${this_local}.all_depth.tsv
