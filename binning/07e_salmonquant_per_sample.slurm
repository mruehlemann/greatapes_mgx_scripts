#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_07e_quant"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt
###########################
####    SAMPLE SELECTION     ########
####################################

s=${all_samples[$SLURM_ARRAY_TASK_ID]}
cd $TMPDIR

this_local=$(grep -w $s $basefolder/groupings.txt | cut -f 2 )


#############################
####    SALMON QUANT  #######
#############################

if [ -e "$workfolder/samples/${s}/${s}.salmon.out" ]; then exit; fi

cd $TMPDIR

source activate binning_env

salmon quant -i $workfolder/allgroups/salmon_index -l IU -1 $datafolder/$s/${s}.clean.R1.fastq.gz -2 $datafolder/$s/${s}.clean.R2.fastq.gz --validateMappings -o out -p ${SLURM_CPUS_PER_TASK} --meta

cp out/quant.sf $workfolder/samples/${s}/${s}.salmon.out

cd $workfolder/samples/${s}

module load R
Rscript  $scriptdir/binning/07e_salmon_tax_profile.R
