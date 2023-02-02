#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_09a_v_ub"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt
###########################
####    SAMPLE SELECTION     ########
####################################

s=${all_samples[$SLURM_ARRAY_TASK_ID]}
cd $TMPDIR

this_local=$(grep -w $s $basefolder/groupings.txt | cut -f 2 )


module load miniconda3
conda activate mgx_tools_env

# cd /work_beegfs/sukmb276/Metagenomes/reference_data
# cat manara_et_al/*fa > manara_catalog.fa
# salmon index -t manara_catalog.fa -i manara_salmon_index


cd $TMPDIR

salmon quant -i /work_beegfs/sukmb276/Metagenomes/reference_data/manara_salmon_index -l IU -1 $datafolder/$s/${s}.clean.R1.fastq.gz -2 $datafolder/$s/${s}.clean.R2.fastq.gz --validateMappings -o ${s}_manara_out -p ${SLURM_CPUS_PER_TASK} --meta

manara_count=$(cut -f 5 ${s}_manara_out/quant.sf | awk 'BEGIN{i=0}{i=i+$0}END{print i}')
ga_count=$(cut -f 5 $workfolder/samples/$s/${s}.salmon.out | awk 'BEGIN{i=0}{i=i+$0}END{print i}')

echo $s $manara_count $ga_count | awk '{print $0"\t"$3/$2}' > /work_beegfs/sukmb276/Metagenomes/reference_data/map_to_ref_output/${s}.salmon_stats.out
