#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete//log/%A_%a.out
#SBATCH --job-name="mb_07e_quant"

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


#############################
####    SALMON QUANT  #######
#############################

if [ -e "$workfolder/samples/${s}/${s}.salmon.out" ]; then exit; fi

cd $TMPDIR
module load miniconda3
source activate salmon_env

salmon quant -i $workfolder/allgroups/salmon_index -l IU -1 $datafolder/$s/${s}.clean.R1.fastq.gz -2 $datafolder/$s/${s}.clean.R2.fastq.gz --validateMappings -o ${s}_out -p ${SLURM_CPUS_PER_TASK} --meta

cp ${s}_out/quant.sf $workfolder/samples/${s}/${s}.salmon.out
cp ${s}_out/logs/salmon_quant.log $workfolder/samples/${s}/${s}.salmon.log
rm -r ${s}_out

cd $workfolder/samples/${s}

source activate r_microbiome_env

Rscript $scriptdir/binning/07e_salmon_tax_profile.R ${s}.salmon.out

grep "Mapping rate" $s.salmon.log | awk -v s=$s '{print s"\tmapping_rate\t"$NF}' | tr -d '%' > $s.salmon_stats.tsv
grep "total reads in the equivalence classes" $s.salmon.log | awk -v s=$s '{print s"\tmapped_reads\t"$(NF-6)}' | tr -d ',' >> $s.salmon_stats.tsv




