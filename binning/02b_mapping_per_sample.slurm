#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=120Gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_02b_map"

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

###########################
####    MAPPING     #######
###########################

source activate metagenome_env

mkdir -p bams

#for s in ${local_samples}; do
#ava=$(grep $s $workfolder/samples_available.txt)
if [ -e "$workfolder/subgroups/${this_local}/${this_local}.${s}.depth.txt" ]; then exit; fi
echo $s
minimap2 -t ${SLURM_CPUS_PER_TASK} -N 50 -ax sr  $workfolder/subgroups/${this_local}/${this_local}.catalogue.mmi $datafolder/$s/${s}.clean.R1.fastq.gz $datafolder/$s/${s}.clean.R2.fastq.gz > $s.minimap.sam
samtools view -F 3584 -b --threads ${SLURM_CPUS_PER_TASK} $s.minimap.sam | samtools sort > bams/$s.minimap.bam
rm $s.minimap.sam
#done

jgi_summarize_bam_contig_depths --outputDepth $workfolder/subgroups/${this_local}/${this_local}.${s}.depth.txt bams/*
