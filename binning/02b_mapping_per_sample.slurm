#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40Gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
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

this_local=$(grep -w $s $groupingfile | cut -f 2 )

###########################
####    MAPPING     #######
###########################
module load miniconda3
source activate metagenome_env
module load samtools

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


cat $datafolder/${s}/$s.spades_contigs.fasta  | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > ${s}.sample.fasta
minimap2 -I100G -d ${s}.catalogue.mmi ${s}.sample.fasta # make index
minimap2 -t ${SLURM_CPUS_PER_TASK} -N 50 -ax sr  ${s}.catalogue.mmi $datafolder/$s/${s}.clean.R1.fastq.gz $datafolder/$s/${s}.clean.R2.fastq.gz > $s.self.minimap.sam
samtools view -F 3584 -b --threads ${SLURM_CPUS_PER_TASK} $s.self.minimap.sam | samtools sort > $datafolder/${s}/${s}_final.bam

