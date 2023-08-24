#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_04a_metabat"

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
####    METABAT     #######
###########################
module load miniconda3

source activate metagenome_env

#for s in ${local_samples}; do
#ava=$(grep $s $workfolder/samples_available.txt)

if [ -e "$workfolder/samples/${s}/binning/${s}.metabat2.contigs_to_bin.tsv" ]; then exit; fi
echo $s
mkdir -p $workfolder/samples/${s}
mkdir -p metabat
metabat2 -i $datafolder/${s}/$s.spades_contigs.fasta -a $workfolder/samples/${s}/${s}.depth.txt -o ${s}_metabat/${s}.metabat2 -t ${SLURM_CPUS_PER_TASK} -m 2000
for fa in $(ls ${s}_metabat/${s}.metabat2.*.fa); do mv $fa ${fa}sta; done

grep '>' ${s}_metabat/*fasta | cut -d '/' -f 2 | sed 's/:>/\ /' > $workfolder/samples/${s}/binning/${s}.metabat2.contigs_to_bin.tsv
