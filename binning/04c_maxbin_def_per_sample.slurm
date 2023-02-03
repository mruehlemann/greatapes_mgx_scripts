#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_04c_maxbin"

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
####    MAXBIN2    #######
###########################
module load miniconda3

source activate binning_env
module load gcc

if [ -e "$workfolder/samples/${s}/${s}.maxbin2.contigs_to_bin.tsv" ]; then exit; fi

awk '{if(NR>1) print $1"\t"$3}' $workfolder/samples/${s}/${s}.depth.txt > ${s}.abu

mkdir -p ${s}_maxbin2_out
perl /work_beegfs/sukmb276/software/MaxBin-2.2.4/run_MaxBin_prodigal.pl -contig $datafolder/${s}/$s.spades_contigs.fasta -prodigal $workfolder/samples/${s}/${s}.prodigal.faa -out ${s}_maxbin2_out/${s}.maxbin2 -abund ${s}.abu -thread ${SLURM_CPUS_PER_TASK} -min_contig_length 2000

grep '>' ${s}_maxbin2_out/*fasta | cut -d '/' -f 2 | sed 's/:>/\ /' > $workfolder/samples/${s}/${s}.maxbin2.contigs_to_bin.tsv
rm -r binning/maxbin_multi_out/
