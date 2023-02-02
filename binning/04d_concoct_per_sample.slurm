#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_04d_concoct"

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
####    concoct       #######
#############################

source activate binning_env
module load samtools/1.9


if [ -e "$workfolder/samples/${s}/${s}.concoct.contigs_to_bin.tsv" ]; then exit; fi
echo $s
mkdir -p $workfolder/samples/${s}
cut_up_fasta.py $datafolder/${s}/$s.spades_contigs.fasta -c 10000 -o 0 --merge_last -b $s.minimap.bed > ${s}.filtered.10k.fna
cp $datafolder/${s}/${s}_final.bam $s.sorted.bam
samtools index ${s}.sorted.bam
concoct_coverage_table.py $s.minimap.bed $s.sorted.bam > ${s}.coverage_table.tsv
concoct --composition_file ${s}.filtered.10k.fna --coverage_file ${s}.coverage_table.tsv -c 1000 -r 151 -t ${SLURM_CPUS_PER_TASK} -l 2000 -s 1234 -i 500 -b ${s}

merge_cutup_clustering.py ${s}_clustering_gt2000.csv > ${s}_clustering_merged.csv

awk -F',' -v sample=$s '{if(NR>1) print sample"_concoct_bin_"$2".fasta\t"$1}'  ${s}_clustering_merged.csv > $workfolder/samples/${s}/${s}.concoct.contigs_to_bin.tsv

#done
