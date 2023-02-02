#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_05b_refine"

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
####    MAGSCOT #######
#############################


if [ -e "$workfolder/samples/${s}/${s}.refined.out" ]; then exit; fi

cd $workfolder/samples/${s}

grep $s $workfolder/subgroups/${this_local}/VAMB_out/clusters.tsv > ${s}.vamb.contigs_to_bin.tsv

awk '{print $1"\t"$2"\tvamb"}'  ${s}.vamb.contigs_to_bin.tsv > ${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tconcoct"}'  ${s}.concoct.contigs_to_bin.tsv >> ${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tmetabat2"}'  ${s}.metabat2.contigs_to_bin.tsv >> ${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tmaxbin2"}'  ${s}.maxbin2.contigs_to_bin.tsv >> ${s}.contigs_to_bin.tsv

source activate biscoreto_env

Rscript $MAGSCOT_EXEC -i ${s}.contigs_to_bin.tsv --hmm ${s}.hmm -o ${s}

source activate binning_env
mkdir -p cleanbins
cat ${s}.refined.contig_to_bin.out | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $datafolder/${s}/$s.spades_contigs.fasta /dev/stdin --output_path cleanbins
