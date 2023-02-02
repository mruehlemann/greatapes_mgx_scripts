#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_03a_cds"

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

######################################################
######### CDS CALLING AND ANNOTATION     ############
######################################################

module load miniconda3
source activate gapseq-dev

cd $TMPDIR

cp $workfolder/samples/${s}/cleanbins/*.fa .

#parallel -j $SLURM_CPUS_PER_TASK cp {} {/.}.fna ::: $thisbatch
parallel -j $SLURM_CPUS_PER_TASK gzip {} ::: $(ls *.fa)

# Reaction & Pathway prediction
parallel -j $SLURM_CPUS_PER_TASK /work_ifs/sukmb276/github/gapseq/gapseq find -b 200 -v 1 -p all -k -t Bacteria {}.fa.gz ::: $(ls *.fa.gz | cut -d '.' -f 1)

# Transporter prediction
parallel -j $SLURM_CPUS_PER_TASK /work_ifs/sukmb276/github/gapseq/gapseq find-transport -b 200 -k {}.fna.gz ::: $(ls *.fa.gz | cut -d '.' -f 1)

# Draft network reconstruction
parallel -j $SLURM_CPUS_PER_TASK /work_ifs/sukmb276/github/gapseq/gapseq draft -r {}-all-Reactions.tbl -t {}-Transporter.tbl -b Bacteria -c {}.fna.gz -p {}-all-Pathways.tbl -u 200 -l 100 ::: $(ls *.fa.gz | cut -d '.' -f 1)

# gapfill/growth medium prediction
parallel -j $SLURM_CPUS_PER_TASK /work_ifs/sukmb276/github/gapseq/gapseq medium -m {}-draft.RDS -p {}-all-Pathways.tbl -c "cpd00007:0" ::: $(ls *.fa.gz | cut -d '.' -f 1)

# gap-filling
parallel -j $SLURM_CPUS_PER_TASK /work_ifs/sukmb276/github/gapseq/gapseq fill -m {}-draft.RDS -n {}-medium.csv -c {}-rxnWeights.RDS -g {}-rxnXgenes.RDS -b 100 ::: $(ls *.fa.gz | cut -d '.' -f 1)

mkdir -p $workfolder/samples/${s}/gapseq_out

for i in $(ls *.fa.gz | cut -d '.' -f 1); do
mkdir -p $workfolder/samples/${s}/gapseq_out/${i}
cp ${i}*.csv ${i}*.RDS ${i}*.tbl $workfolder/samples/${s}/gapseq_out/${i}
done
