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

source activate PROKKA

cd $workfolder/samples/${s}
mkdir prokka_annot
for cb in $(ls cleanbins | cut -d '.' -f 1)
read domain phylum class order family genus species <<< $(grep $cb $workfolder/allgroups/GreatApes.cleanbins_tax.tsv | cut -f 2 | tr ' ' '_' | tr ';' ' ')
prokka --cpus $SLURM_CPUS_PER_TASK --outdir prokka_annot/$cb/ --locustag ${cb}_CDS --genus $genus --species $species --prefix ${cb} cleanbins/${cb}.fa
done
