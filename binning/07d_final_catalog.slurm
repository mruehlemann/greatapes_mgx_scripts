#!/bin/bash
#SBATCH --reservation=sukmb412_13
#SBATCH --gres=gpu:1
#SBACTH --nodelist=medg01
#SBATCH -c 16
#SBATCH --mem=500gb
#SBATCH --time=3-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_04b_vamb"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME
scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

###########################
####    VAMB        #######
###########################

source activate metagenome_env
cd $TMPDIR

while read line; do
echo $line
read path genome <<< $(echo $line)
awk -v genome=$genome '/^>/ {print ">"genome"_contig_" ++i; next}{print}' $path > ${genome}.fna
done < $workfolder/allgroups/${PROJECTID}.representatives_list

cat ${PROJECTID}_cluster95_*.fna >  ${PROJECTID}.catalogue.fna

minimap2 -I100G -d ${PROJECTID}.catalogue.mmi ${PROJECTID}.catalogue.fna # make index

conda deactivate
conda activate binning_env

salmon index -t ${PROJECTID}.catalogue.fna -i ${PROJECTID}_salmon_index

gzip ${PROJECTID}.catalogue.fna

cp ${PROJECTID}.catalogue.mmi ${PROJECTID}.catalogue.fna.gz  $workfolder/subgroups/allgroups
cp -r ${PROJECTID}_salmon_index $workfolder/allgroups/salmon_index


#
#
# source activate mashtee_env
# cut -f 1 GMbC.representatives_list > GMbC.representatives_genomes
# mashtree --numcpus 24 --file-of-files GMbC.representatives_genomes --outtree GMbC.mashtree.tre
