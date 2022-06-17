#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_04b_vamb"


###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME
scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

####################################
####    SAMPLE SELECTION     ########
####################################

mkdir -p $workfolder/allgroups

###########################
####    DREP all        #######
###########################

module load miniconda3
source activate binning_env

head -n 1 $workfolder/samples/$(ls $workfolder/samples | head -n 1)/*.refined.out > $workfolder/allgroups/GMbC.rep97.refined.out

for this_local in ${all_local[@]};do

cd $workfolder/subgroups/$this_local/

cut -f 1 ${this_local}.dRep_cluster97_representatives.tsv  | xargs -I {} grep {} ${this_local}.cleanbins_list >>  $workfolder/allgroups/${PROJECTID}.rep97.cleanbins_list
cut -f 1 ${this_local}.dRep_cluster97_representatives.tsv  | xargs -I {} grep {} ${this_local}.refined.out  >>  $workfolder/allgroups/${PROJECTID}.rep97.refined.out

done


cd $TMPDIR
cut -f 1 $workfolder/allgroups/${PROJECTID}.rep97.cleanbins_list | sed '/^$/d' > genomes_list
dRep compare -g genomes_list -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -sa 0.95 -d .
cp data_tables/Cdb.csv $workfolder/allgroups/${PROJECTID}.rep97.dRep_95.csv
$PYTHON $scriptdir/binning/07a_drep_fastANI99.py
cp data_tables/Cdb_97.csv $workfolder/allgroups/${PROJECTID}.rep97.dRep_97.csv
cp data_tables/Cdb_99.csv $workfolder/allgroups/${PROJECTID}.rep97.dRep_99.csv
conda deactivate

cd /$workfolder/allgroups/

module load R
Rscript $scriptdir/binning/binning/07b_drep_all_get_representative.R


while read line; do
read genome cluster_95 domain score cluster_95_final <<< $(echo $line)
grep -w $genome ${PROJECTID}.rep97.cleanbins_list | awk -v cl=$cluster_95_final '{print $1"\t"cl}'
done < ${PROJECTID}.dRep_cluster95_representatives.tsv | tee ${PROJECTID}.representatives_list
