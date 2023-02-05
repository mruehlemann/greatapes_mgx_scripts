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

scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt


####################################
####    SAMPLE SELECTION     ########
####################################

## selects all samples that belong to the respective subgroup (locality = column 3 of the metadata sheet)
this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples_indices=$(awk -veth=${this_local} '{if($2==eth) print (NR-1)}' $groupingfile | tr '\n' ',' | sed 's/[,]$//')
local_samples=$(awk -veth=${this_local} '{if($2==eth) print $1}' $groupingfile)

###########################
####    drep        #######
###########################

module load miniconda3
source activate binning_env

cd $workfolder/subgroups/$this_local/
mkdir -p genome_links

for i in $local_samples; do ls $workfolder/samples/$i/cleanbins/*.fa; done | awk -F '/' '{print $0"\t"$NF}' | sed 's/[.]fa$//' >  $workfolder/subgroups/$this_local/$this_local.cleanbins_list

head -n 1 $workfolder/samples/$(ls $workfolder/samples  | head -n 1 | cut -d ' ' -f 1)/*.refined.out > ${this_local}.refined.out
for i in $local_samples; do cat $workfolder/samples/$i/$i.refined.out | grep -v "uniqueSCGs.gtdb_rel207_bac120"; done >> ${this_local}.refined.out

cd $TMPDIR
cut -f 1 $workfolder/subgroups/$this_local/$this_local.cleanbins_list > genomes_list
dRep compare -g genomes_list -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -sa 0.95 -d .
cp data_tables/Cdb.csv ${this_local}.dRep_fastANI_95.csv
python $scriptdir/binning/07a_drep_fastANI99.py
cp data_tables/Cdb_97.csv ${this_local}.dRep_fastANI_97.csv
cp data_tables/Cdb_99.csv ${this_local}.dRep_fastANI_99.csv
cp ${this_local}.dRep_fastANI_95.csv ${this_local}.dRep_fastANI_97.csv ${this_local}.dRep_fastANI_99.csv $workfolder/subgroups/$this_local
conda deactivate

cd $workfolder/subgroups/$this_local/
source activate r_microbiome_env
Rscript  $scriptdir/binning/07a_drep_get_representative.R
