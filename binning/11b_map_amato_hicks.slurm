#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_09a_v_ub"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt
###########################
####    SAMPLE SELECTION     ########
####################################

module load miniconda3
conda activate mgx_tools_env

# cd /work_beegfs/sukmb276/Metagenomes/reference_data
# cut -f 3 manara_meta_data_representatives.tsv | xargs -I {} cat manara_et_al/{}.fa > manara_catalog.fa
# salmon index -t manara_catalog.fa -i manara_salmon_index

allga=($(cat /work_beegfs/sukmb276/Metagenomes/reference_data/greatapes_files.tsv))
thisga=${allga[$SLURM_ARRAY_TASK_ID]}
thiss=$(echo $thisga | awk -F '/' '{print $NF}' | cut -d '_' -f 1 | cut -d '.' -f 1)


cd $TMPDIR

if [ -e "$(echo $thisga | sed 's/R1/R2/')" ] && [ "$(echo $thisga | sed 's/R1/R2/')" != "$thisga" ]; then 
echo "yes"; 
salmon quant -i $workfolder/allgroups/salmon_index -l IU -1 $thisga -2 $(echo $thisga | sed 's/R1/R2/') --validateMappings -o ${thiss}_out -p ${SLURM_CPUS_PER_TASK} --meta
salmon quant -i /work_beegfs/sukmb276/Metagenomes/reference_data/manara_salmon_index -l IU -1 $thisga -2 $(echo $thisga | sed 's/R1/R2/') --validateMappings -o ${thiss}_manara_out -p ${SLURM_CPUS_PER_TASK} --meta

else 
echo no; 
salmon quant -i $workfolder/allgroups/salmon_index -l IU -r $thisga  --validateMappings -o ${thiss}_out -p ${SLURM_CPUS_PER_TASK} --meta
salmon quant -i /work_beegfs/sukmb276/Metagenomes/reference_data/manara_salmon_index -l IU -r $thisga  --validateMappings -o ${thiss}_manara_out -p ${SLURM_CPUS_PER_TASK} --meta

fi

manara_count=$(cut -f 5 ${thiss}_manara_out/quant.sf | awk 'BEGIN{i=0}{i=i+$0}END{print i}')
ga_count=$(cut -f 5 ${thiss}_out/quant.sf | awk 'BEGIN{i=0}{i=i+$0}END{print i}')

echo $thiss $manara_count $ga_count | awk '{print $0"\t"$3/$2}' > /work_beegfs/sukmb276/Metagenomes/reference_data/map_to_ref_output/${thiss}.salmon_stats.out

