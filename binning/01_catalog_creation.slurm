#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=100gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_01_catalog"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

cd $workfolder

source activate metagenome_env

#ls $workfolder/assemblies > $workfolder/samples_available.txt

####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($2==eth) print $1}' $groupingfile)
cd $TMPDIR

####################################
####    CATALOG CREATION    ########
####################################

for s in ${local_samples}; do
#ava=$(grep $s $workfolder/samples_available.txt)
ava=$(ls -d $datafolder/${s})
if [ "$ava" == "" ]; then continue; fi
echo $s
cat $datafolder/${s}/$s.spades_contigs.fasta  | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > ${s}.sample.fasta
done

cat *.sample.fasta > ${this_local}.catalogue.fna

minimap2 -I100G -d ${this_local}.catalogue.mmi ${this_local}.catalogue.fna # make index

gzip ${this_local}.catalogue.fna

cp ${this_local}.catalogue.mmi ${this_local}.catalogue.fna.gz  $workfolder/subgroups/${this_local}
