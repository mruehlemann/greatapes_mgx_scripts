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

s=${all_samples[$SLURM_ARRAY_TASK_ID]}
cd $TMPDIR

this_local=$(grep -w $s $basefolder/groupings.txt | cut -f 2 )

#############################
####   PHAMB         #######
#############################
#
# This pipleine follows the PhAMB pipeline, however a bit smarter, as contigs binned by the refinement procedure are
# not included in the processing => PhAMB works based on bins and would classify these as bacterial anyway
# The clean bins are processed separately with an approach not relying on PhAMB
# In addition, not only the VAMB input, but also CONCOCT and METABAT2 output are used as sources for potential viral bins
cd $TMPDIR

source activate phamb-1.0.2_env
PHAMB_DIR="/work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources"

ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/virus_detection/binned_contigs/${s}.checkv.out" ]; then exit; fi

cat $workfolder/samples/${s}/cleanbins/*.fa | cut -d ' ' -f 1  > allbinned.fa
grep '>' $workfolder/samples/${s}/cleanbins/*.fa | cut -d ' ' -f 1 | tr -d '>' | awk -F "/" '{print $NF}' | sed -r 's/[.]fa[:]/\t/' > contig_to_bin.tsv

python3 $PHAMB_DIR/DeepVirFinder/dvf.py -i allbinned.fa -o DVF -l 10000

awk '{if($NF <= 0.01) print $1}' DVF/allbinned.fa_gt10000bp_dvfpred.txt | xargs -I {} grep -w {} contig_to_bin.tsv > contig_to_bin.cand.tsv

# for i in $(cut -f 2 contig_to_bin.cand.tsv | uniq); do
#   egrep -w $i contig_to_bin.cand.tsv | cut -f 2  |
#   echo ">${i}_viral"
#   grep -v '>' tmp_${i}_viral.fna | tr -d "\n" | awk '{print}'
# done > ${s}.viral_candidates.fa

cut -f 2 contig_to_bin.cand.tsv  | /work_ifs/sukmb276/software/bin/faSomeRecords allbinned.fa /dev/stdin ${s}.viral_candidates.fa

checkv contamination ${s}.viral_candidates.fa checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv completeness ${s}.viral_candidates.fa checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv complete_genomes ${s}.viral_candidates.fa checkv_out
checkv quality_summary ${s}.viral_candidates.fa checkv_out

mkdir -p $workfolder/samples/${s}/virus_detection/binned_contigs

cp checkv_out/quality_summary.tsv $workfolder/samples/${s}/virus_detection/binned_contigs/${s}.checkv.out
head -n 1 checkv_out/quality_summary.tsv > $workfolder/samples/${s}/virus_detection/binned_contigs/${s}.checkv.proviruses.out
cat checkv_out/quality_summary.tsv | awk -F '\t' '{if($3=="Yes") print}' >> $workfolder/samples/${s}/virus_detection/binned_contigs/${s}.checkv.proviruses.out

touch $workfolder/samples/${s}/virus_detection/binned_contigs/${s}.checkv.proviruses.fna
if [ -f "checkv_out/proviruses.fna" ]; then
cp checkv_out/proviruses.fna $workfolder/samples/${s}/virus_detection/binned_contigs/${s}.checkv.proviruses.fna
fi
