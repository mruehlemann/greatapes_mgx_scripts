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

module load miniconda3
source activate gtdbtk_2.1.0
GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/databases/GTDB-TK/release207_v2

cd $workfolder/allgroups

gtdbtk classify_wf --batchfile $workfolder/allgroups/${PROJECTID}.representatives_list --cpus $SLURM_CPUS_PER_TASK -x fa --pplacer_cpus $SLURM_CPUS_PER_TASK --full_tree --out_dir gtdbtk_r207v2_out --prefix ${PROJECTID}_fulltree

head -n 1 gtdbtk_r207v2_out/${PROJECTID}_fulltree.ar53.summary.tsv | cut -f 1-2 > ${PROJECTID}.cluster_final_tax.tsv
tail -n+2 gtdbtk_r207v2_out/${PROJECTID}_fulltree.ar53.summary.tsv | cut -f 1-2 >>  ${PROJECTID}.cluster_final_tax.tsv
tail -n+2 gtdbtk_r207v2_out/${PROJECTID}_fulltree.bac120.summary.tsv | cut -f 1-2 >>  ${PROJECTID}.cluster_final_tax.tsv

cd gtdbtk_r207v2_out/align
gtdbtk infer --msa_file ${PROJECTID}_fulltree.ar53.user_msa.fasta.gz --cpus $SLURM_CPUS_PER_TASK --out_dir ar53_infer_tree
gtdbtk infer --msa_file ${PROJECTID}_fulltree.bac120.user_msa.fasta.gz --cpus $SLURM_CPUS_PER_TASK --out_dir bac120_infer_tree

cd ..

gtdbtk root --input_tree gtdbtk_r207v2_out/align/bac120_infer_tree/gtdbtk.unrooted.tree --outgroup_taxon p__Patescibacteria --output_tree gtdbtk_r207v2_out/align/bac120_infer_tree/gtdbtk.rooted.tree --gtdbtk_classification_file GreatApes.cluster_final_tax.tsv
gtdbtk root --input_tree gtdbtk_r207v2_out/align/ar53_infer_tree/gtdbtk.unrooted.tree --outgroup_taxon p__Thermoplasmatota --output_tree gtdbtk_r207v2_out/align/ar53_infer_tree/gtdbtk.rooted.tree --gtdbtk_classification_file GreatApes.cluster_final_tax.tsv

echo "(" > ${PROJECTID}.bac120+ar53.rooted.tre
cat gtdbtk_r207v2_out/align/bac120_infer_tree/gtdbtk.rooted.tree | tr ";" "," >> ${PROJECTID}.bac120+ar53.rooted.tre
cat gtdbtk_r207v2_out/align/ar53_infer_tree/gtdbtk.rooted.tree | tr -d ";" >> ${PROJECTID}.bac120+ar53.rooted.tre
echo ")1:1;" >> ${PROJECTID}.bac120+ar53.rooted.tre

conda deactivate


for this_local in Gberingei Ggorilla HsA1 HsA2 HsG Ppan Ptschw Pttrog Ptver;do
cat $workfolder/subgroups/$this_local/$this_local.cleanbins_list
done > ${PROJECTID}.cleanbins_list

gtdbtk identify --batchfile $workfolder/allgroups/${PROJECTID}.cleanbins_list --cpus $SLURM_CPUS_PER_TASK -x fa  --out_dir ${PROJECTID}_cleanbins --prefix ${PROJECTID}_cleanbins
gtdbtk align --identify_dir ${PROJECTID}_cleanbins --out_dir ${PROJECTID}_cleanbins --prefix ${PROJECTID}_cleanbins --skip_gtdb_refs --cpus $SLURM_CPUS_PER_TASK
gtdbtk align --identify_dir ${PROJECTID}_cleanbins --out_dir ${PROJECTID}_cleanbins --prefix ${PROJECTID}_cleanbins --cpus $SLURM_CPUS_PER_TASK
gtdbtk classify --align_dir ${PROJECTID}_cleanbins --batchfile $workfolder/allgroups/${PROJECTID}.cleanbins_list --out_dir ${PROJECTID}_cleanbins --prefix ${PROJECTID}_cleanbins --cpus $SLURM_CPUS_PER_TASK --pplacer_cpus $SLURM_CPUS_PER_TASK

gunzip GreatApes_cleanbins/align/GreatApes_cleanbins.ar53.user_msa.fasta.gz
gunzip GreatApes_cleanbins/align/GreatApes_cleanbins.bac120.user_msa.fasta.gz

head -n 1 GreatApes_cleanbins/${PROJECTID}_cleanbins.ar53.summary.tsv | cut -f 1-2 > ${PROJECTID}.cleanbins_tax.tsv
tail -n+2 GreatApes_cleanbins/${PROJECTID}_cleanbins.ar53.summary.tsv | cut -f 1-2 >>  ${PROJECTID}.cleanbins_tax.tsv
tail -n+2 GreatApes_cleanbins/${PROJECTID}_cleanbins.bac120.summary.tsv | cut -f 1-2 >>  ${PROJECTID}.cleanbins_tax.tsv




#
#
# source activate mashtee_env
# cut -f 1 GMbC.representatives_list > GMbC.representatives_genomes
# mashtree --numcpus 24 --file-of-files GMbC.representatives_genomes --outtree GMbC.mashtree.tre
