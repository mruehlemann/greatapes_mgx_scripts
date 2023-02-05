#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=3-00:00
#SBATCH --job-name="mb_07b_emapper"
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME
scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

###########################
####    VAMB        #######
###########################

module load miniconda3
conda activate mgx_tools_env
btchs=$(ls $workfolder/allgroups/protein_catalogs/emapper_batches/ | head -n 1 | tail -n 1)
cp $workfolder/allgroups/protein_catalogs/emapper_batches/$btchs $btchs.emapper.hits
mkdir $TMPDIR/${btchs}_tmp
emapper.py --cpu  ${SLURM_CPUS_PER_TASK} -i ${PROJECTID}.prot95_rep_seq.fasta --tax_scope 1,2,2157,10239 --itype proteins --data_dir /work_beegfs/ikmb_repository/databases/eggnog_mapper/220425/ --output $btchs --temp_dir $TMPDIR/${btchs}_tmp --dbmem -m mmseqs --resume
cp $btchs.emapper.seed_orthologs $btchs.emapper.annotations $workfolder/allgroups/protein_catalogs/emapper_batches/
