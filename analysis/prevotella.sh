#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="ca_02"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

###########################
####    SUbGROUP selection    ########
####################################

cd $workfolder
mkdir -p prevotella_analysis
grep Prevotella allgroups/GreatApes.cleanbins_tax.tsv > prevotella_analysis/prevotella_tax.tsv

while read line; do
read genome tax <<< $(echo $line)
sample=$(echo $genome | cut -d '_' -f 1)
tail -n+2 samples/$sample/prokka_annot/$genome/${genome}.tsv | awk -v g=$genome '{print g"\t"$0}'
done < prevotella_analysis/prevotella_tax.tsv > prevotella_analysis/prevotella_prokka.tsv

cd prevotella_analysis

mkdir cyda_hmm_annot

source activate metagenome_env

while read line; do
read genome tax <<< $(echo $line)
echo $genome
sample=$(echo $genome | cut -d '_' -f 1)
hmmsearch --tblout cyda_hmm_annot/$genome.hmmout /work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/201125_analysis_new/07_gapseq/cyda/uniref90-gene_cyda.filt.aln.hmm $workfolder/samples/$sample/prokka_annot/$genome/${genome}.faa
done < prevotella_tax.tsv

cat cyda_hmm_annot/*.hmmout | grep -v "^#" > prevotella_cyda_hmm.out
