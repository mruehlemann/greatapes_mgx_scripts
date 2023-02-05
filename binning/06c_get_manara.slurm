#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_beegfs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_08a"

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

source activate metagenome_env
s="ManaraNHP"

cd $workfolder/subgroups
ln -s /work_beegfs/sukmb276/Metagenomes/reference_data/manara_et_al_sgb_new $s
cd species_rep_genomes
/work_beegfs/sukmb276/software/bin/parallel -j 24 prodigal -i {}.fna -a ../species_rep_prodigal/{}.faa ::: $(ls *.fna | cut -d . -f 1)
grep '^>' *.fna | sed 's/[.]fna[:][>]/\t/' | awk '{print $0"\tManaraNHP"}' > ../${s}.contigs_to_bin.tsv

cd ..


#mkdir species_rep_proteins/
#cd species_rep_proteins/
#~/Isilon/software/bin/parallel -j 10 curl -O {} ::: $(awk -F"\t" '{if($1==$14) print $NF}' ../genomes-all_metadata.tsv | sed 's/[.]gff[.]gz/.fna/' | awk -F'/' '{print "ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/species_catalogue/"$11"/"$12"/genome/"$12".faa"}')
#cd ..

cat species_rep_prodigal/*.faa | cut -d ' ' -f 1 > ${s}.prodigal.faa

mkdir hmm
source activate MAGScoT_env

hmmsearch -o hmm/${s}.hmm.tigr.out --tblout hmm/${s}.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} //work_beegfs/sukmb276/github/MAGScoT/hmm/gtdbtk_rel207_tigrfam.hmm ${s}.prodigal.faa
#pfam_scan.pl -cpu ${SLURM_CPUS_PER_TASK} -outfile ${sample}.hmm.pfam.hit.out -fasta $outdir/$sample/${sample}.prodigal.faa -dir /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/
hmmsearch -o hmm/${s}.hmm.pfam.out --tblout hmm/${s}.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_beegfs/sukmb276/github/MAGScoT/hmm/gtdbtk_rel207_Pfam-A.hmm ${s}.prodigal.faa

cat hmm/${s}.hmm.tigr.hit.out | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > hmm/${s}.tigr
cat hmm/${s}.hmm.pfam.hit.out | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > hmm/${s}.pfam
cat hmm/${s}.pfam hmm/${s}.tigr > ${s}.hmm

Rscript /work_beegfs/sukmb276/github/MAGScoT/MAGScoT.R -i ${s}.contigs_to_bin.tsv --hmm ${s}.hmm --score_only -o $s
ls $PWD/species_rep_genomes/*fna | awk -F'/' '{print $0"\t"$NF}' | sed 's/[.]fna$//' > ${s}.cleanbins_list

awk 'BEGIN{print "genome\tcluster_97\tdomain\tscore\tcluster_97_final"};{if(NR==1) next;if($16==$8){tax="Bacteria"}else{tax="Archaea"}; print $1"\t"$1"\t"tax"\t"$16"\t"$1}' ManaraNHP.scores.out | tr -d '\"' >$s.dRep_cluster97_representatives.tsv
awk 'BEGIN{print "genome\tcluster_95\tcluster_97\tcluster_99\tdomain\tscore\tcluster_95_final\tcluster_97_final\tcluster_99_final"};{if(NR==1) next;if($16==$8){tax="Bacteria"}else{tax="Archaea"}; print $1"\t"$1"\t"$1"\t"$1"\t"tax"\t"$16"\t"$1"\t"$1"\t"$1}' ManaraNHP.scores.out | tr -d '\"' > $s.dRep_cluster.tsv
