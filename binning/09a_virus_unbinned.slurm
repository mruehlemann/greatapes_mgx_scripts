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

if [ -e "$workfolder/samples/${s}/virus_detection/unbinned_contigs_new/${s}.checkv.out" ]; then exit; fi

mkdir phamb_discovery

### Identification of short and binned contigs for exclusion
grep '>' $workfolder/samples/${s}/cleanbins/*fa | cut -d '>' -f 2 | cut -d ' ' -f 1 > phamb_discovery/bact_contigs.txt
cat $datafolder/${s}/$s.spades_contigs.fasta | cut -d ' ' -f 1 | awk '/^>/ {print "\n", $1; next}{printf $1}'| awk 'BEGIN {RS = ">" ; ORS = ""} length($2) < 2000 {print $1"\n"}' > phamb_discovery/short_contigs.txt
cat phamb_discovery/bact_contigs.txt phamb_discovery/short_contigs.txt | cut -d ' ' -f 1 > phamb_discovery/bad_contigs.txt

/work_ifs/sukmb276/software/bin/faSomeRecords $datafolder/${s}/$s.spades_contigs.fasta -exclude phamb_discovery/bad_contigs.txt /dev/stdout | cut -d ' ' -f 1 > phamb_discovery/contigs.fna


### Gene calls from unbinned contigs are retained for annotation
grep '>' phamb_discovery/contigs.fna | tr -d '>' | cut -f 1 -d ' ' > phamb_discovery/good_contigs.txt
#cat phamb_discovery/bad_contigs.txt | cut -f 1 | awk '{for(i=1;i<=1000;++i) print $1"_"i}' | /work_ifs/sukmb276/software/bin/faSomeRecords $workfolder/samples/${s}/${s}.prodigal.faa -exclude /dev/stdin phamb_discovery/proteins.faa
cat phamb_discovery/good_contigs.txt | cut -f 1 | awk '{for(i=1;i<=1000;++i) print $1"_"i}' | /work_ifs/sukmb276/software/bin/faSomeRecords $workfolder/samples/${s}/${s}.prodigal.faa /dev/stdin phamb_discovery/proteins.faa

#### annotation of contigs using DeepVirFinder
cd phamb_discovery
mkdir -p annotations
python3 $PHAMB_DIR/DeepVirFinder/dvf.py -i contigs.fna -o DVF -l 2000
mv DVF/contigs.fna_gt2000bp_dvfpred.txt annotations/all.DVF.predictions.txt

### annotation of gene calls (proteins sequences) with bacterial and viral ortholog databases (VOG annotation is slow because of the large database size)
hmmsearch --cpu ${SLURM_CPUS_PER_TASK} -E 1.0e-05 -o output.txt --tblout annotations/all.hmmMiComplete105.tbl $PHAMB_DIR/databases/Bact105.hmm proteins.faa
hmmsearch --cpu ${SLURM_CPUS_PER_TASK} -E 1.0e-05 -o output.txt --tblout annotations/all.hmmVOG.tbl $PHAMB_DIR/databases/AllVOG.hmm proteins.faa

### reformatting of contig to bin files from differnt binners and only retaining bins with the so-far unbinned contigs
#cat $workfolder/samples/${s}/${s}.concoct.contigs_to_bin.tsv | cut -f 1-2 | sed 's/[.]fasta//' | grep -Fwf good_contigs.txt > concoct_clusters.tsv
cat $workfolder/samples/${s}/${s}.metabat2.contigs_to_bin.tsv | awk '{print $1"\t"$2}' | sed 's/[.]fasta//' | grep -Fwf good_contigs.txt > metabat_clusters.tsv
cat $workfolder/samples/${s}/${s}.vamb.contigs_to_bin.tsv | awk '{print "VAMB"$1"\t"$2}' | sed 's/[.]fa//' | grep -Fwf good_contigs.txt > vamb_clusters.tsv

### link to the phamb workflow
ln -s $PHAMB_DIR/phamb/workflows/mag_annotation mag_annotation
#
# ### run PhAMB for the different binners; has to be done separately as every contig can only occure once, otherwise phamb throws errors
gzip contigs.fna
python $PHAMB_DIR/phamb/workflows/mag_annotation/scripts/run_RF.py contigs.fna.gz metabat_clusters.tsv annotations resultsdir_metabat
# python $PHAMB_DIR/phamb/workflows/mag_annotation/scripts/run_RF.py contigs.fna.gz concoct_clusters.tsv annotations resultsdir_concoct
python $PHAMB_DIR/phamb/workflows/mag_annotation/scripts/run_RF.py contigs.fna.gz vamb_clusters.tsv annotations resultsdir_vamb
#
# ### copy output to
mkdir -p $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output
# cp resultsdir_concoct/vambbins_aggregated_annotation.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/concoct_aggregated_annotation.txt
# cp resultsdir_concoct/vambbins_RF_predictions.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/concoct_RF_predictions.txt
cp resultsdir_metabat/vambbins_aggregated_annotation.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/metabat_aggregated_annotation.txt
cp resultsdir_metabat/vambbins_RF_predictions.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/metabat_RF_predictions.txt
cp resultsdir_vamb/vambbins_aggregated_annotation.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/vamb_aggregated_annotation.txt
cp resultsdir_vamb/vambbins_RF_predictions.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/vamb_RF_predictions.txt
#

grep viral $workfolder/samples/$s/virus_detection/unbinned_contigs/phamb_output/metabat_RF_predictions.txt | cut -f 1 | xargs -I {} grep -w {} metabat_clusters.tsv  > viral_candidate_clusters.tsv
grep viral $workfolder/samples/$s/virus_detection/unbinned_contigs/phamb_output/vamb_RF_predictions.txt | cut -f 1 | xargs -I {} grep -w {} vamb_clusters.tsv >> viral_candidate_clusters.tsv
#grep viral $workfolder/samples/$s/virus_detection/unbinned_contigs/phamb_output/concoct_RF_predictions.txt | cut -f 1 | xargs -I {} grep -w {} concoct_clusters.tsv  >> viral_candidate_clusters.tsv

cut -f 2 viral_candidate_clusters.tsv | sort | uniq | /work_ifs/sukmb276/software/bin/faSomeRecords contigs.fna.gz /dev/stdin viral_candidates_split.fna

checkv contamination viral_candidates_split.fna checkv_split_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv completeness viral_candidates_split.fna checkv_split_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv complete_genomes viral_candidates_split.fna checkv_split_out
checkv quality_summary viral_candidates_split.fna checkv_split_out

awk '{cont=0;if($3!=0){cont=$5/($3+0.0000001)}; if($4>0 || cont<0.49999) print $0,cont}' checkv_split_out/contamination.tsv | cut -f 1 > good_viral_contigs.txt

#cat concoct_clusters.tsv | cut -f 1-2 | sed 's/[.]fasta//' | grep -Fwf good_viral_contigs.txt > concoct_viral_clusters.tsv
cat metabat_clusters.tsv | awk '{print $1"\t"$2}' | sed 's/[.]fasta//' | grep -Fwf good_viral_contigs.txt > metabat_viral_clusters.tsv
cat vamb_clusters.tsv | awk '{print $1"\t"$2}' | sed 's/[.]fa//' | grep -Fwf good_viral_contigs.txt > vamb_viral_clusters.tsv

source activate binning_env
mkdir -p candidate_bins
gunzip contigs.fna.gz
cat vamb_viral_clusters.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | extract_fasta_bins.py contigs.fna /dev/stdin --output_path candidate_bins
#cat concoct_viral_clusters.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | extract_fasta_bins.py contigs.fna /dev/stdin --output_path candidate_bins
cat metabat_viral_clusters.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | extract_fasta_bins.py contigs.fna /dev/stdin --output_path candidate_bins

for i in $(ls candidate_bins | sed 's/[.]fa$//'); do echo $i $(cat candidate_bins/$i.fa | grep -v '>' | tr -d '\n' | wc -c); done | awk '{if($2>=10000 && $2<500000) print}' > keep.tsv

for i in $(cut -d ' ' -f 1 keep.tsv); do
echo ">"$i
cat candidate_bins/$i.fa | grep -v '>' | tr -d '\n'
echo -e "\n"
done > viral_candidates.fna

source activate phamb-1.0.2_env

checkv contamination viral_candidates.fna checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv completeness viral_candidates.fna checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv complete_genomes viral_candidates.fna checkv_out
checkv quality_summary viral_candidates.fna checkv_out

mkdir $workfolder/samples/${s}/virus_detection/unbinned_contigs_new

cp checkv_out/quality_summary.tsv $workfolder/samples/${s}/virus_detection/unbinned_contigs_new/${s}.checkv.out
cat checkv_out/*.fna > $workfolder/samples/${s}/virus_detection/unbinned_contigs_new/${s}.checkv.candidates.fna
