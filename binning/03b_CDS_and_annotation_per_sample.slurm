#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_03a_cds"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt

###########################
####    SAMPLE SELECTION     ########
####################################

s=${all_samples[$SLURM_ARRAY_TASK_ID]}
cd $TMPDIR

this_local=$(grep -w $s $groupingfile | cut -f 2 )

######################################################
######### CDS CALLING AND ANNOTATION     ############
######################################################

source activate metagenome_env

if [ -e "$workfolder/samples/${s}/${s}.hmm" ]; then exit; fi
echo $s
mkdir -p $workfolder/samples/${s}
mkdir -p ${s}_tmp_workfolder
cat $datafolder/${s}/$s.spades_contigs.fasta | $PARALLEL -j ${SLURM_CPUS_PER_TASK} --block 999k --recstart '>' --pipe prodigal -p meta -a ${s}_tmp_workfolder/${s}.{#}.faa -d ${s}_tmp_workfolder/${s}.{#}.ffn -o tmpfile
cat ${s}_tmp_workfolder/$s.*.faa > $workfolder/samples/${s}/${s}.prodigal.faa
cat ${s}_tmp_workfolder/$s.*.ffn > $workfolder/samples/${s}/${s}.prodigal.ffn

hmmsearch -o ${s}.hmm.tigr.out --tblout ${s}.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} //work_ifs/sukmb276/github/MAGScoT/hmm/gtdbtk_rel207_tigrfam.hmm $workfolder/samples/${s}/${s}.prodigal.faa
#pfam_scan.pl -cpu ${SLURM_CPUS_PER_TASK} -outfile ${sample}.hmm.pfam.hit.out -fasta $outdir/$sample/${sample}.prodigal.faa -dir /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/
hmmsearch -o ${s}.hmm.pfam.out --tblout ${s}.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/sukmb276/github/MAGScoT/hmm/gtdbtk_rel207_Pfam-A.hmm $workfolder/samples/${s}/${s}.prodigal.faa

cat ${s}.hmm.tigr.hit.out | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > ${s}.tigr
cat ${s}.hmm.pfam.hit.out | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > ${s}.pfam
cat ${s}.pfam ${s}.tigr > $workfolder/samples/${s}/${s}.hmm

cp $datafolder/${s}/${s}_depth.txt $workfolder/samples/${s}/${s}.depth.txt

#jgi_summarize_bam_contig_depths --outputDepth $workfolder/samples/${s}/${s}.depth.txt $datafolder/${s}/${s}_final.bam
