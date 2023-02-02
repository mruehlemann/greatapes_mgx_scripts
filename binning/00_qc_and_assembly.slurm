#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=120Gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="mb_02b_map"

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

this_local=$(grep -w $s $basefolder/groupings.txt | cut -f 2 )

###########################
####    MAPPING     #######
###########################

source activate metagenome_env


ADAPTERS="/work_beegfs/ikmb_repository/software/bbmap/38.57/resources/nextera.fa.gz"
genome_index="/work_beegfs/sukmb276/references/bbmask_references/hg19/"
READMINLEN=60

allsamples=($(ls $datafolder/*R1_001.fastq.gz | awk -F '/' '{print $NF}' | cut -d '_' -f 1))

sample=${allsamples[$SLURM_ARRAY_TASK_ID]}

task_cpus=16

left=$(ls $datafolder/${sample}_*_R1_001.fastq.gz)
right=$(ls $datafolder/${sample}_*_R2_001.fastq.gz)

cd $TMPDIR

left_trimmed=tmp_${sample}.trimmed.R1.fastq.gz
right_trimmed=tmp_${sample}.trimmed.R2.fastq.gz
unpaired_trimmed=tmp_${sample}.trimmed.RU.fastq.gz
left_nophix=tmp_${sample}.nophix.R1.fastq.gz
right_nophix=tmp_${sample}.nophix.R2.fastq.gz
unpaired_nophix=tmp_${sample}.nophix.RU.fastq.gz

###QC
### basic stats of input files (can be omitted)
#reformat.sh threads=${task_cpus} in=${left} in2=${right}

### bbduk does adaptor trimming and artifact removal
bbduk.sh stats=$bbduk_adapter_stats threads=${task_cpus} in=${left} in2=${right} out1=${left_trimmed} out2=${right_trimmed} outs=${unpaired_trimmed} ref=${ADAPTERS} ktrim=r k=23 mink=11 hdist=1 minlength=${READMINLEN} tpe tbo
bbduk.sh stats=$bbduk_artifact_stats threads=${task_cpus} in=${left_trimmed} in2=${right_trimmed} k=31 ref=artifacts,phix ordered cardinality out1=${left_nophix} out2=${right_nophix} minlength=${READMINLEN}
bbduk.sh threads=${task_cpus} in=${unpaired_trimmed}  k=31 ref=artifacts,phix ordered cardinality out1=${unpaired_nophix} minlength=${READMINLEN}


### removal of host reads (human is used as reference), see for details: http://seqanswers.com/forums/archive/index.php/t-42552.html
bbwrap.sh -Xmx23g threads=${task_cpus} minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=20 minlength=${READMINLEN} in=${left_nophix} in2=${right_nophix} path=${genome_index} outu1=${left_decon} outu2=${right_decon}

##SPADES ASSEMBLY
### assembly
spades.py --meta --pe1-1 ${sample}.clean.R1.fastq.gz --pe1-2 ${sample}.clean.R2.fastq.gz --pe1-s ${sample}.clean.RU.fastq.gz -k 21,33,55 -o ${sample}_spades_out -t ${task_cpus}

mkdir -p $datafolder/${sample}

### renaming of contigs
outcontigs=${sample}.spades_contigs.fasta
awk -v id=$sample '/^>/{print ">Spades_"id"_contig_"++i; next}{print}' ${sample}_spades_out/contigs.fasta > $datafolder/${sample}/$outcontigs

cp $(ls ${sample}_spades_out/corrected/*.R1.*) $datafolder/${sample}/${sample}.clean.R1.fastq.gz
cp $(ls ${sample}_spades_out/corrected/*.R2.*) $datafolder/${sample}/${sample}.clean.R2.fastq.gz
zcat $(ls ${sample}_spades_out/corrected/*fastq.gz | egrep '([.]RU[.])|(R_unpaired)') | gzip -c > $datafolder/${sample}/${sample}.clean.RU.fastq.gz


exit