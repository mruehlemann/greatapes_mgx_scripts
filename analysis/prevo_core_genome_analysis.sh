#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80Gb
#SBATCH --time=0-04:00
#SBATCH --output=/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/log/%A_%a.out
#SBATCH --job-name="ranger"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
workfolder="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final"
source $scriptdir/00_sources.txt
#####
module load miniconda3
conda activate r_microbiome_env

all_genes=($(cut -f 2 $workfolder/analyses/host_core_genome_testset.tsv | sort | uniq))

gene=${all_genes[$SLURM_ARRAY_TASK_ID]}

mkdir -p $workfolder/analyses/core_genome_analysis/$gene
cd $workfolder/analyses/core_genome_analysis/$gene

awk -vgene=$gene '{if($(NF-2) == gene) print}' $workfolder/analyses/analysis_gene_set_complete.tsv > geneset
cut -f 4 geneset | tr -d '\"'  | /work_beegfs/sukmb276/software/bin/faSomeRecords /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/protein_catalogs/GreatApes.allprot.faa /dev/stdin $gene.faa

/work_beegfs/sukmb276/software/bin/clustalo -i $gene.faa -o $gene.aln.faa --force
/work_beegfs/sukmb276/software/bin/iqtree2 -s $gene.aln.faa -T $SLURM_CPUS_PER_TASK --alrt 1000 -B 1000 -m WAG



R --vanilla <<< '''
library(ape)
library(phangorn)
library(dplyr)
genetree = read.tree(list.files(pattern=".aln.faa.contree"))
rena=read.table("../../rename_list.tsv", head=T, stringsAsFactors=F)
genes=read.table("geneset",head=F,stringsAsFactors=F)
genes = genes %>% left_join(rena, by=c("V1"="cluster")) %>% mutate(genename=paste0(rename,"_gene1")) %>% select(geneid=V4, genename, genome=V2, genome_rename=rename)

genetree$tip.label = genes[match(genetree$tip.label, genes$geneid),"genename"]
genetree$edge.length = NULL
genetree$node.label = (genetree$node.label %>% as.numeric)
genetree$node.label = ifelse(is.na(genetree$node.label),"",genetree$node.label)
write.tree(genetree, "genetree_ranger.tre")
'''


cat ../markertree_ranger.tre genetree_ranger.tre > ranger_input.tre
../../cydab/Linux/SupplementaryPrograms/OptResolutions.linux -i ranger_input.tre -B 60 -N 100 | tee ranger_optresolutions.log
grep Optimal ranger_optresolutions.log -A 1 | tail -n 1 > $gene.ranger_optim.tre

cat ../markertree_ranger.tre $gene.ranger_optim.tre > ranger_input_optim.tre
mkdir -p recon

/work_beegfs/sukmb276/software/bin/parallel -j $SLURM_CPUS_PER_TASK ../../cydab/Linux/CorePrograms/Ranger-DTL.linux --seed {} -i ranger_input_optim.tre -o recon/rangerOutput{} ::: $(seq 1 100) | tee ranger_core.log
../../cydab/Linux/CorePrograms/AggregateRanger.linux recon/rangerOutput >> ${gene}.ranger_AggregateOutput.txt
grep -E "is: [0-9]+" recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Cost",i, c, c/i}' > ${gene}.ranger_final_stats.tsv
grep -E "Duplications: [0-9]+" recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Duplications",i, c, c/i}' >> ${gene}.ranger_final_stats.tsv
grep -E "Transfers: [0-9]+" recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Transfers",i, c, c/i}' >> ${gene}.ranger_final_stats.tsv
grep -E "Losses: [0-9]+" recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Losses",i, c, c/i}' >> ${gene}.ranger_final_stats.tsv

grep Transfer, recon/* | tr ":" "," | cut -d "," -f 1,5-6 | sed 's/ Mapping --> //' | sed 's/ Recipient --> //' | sed 's/recon\///' | awk '{print "transfer,"$0}' > $gene.events.csv
grep Speciation, recon/* | tr ":" "," | cut -d "," -f 1,5 | sed 's/ Mapping --> //'  | sed 's/recon\///' | awk '{print "speciation,"$0","}' >> $gene.events.csv

grep "Species Tree" recon/rangerOutput1 -A 1 | tail -n 1 > species_tree.nodenames.tre

exit



