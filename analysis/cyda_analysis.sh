###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

scriptdir="/work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/greatapes_mgx_scripts"
source $scriptdir/00_sources.txt
#####
module load miniconda3
conda activate r_microbiome_env

cd $workfolder/allgroups

Rscript gyrb_cyda_tax_annot.R

cd $workfolder/analyses/cydab

cat cyda_prot_to_genome.tsv  | cut -f 5 | /work_beegfs/sukmb276/software/bin/faSomeRecords /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/protein_catalogs/GreatApes.allprot.faa.gz /dev/stdin GreatApes.cyda.faa
cat GreatApes.cyda.faa | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' > GreatApes.cyda.filt.faa

grep -w -E "g__Prevotella|g__Paraprevotella" ../../allgroups/${PROJECTID}.cluster_final_tax.tsv | cut -f 1 \
 | xargs -I {} grep -w {} ../../allgroups/GreatApes.dRep_cluster95_representatives.tsv | cut -f 1 \
 | xargs -I {} awk -vmag={} '{if($3==mag) print}' cyda_prot_to_genome.tsv | cut -f 5 \
 | /work_beegfs/sukmb276/software/bin/faSomeRecords /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/protein_catalogs/GreatApes.allprot.faa.gz /dev/stdin GreatApes.cyda_sgb.faa
cat GreatApes.cyda_sgb.faa | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 200 {print ">"$0}' > GreatApes.cyda_sgb.filt.faa

grep -w -E "g__Prevotella|g__Paraprevotella" ../../allgroups/${PROJECTID}.cluster_final_tax.tsv | cut -f 1 | xargs -I {} grep -w {} cyda_prot_to_genome.tsv 

/work_beegfs/sukmb276/software/bin/clustalo -i GreatApes.cyda_sgb.filt.faa -o GreatApes.cyda_sgb.aln.faa --force
/work_beegfs/sukmb276/software/bin/iqtree2 -s GreatApes.cyda_sgb.aln.faa -T 10 --alrt 1000 -B 1000 -m WAG


grep -w -E "g__Prevotella|g__Paraprevotella" ../../allgroups/${PROJECTID}.cluster_final_tax.tsv | cut -f 1  | /work_beegfs/sukmb276/software/bin/faSomeRecords /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/gtdbtk_r207v2_out/align/GreatApes.bac120.user_msa.fasta.gz /dev/stdin GreatApes.bac120.user_msa.paraprevo.fasta
/work_beegfs/sukmb276/software/bin/iqtree2 -s GreatApes.bac120.user_msa.paraprevo.fasta -T 24 --alrt 1000 -B 1000 -m WAG


#treefixDTL -s GreatApes.cyda_sgb.aln.faa.contree -S smap.txt -A .fasta -o .raxmlTree.rooted -n .treefixDTL.tree -e ‘‘-m PROTGAMMAJTT’’ -V1 -l treefixDTL.log geneFamily1.raxmlTree.rooted

R --vanilla < '
library(tidyverse)
### representatives list to filter 
sgb_reps_full = read.table(file = "../../allgroups/GreatApes.representatives_list")
sgb_reps_full$genome = sapply(sgb_reps_full$V1, function(x) rev(strsplit(x, split="/")[[1]])[1]) %>% gsub("[.]f[n]*a$","",.)
sgb_reps_full = sgb_reps_full %>% select(sgb=V2, genome)

### taxonomy assigments of sgbs
tax=read.table("../../allgroups/GreatApes.cluster_final_tax.tsv", head=T, stringsAsFactors=F, sep="\t")

### cyda to genome file
cyda_genomes = read.table("cyda_prot_to_genome.tsv", head=T, stringsAsFactors = F)
### filter and rename for RANGER-DTL input 
cyda_genomes = cyda_genomes %>% filter(mag_id %in% sgb_reps_full$genome) %>% group_by(GreatApes_95_final) %>% mutate(co = seq_along(prot_id), newlab = paste0(gsub("_cluster95_","", GreatApes_95_final),"_cyda",co)) %>% ungroup() 

library(ape)
#### import of cydA gene tree; rooting on Paraprevotella outgroup; Renaming for RANGER-DTL input
cyda_prevo_tree_full = read.tree("GreatApes.cyda_sgb.aln.faa.contree")
cyda_prevo_tree_root = cyda_prevo_tree_full %>% root(., outgroup=(cyda_genomes %>% filter(genus=="g__Paraprevotella", prot_id %in% cyda_prevo_tree_full$tip.label  ) %>% pull(prot_id)), resolve.root=T)
cyda_prevo_tree_root_renam=cyda_prevo_tree_root
cyda_prevo_tree_root_renam$tip.label = cyda_genomes$newlab[match(cyda_prevo_tree_root_renam$tip.label, cyda_genomes$prot_id)]
write.tree(cyda_prevo_tree_root_renam, "GreatApes.cyda_sgb.prevo.rooted.renamed.tre")

#### Remove branch lengths for RANGER-DTL input
cyda_prevo_tree_ranger = cyda_prevo_tree_root_renam
cyda_prevo_tree_ranger$edge.length = NULL
write.tree(cyda_prevo_tree_ranger, "GreatApes.cyda_sgb.prevo.rooted.ranger.tre")
cyda_prevo_tree_ranger$node.label = (cyda_prevo_tree_ranger$node.label %>% as.numeric)
cyda_prevo_tree_ranger$node.label = ifelse(is.na(cyda_prevo_tree_ranger$node.label),"",cyda_prevo_tree_ranger$node.label)
write.tree(cyda_prevo_tree_ranger, "GreatApes.cyda_sgb.prevo.rooted.ranger_bs.tre")


### import GTDB marker gene based tree; rooting on Paraprevotella outgroup; Renaming for RANGER-DTL input
marker_prevo_tree_full = read.tree("GreatApes.bac120.user_msa.paraprevo.fasta.contree")
marker_prevo_tree_root = marker_prevo_tree_full %>% root(., outgroup=(tax %>% filter(grepl("g__Paraprevotella", classification), user_genome %in% marker_prevo_tree_full$tip.label) %>% pull(user_genome)), resolve.root=T)
marker_prevo_tree_root_renam=marker_prevo_tree_root
marker_prevo_tree_root_renam$tip.label = gsub("_cluster95_","", marker_prevo_tree_root_renam$tip.label)
write.tree(marker_prevo_tree_root_renam, "GreatApes.gtdb_sgb.prevo.rooted.renamed.tre")


#### Remove branch lengths for RANGER-DTL input
marker_prevo_tree_ranger = marker_prevo_tree_root_renam
marker_prevo_tree_ranger$edge.length = NULL
marker_prevo_tree_ranger$node.label = NULL
write.tree(marker_prevo_tree_ranger, "GreatApes.gtdb_sgb.prevo.rooted.ranger.tre")

'

cat GreatApes.gtdb_sgb.prevo.rooted.ranger.tre GreatApes.cyda_sgb.prevo.rooted.ranger_bs.tre > ranger_input.tre
./Linux/SupplementaryPrograms/OptResolutions.linux -i ranger_input.tre -B 60 | tee ranger_optresolutions.log
grep Optimal ranger_optresolutions.log -A 1 | tail -n 1 > GreatApes.cyda_sgb.prevo.rooted.ranger_optim.tre

cat GreatApes.gtdb_sgb.prevo.rooted.ranger.tre GreatApes.cyda_sgb.prevo.rooted.ranger_optim.tre > ranger_input_optim.tre
mkdir -p cyda_recon

/work_beegfs/sukmb276/software/bin/parallel -j 20 ./Linux/CorePrograms/Ranger-DTL.linux --seed {} -i ranger_input_optim.tre -o cyda_recon/rangerOutput{} ::: $(seq 1 1000) | tee ranger_core.log
./Linux/CorePrograms/AggregateRanger.linux cyda_recon/rangerOutput >> cyda_ranger_AggregateOutput.txt
grep -E "is: [0-9]+" cyda_recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Cost",i, c, c/i}' > ranger_final_stats.tsv
grep -E "Duplications: [0-9]+" cyda_recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Duplications",i, c, c/i}' >> ranger_final_stats.tsv
grep -E "Transfers: [0-9]+" cyda_recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Transfers",i, c, c/i}' >> ranger_final_stats.tsv
grep -E "Losses: [0-9]+" cyda_recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+$2}END{print "Losses",i, c, c/i}' >> ranger_final_stats.tsv
grep -E "Duplications: [0-9]+" cyda_recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+($2-2.024)^2}END{print "Duplications_sd",i, c, (c/i)^0.5}' >> ranger_final_stats.tsv
grep -E "Transfers: [0-9]+" cyda_recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+($2-38.976)^2}END{print "Transfers_sd",i, c, (c/i)^0.5}' >> ranger_final_stats.tsv
grep -E "Losses: [0-9]+" cyda_recon/rangerOutput* -o | awk 'BEGIN{i=0; c=0}{i=i+1; c=c+($2-34.024)^2}END{print "Losses_sd",i, c, (c/i)^0.5}' >> ranger_final_stats.tsv
