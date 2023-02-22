conda activate r_microbiome_env

Rscript gyrb_cyda_tax_annot.R
cd /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/gyrb

cat gyrb_prot_to_genome.tsv  | cut -f 3 | /work_beegfs/sukmb276/software/bin/faSomeRecords /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/protein_catalogs/GreatApes.allprot.faa.gz /dev/stdin GreatApes.gyrb.faa
grep f__Bifidobacteriaceae gyrb_prot_to_genome.tsv | cut -f 3 | /work_beegfs/sukmb276/software/bin/faSomeRecords GreatApes.gyrb.faa /dev/stdin GreatApes.bifido.gyrb.faa
grep c__Coriobacteriia gyrb_prot_to_genome.tsv | cut -f 3 | /work_beegfs/sukmb276/software/bin/faSomeRecords GreatApes.gyrb.faa /dev/stdin GreatApes.corio.gyrb.faa

cat GreatApes.bifido.gyrb.faa | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' > GreatApes.bifido.gyrb.filt.faa
cat GreatApes.corio.gyrb.faa | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' > GreatApes.corio.gyrb.filt.faa
head -n 2 GreatApes.corio.gyrb.filt.faa | tail -n 1 | awk 'BEGIN{print ">outgroup"}{print}' >  GreatApes.bifido.gyrb.outgroup.faa
cat GreatApes.bifido.gyrb.filt.faa GreatApes.bifido.gyrb.outgroup.faa > GreatApes.bifido.gyrb.filt+og.faa

/work_beegfs/sukmb276/software/bin/clustalo -i GreatApes.bifido.gyrb.filt+og.faa -o GreatApes.bifido.gyrb.filt.aln.faa --force
/work_beegfs/sukmb276/software/bin/FastTree -wag GreatApes.bifido.gyrb.filt.aln.faa > GreatApes.bifido.gyrb.filt.aln.tre
/work_beegfs/sukmb276/software/bin/iqtree2 -s GreatApes.bifido.gyrb.filt.aln.faa --alrt 1000 -B 1000 -m WAG


grep f__Bifidobacteriaceae gyrb_prot_to_genome.tsv  | cut -f 2 | /work_beegfs/sukmb276/software/bin/faSomeRecords /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/GreatApes_cleanbins.bac120.user_msa.fasta.gz /dev/stdin GreatApes.bac120.user_msa.bifido.fasta
cp /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/fam_trees/f__Bifidobacteriaceae/f__Bifidobacteriaceae.faa.treefile .


#### moeller ref
mkdir sp
mv ../gyrb_bifido_moeller.fna .
awk '/^>/ {if(x>0) close(outname); x++; outname=sprintf("_%d.fa",x); print > outname;next;} {if(x>0) print >> outname;}' gyrb_bifido_moeller.fna

cd ..
blastx -db GreatApes.bifido.gyrb.filt+og.faa -query sp/gyrb_bifido_moeller.fna -evalue 0.000000001 -outfmt "6 qseqid qseq" -num_alignments 1 | awk '{if(length($2)>=60) printf ">"$1"\n"$2"\n"}' > gyrb_bifido_moeller.translated.fna
cat GreatApes.bifido.gyrb.filt.faa GreatApes.bifido.gyrb.outgroup.faa gyrb_bifido_moeller.translated.fna > GreatApes.bifido.gyrb.filt+og+moeller.faa
/work_beegfs/sukmb276/software/bin/clustalo -i GreatApes.bifido.gyrb.filt+og+moeller.faa -o GreatApes.bifido.gyrb.filt+moeller.aln.faa --force
/work_beegfs/sukmb276/software/bin/iqtree2 -s GreatApes.bifido.gyrb.filt+moeller.aln.faa --alrt 1000 -B 1000 -m WAG -T 10

