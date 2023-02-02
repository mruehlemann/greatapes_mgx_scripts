conda activate r_microbiome_env

Rscript gyrb_tax_annot.R
cd /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/analyses/gyrb

cat gyrb_prot_to_genome.tsv  | cut -f 3 | /work_beegfs/sukmb276/software/bin/faSomeRecords /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/protein_catalogs/GreatApes.allprot.faa /dev/stdin GreatApes.gyrb.faa
grep f__Bifidobacteriaceae gyrb_prot_to_genome.tsv | cut -f 3 | /work_beegfs/sukmb276/software/bin/faSomeRecords GreatApes.gyrb.faa /dev/stdin GreatApes.bifido.gyrb.faa
grep c__Coriobacteriia gyrb_prot_to_genome.tsv | cut -f 3 | /work_beegfs/sukmb276/software/bin/faSomeRecords GreatApes.gyrb.faa /dev/stdin GreatApes.corio.gyrb.faa

cat GreatApes.bifido.gyrb.faa | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' > GreatApes.bifido.gyrb.filt.faa
cat GreatApes.corio.gyrb.faa | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' > GreatApes.corio.gyrb.filt.faa
head -n 2 GreatApes.corio.gyrb.filt.faa | tail -n 1 | awk 'BEGIN{print ">outgroup"}{print}' >  GreatApes.bifido.gyrb.outgroup.faa
cat GreatApes.bifido.gyrb.filt.faa GreatApes.bifido.gyrb.outgroup.faa > GreatApes.bifido.gyrb.filt+og.faa

/work_beegfs/sukmb276/software/bin/clustalo -i GreatApes.bifido.gyrb.filt+og.faa -o GreatApes.bifido.gyrb.filt.aln.faa --force
/work_beegfs/sukmb276/software/bin/FastTree -wag GreatApes.bifido.gyrb.filt.aln.faa > GreatApes.bifido.gyrb.filt.aln.tre
 /work_beegfs/sukmb276/software/bin/iqtree2 -s GreatApes.bifido.gyrb.filt.aln.faa -B 1000


grep f__Bifidobacteriaceae gyrb_prot_to_genome.tsv  | cut -f 2 | /work_beegfs/sukmb276/software/bin/faSomeRecords /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/GreatApes_cleanbins.bac120.user_msa.fasta.gz /dev/stdin GreatApes.bac120.user_msa.bifido.fasta
cp /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/fam_trees/f__Bifidobacteriaceae/f__Bifidobacteriaceae.faa.treefile .

#### Bacteroidaceae

grep f__Bacteroidaceae gyrb_prot_to_genome.tsv | cut -f 3 | /work_beegfs/sukmb276/software/bin/faSomeRecords GreatApes.gyrb.faa /dev/stdin GreatApes.f_bact.gyrb.faa
grep o__Bacteroidales gyrb_prot_to_genome.tsv | grep -v f__Bacteroidaceae | cut -f 3 | /work_beegfs/sukmb276/software/bin/faSomeRecords GreatApes.gyrb.faa /dev/stdin GreatApes.o_bact.gyrb.faa

cat GreatApes.f_bact.gyrb.faa | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' > GreatApes.f_bact.gyrb.filt.faa
cat GreatApes.o_bact.gyrb.faa | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 100 {print ">"$0}' > GreatApes.o_bact.gyrb.filt.faa
head -n 2 GreatApes.o_bact.gyrb.filt.faa | tail -n 1 | awk 'BEGIN{print ">outgroup"}{print}' >  GreatApes.f_bact.gyrb.outgroup.faa
cat GreatApes.f_bact.gyrb.filt.faa GreatApes.f_bact.gyrb.outgroup.faa > GreatApes.f_bact.gyrb.filt+og.faa

/work_beegfs/sukmb276/software/bin/clustalo -i GreatApes.f_bact.gyrb.filt+og.faa -o GreatApes.f_bact.gyrb.filt.aln.faa --force
/work_beegfs/sukmb276/software/bin/FastTree -wag GreatApes.f_bact.gyrb.filt.aln.faa > GreatApes.f_bact.gyrb.filt.aln.tre
 /work_beegfs/sukmb276/software/bin/iqtree2 -s GreatApes.f_bact.gyrb.filt.aln.faa -B 1000

grep f__Bacteroidaceae gyrb_prot_to_genome.tsv  | cut -f 2 | /work_beegfs/sukmb276/software/bin/faSomeRecords /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/GreatApes_cleanbins.bac120.user_msa.fasta.gz /dev/stdin GreatApes.bac120.user_msa.bifido.fasta
cp /work_beegfs/sukmb276/Metagenomes/projects/ApesComplete/output/220616_analysis_final/allgroups/fam_trees/f__Bacteroidaceae/f__Bacteroidaceae.faa.treefile .


#### moeller ref
mkdir sp
mv ../gyrb_bifido_moeller.fna .
awk '/^>/ {if(x>0) close(outname); x++; outname=sprintf("_%d.fa",x); print > outname;next;} {if(x>0) print >> outname;}' gyrb_bifido_moeller.fna

for f in _*.fa
do
    sixpack -sequence $f -outfile $f.sixpack.out -outseq $f.sixpack.fa -orfminsize 80 -width 200
done

cat *.sixpack.fa | awk -v RS=">" -v FS="\n" -v ORS="\n" -v OFS="" '$0 {$1=">"$1"\n"; print}' | cut -f 1 -d " " | sed 's/\t/\n/' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 80 {print ">"$0}' > bigfasta.fasta

cd ..
cat GreatApes.bifido.gyrb.filt.faa sp/bigfasta.fasta > aaa.faa
/work_beegfs/sukmb276/software/bin/clustalo -i aaa.faa -o aaa.aln.faa --force
/work_beegfs/sukmb276/software/bin/FastTreeMP -wag aaa.aln.faa > aaa.aln.tre

#####
mkdir btra
cp GreatApes.bifido.gyrb.filt+og.faa btra
cd btra

awk '/^>/ {if(x>0) close(outname); x++; outname=sprintf("_%d.fa",x); print > outname;next;} {if(x>0) print >> outname;}' GreatApes.bifido.gyrb.filt+og.faa

for f in _*.fa
do
    backtranambig -sequence $f -outfile $f.backtranseq.fa -osformat2 fasta 
done

cat *.backtranseq.fa | awk -v RS=">" -v FS="\n" -v ORS="\n" -v OFS="" '$0 {$1=">"$1"\n"; print}' | cut -f 1 -d " " | sed 's/\t/\n/' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 80 {print ">"$0}' > bigfasta.fasta

cd ..
cat sp/gyrb_bifido_moeller.fna btra/bigfasta.fasta > aaa.fna
muscle -super5 aaa.fna -output aaa.aln.fna
/work_beegfs/sukmb276/software/bin/FastTreeMP -gtr -nt aaa.aln.fna > aaa.aln.fna.tre
