library(dplyr)

this=rev(strsplit(getwd(),split="/")[[1]])[1]

depth = read.table(paste0(this,".depth.txt"),head=T, stringsAsFactors=F)

total_contigs = nrow(depth)
total_contigs_2k = nrow(depth %>% filter(contigLen >= 2000))

total_assembly_size = depth %>% pull(contigLen) %>% sum
total_assembly_size_2k = depth %>% filter(contigLen >= 2000) %>% pull(contigLen) %>% sum

total_mapped_bases = depth %>% mutate(mapped=contigLen * totalAvgDepth) %>% pull(mapped) %>% sum
total_mapped_bases_2k = depth %>% filter(contigLen >= 2000) %>% mutate(mapped=contigLen * totalAvgDepth) %>% pull(mapped) %>% sum

N50 = depth %>% arrange(-contigLen) %>% mutate(cumlen = cumsum(contigLen), rel_cumlen = cumlen/total_assembly_size) %>% filter(rel_cumlen < .5) %>% tail(1) %>% pull(contigLen)
N50_2k = depth  %>% filter(contigLen >= 2000) %>% arrange(-contigLen) %>% mutate(cumlen = cumsum(contigLen), rel_cumlen = cumlen/total_assembly_size_2k) %>% filter(rel_cumlen < .5) %>% tail(1) %>% pull(contigLen)

L50 = depth %>% arrange(-contigLen) %>% mutate(cumlen = cumsum(contigLen), rel_cumlen = cumlen/total_assembly_size) %>% filter(rel_cumlen < .5) %>% nrow()
L50_2k = depth  %>% filter(contigLen >= 2000) %>% arrange(-contigLen) %>% mutate(cumlen = cumsum(contigLen), rel_cumlen = cumlen/total_assembly_size_2k) %>% filter(rel_cumlen < .5) %>% nrow()

if(file.exists(paste0(this,".refined.contig_to_bin.out"))){
c2b = read.table(paste0(this,".refined.contig_to_bin.out"), head=T, stringsAsFactors=F)
}else{
c2b = data.frame(contig=c(0), binnew=NA)
}

total_contigs_in_bins = nrow(depth %>% filter(contigName %in% c2b$contig))
total_assembly_in_bins = depth %>% filter(contigName %in% c2b$contig) %>% pull(contigLen) %>% sum
total_mapped_in_bins = depth  %>% mutate(mapped=contigLen * totalAvgDepth) %>% filter(contigName %in% c2b$contig) %>% pull(mapped) %>% sum

num_bins = length(unique(na.omit(c2b$binnew)))

hmm = read.table(paste0(this,".hmm"), head=F, stringsAsFactors=F)

scg_count_max = table(hmm$V2) %>% max
scg_count_q90 = table(hmm$V2) %>% quantile(probs=.9) %>% as.numeric()

smry_data = data.frame(this,total_contigs,total_assembly_size, total_mapped_bases, N50, L50, total_contigs_2k, total_assembly_size_2k, total_mapped_bases_2k, N50_2k, L50_2k, num_bins, total_contigs_in_bins, total_assembly_in_bins, total_mapped_in_bins, scg_count_max, scg_count_q90)

write.table(smry_data, paste0(this,".genome_summary.tsv"), row.names=F, sep="\t")
