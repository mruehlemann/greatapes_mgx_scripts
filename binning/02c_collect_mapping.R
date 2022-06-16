library(tidyverse)

files=list.files()
samples=files[grep(".depth.txt", files)]
base=strsplit(samples[1], split="[.]")[[1]][1]

all_abu = sapply(samples, function(s) read.table(s, head=T, stringsAsFactors=F), simplify=F)

for(i in 1:length(all_abu)){if(i==1){all_abu_frame=all_abu[[1]]; next}; all_abu_frame=all_abu_frame %>% left_join(all_abu[[i]] %>% select(-totalAvgDepth))}

all_abu_frame$totalAvgDepth = all_abu_frame[,grep(".bam$", colnames(all_abu_frame))] %>% rowSums()

print(head(all_abu_frame))

write.table(all_abu_frame, paste0(base,".all_depth.tsv"), sep="\t", quote=F, row.names=F)
