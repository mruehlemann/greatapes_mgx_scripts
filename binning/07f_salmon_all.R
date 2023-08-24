library(tidyverse)

allsamples = list.files("samples")
#allsamples = allsamples[!allsamples %in% c("H07680-L1", "H07609-L1", "G02618-L1","H07610-L1","H07632-L1","H07664-L1","H07690-L1")]
#allsamples = allsamples[grepl("^F",allsamples)==F]
alltax = lapply(allsamples, function(x){tfile=paste0("samples/",x,"/",x,".final_tax.tsv"); if(file.exists(tfile)) read.table(tfile, head=T, stringsAsFactors=F)}) %>% do.call("bind_rows",.) %>% data.frame(stringsAsFactors=F)

for(lv in unique(alltax$level)){
	print(lv)
	alltax %>% filter(level == lv) %>% reshape2::dcast(sample ~ tax, value.var="tpm", fill=0) %>% data.frame(row.names=1) %>% write.table(., paste0("allgroups/abundance_tables/GreatApes_final_table_",lv,".tsv"), sep="\t")
}
