all_tree_dist = function(tree1, tree2) {
all_ret = c(RobFoDist=as.numeric(treedist(tree1, tree2)[1]), ClustInfDist = TreeDist::ClusteringInfoDist(tree1, tree2), QuartDiv=Quartet::QuartetDivergence(Quartet::QuartetStatus(tree1, tree2), similarity=F))
return(all_ret)}

cospec_treebased = function(all_genomes, micro_tree, host_tree, n_perm){
allperm=sapply(1:n_perm, function(y){
print(y)
subtree1_tips=all_genomes %>% group_by(long) %>% sample_n(1)
subtree1 = drop.tip(micro_tree, micro_tree$tip.label[!this_tree$tip.label %in% subtree1_tips$genome])
subtree1$tip.label = unlist(subtree1_tips[match(subtree1$tip.label, subtree1_tips$genome), "long"])
permut = sapply(1:1000, function(x){subtree_shuf=subtree1; subtree_shuf$tip.label = sample(subtree_shuf$tip.label);  all_tree_dist(host_tree, subtree_shuf)})
subtree_p = 1- (colSums((((permut - all_tree_dist(host_tree, subtree1))) %>% t) > 0) / 1000)
return(subtree_p)
})
return(allperm)
}


cospec_hommola = function(all_genomes, micro_distmat, host_distmat, n_perm){
  allperm=sapply(1:n_perm, function(FF){
  subset_1_genomes=all_genomes %>% group_by(long) %>% sample_n(1)
  subset_1 = subset_1_genomes$genome

  micro_distmat.subset1=micro_distmat[subset_1_genomes$genome,subset_1_genomes$genome]

  HP.subset1 = diag(nrow(subset_1_genomes))
  rownames(HP.subset1) = subset_1_genomes$long
  colnames(HP.subset1) = subset_1_genomes$genome
  NLinks.subset1 = sum(HP.subset1)


  ### sorting of rows and columns
  host_distmat.subset1 <- host_distmat[rownames(HP.subset1),rownames(HP.subset1)]
  micro_distmat.subset1 <- micro_distmat.subset1[colnames(HP.subset1),colnames(HP.subset1)]

  ####
  host1=vector()
  par1=vector()

  for(j in 1:nrow(HP.subset1)){
          for(i in 1:ncol(HP.subset1)){
                  if(HP.subset1[j,i]==1){host1=c(host1,j);par1=c(par1,i)}
                  }
          }

  x1=as.vector(unlist(sapply(seq_along(host1[-1]), function(x) as.matrix(host_distmat.subset1)[host1[x],host1[(x+1):length(host1)]])))
  y1=as.vector(unlist(sapply(seq_along(par1[-1]), function(x) as.matrix(micro_distmat.subset1)[par1[x],par1[(x+1):length(par1)]])))

  co1=cor(x1,y1)

  permut=sapply(1:1000, function(z){
          shuf_host1=sample(host1)
          shuf_par1=sample(par1)
          x_shuf1=as.vector(unlist(sapply(seq_along(host1[-1]), function(x) as.matrix(host_distmat.subset1)[shuf_host1[x],shuf_host1[(x+1):length(host1)]])))
          y_shuf1=as.vector(unlist(sapply(seq_along(par1[-1]), function(x) as.matrix(micro_distmat.subset1)[shuf_par1[x],shuf_par1[(x+1):length(par1)]])))
          cor(x_shuf1,y_shuf1)})

  subset_p = 1- (sum((permut - co1) < 0) / 1000)
  return(subset_p)
})
return(allperm)
}



cospec_parafit = function(all_genomes, micro_distmat, host_distmat, n_perm){
  allperm=sapply(1:n_perm, function(FF){
  subset_1_genomes=all_genomes %>% group_by(long) %>% sample_n(1)
  subset_1 = subset_1_genomes$genome

  micro_distmat.subset1=micro_distmat[subset_1_genomes$genome,subset_1_genomes$genome]

  HP.subset1 = diag(nrow(subset_1_genomes))
  rownames(HP.subset1) = subset_1_genomes$long
  colnames(HP.subset1) = subset_1_genomes$genome
  NLinks.subset1 = sum(HP.subset1)


  ### sorting of rows and columns
  host_distmat.subset1 <- host_distmat[rownames(HP.subset1),rownames(HP.subset1)]
  micro_distmat.subset1 <- micro_distmat.subset1[colnames(HP.subset1),colnames(HP.subset1)]

  pf1=tryCatch(parafit(host_distmat.subset1, micro_distmat.subset1, HP.subset1, correction="cailliez", nperm=1000,silent=T), error=function(x) return(data.frame(ParaFitGlobal=NA,p.global=NA)))
  return(pf1$p.global)
})
return(allperm)
}
