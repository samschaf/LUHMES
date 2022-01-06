#Out of the *overlapping genes* (WTaSyn LUHMES differentially methylated and also in RRBS), how many are differentially methylated in mice?

overlap_DMhuman_ms <-function(CpG_list_ms, background.ms, permutation_number){
  #betaResults_sub = betaResults subset to genes differentially methylated in LUHMES
  #background.ms: betaResults_sub
  #CpG_list_ms: betaResults_sub CpG hits
  #taking a random number of CpGs from betaResults_sub, how often are they DM?
  
  Hit_number_ms<-length(CpG_list_ms)
  
  ## Boot strapping (to see if more DM genes in mice than expected)
  bootstrap_DM <-lapply(1:permutation_number, function(x){
    set.seed(x)
    #take a random sample of all probes or genes, same size as your hit list
    rnd_CpGs_ms <- background.ms[sample(1:nrow(background.ms),Hit_number_ms),c("site","gene","threshold")]
    rnd_CpGs_ms <- rnd_CpGs_ms[rnd_CpGs_ms$threshold==TRUE,]
    rnd_genes_ms <- unique(rnd_CpGs_ms$gene[complete.cases(rnd_CpGs_ms$gene)])
    length(rnd_genes_ms)
  })
  #bind together each permutation
  bootstrap_DM<-do.call(rbind, bootstrap_DM)
  
  print("Permutation P values for enrichment and depletion")
  
  real_DM <- background.ms[background.ms$site %in% CpG_list_ms,"gene"]
  real_DM_gene <- unique(real_DM[complete.cases(real_DM)])
  
  #how many iterations of bootstrapping found MORE overlapping probes than the real number? Divide this by the number of permutations to get a p-value
  enrich_p<-length(which(bootstrap_DM>=real_DM))/permutation_number
  #how many iterations of bootstrapping found LESS overlapping probes than the real number? Divide this by the number of permutations to get a p-value
  depletion_p<-length(which(bootstrap_DM<=real_DM))/permutation_number
  print(paste("Enrichment: ", enrich_p, "; Depletion ", depletion_p,sep=""))
}