overlap_human_ms <-function(CpG_list_ms, CpG_list_human, background.ms, background.human, master_df, permutation_number){
  #use "LUHMES_mC_human_to_mouse" as human background object
  
  Hit_number_ms<-length(CpG_list_ms)
  Hit_number_human<-length(CpG_list_human)
  
  ## Boot strapping (to see if hits more in feature than expected)
  bootstrap_overlap<-lapply(1:permutation_number, function(x){
    set.seed(x)
    #take a random sample of all probes or genes, same size as your hit list
    rnd_CpGs_ms <- background.ms[sample(1:nrow(background.ms),Hit_number_ms),c("site","gene")]
    rnd_genes_ms <- unique(rnd_CpGs_ms$gene[complete.cases(rnd_CpGs_ms$gene)])
    
    rnd_CpGs_human <- master_df[sample(1:nrow(master_df),Hit_number_human),]
    rnd_genes_human <- background.human[background.human$refseq_acc %in% rnd_CpGs_human$acc,"mus_external_gene"]
    rnd_genes_human <- unique(rnd_genes_human[complete.cases(rnd_genes_human)])
    
    #calculate the overlap between genes from each dataset
    overlap <- length(rnd_genes_ms[rnd_genes_ms %in% rnd_genes_human])
  })
  #bind together each permutation
  bootstrap_overlap<-do.call(rbind, bootstrap_overlap)
  
  print("Permutation P values for enrichment and depletion")
  gene_list_ms <- background.ms[background.ms$site %in% CpG_list_ms,"gene"]
  gene_list_ms <- unique(gene_list_ms[complete.cases(gene_list_ms)])
  gene_list_human <- master_df[master_df$site %in% CpG_list_human,]
  gene_list_human <- background.human[background.human$refseq_acc %in% gene_list_human$acc,"mus_external_gene"]
  gene_list_human <- unique(gene_list_human[complete.cases(gene_list_human)])
  real_overlap <- length(gene_list_ms[gene_list_ms %in% gene_list_human])
  #how many iterations of bootstrapping found MORE overlapping probes than the real number? Divide this by the number of permutations to get a p-value
  enrich_p<-length(which(bootstrap_overlap>=real_overlap))/permutation_number
  #how many iterations of bootstrapping found LESS overlapping probes than the real number? Divide this by the number of permutations to get a p-value
  depletion_p<-length(which(bootstrap_overlap<=real_overlap))/permutation_number
  print(paste("Enrichment: ", enrich_p, "; Depletion ", depletion_p,sep=""))
}