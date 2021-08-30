# Author: Nicole Gladish: May 2019

# Here are the arguments you can add (and their defaults): 
# PCA_Plot(PCA_Object, type = c("All", "Sparse"), nPCs = 10, MTC = F, Discrete = T, label.y_size = 12, label.x_size = 12, angle.x = 30, vjust.x = 0.5)
# Here is what some arguments do: 
# type = c("All", "Sparse") - this is to tell it whether your PCA_Object results from a sparse or regular PCA (it does not do it for you it just tells the function what you have).
# nPCs - how many PCs do you want to plot
# MTC = F - do you want the associations multiple test corrected? BH is the method I include here
# Discrete = T, do you want a continuous display of significane values or discrete?
# The rest are to size and orient your labels on the graph itself.... 
# UPDATE 14/11/2019 (YC) - changed "Loadings" to "Rotation". prcomp creates rotation, while princomp creates loadings. Both are the same. If you are running princomp, change rotation into loadings; if prcomp, use rotation. Most people use prcomp.

PCA_Plot <- function(PCA_Object, type = c("All", "Sparse"), nPCs = 10, MTC = F, Discrete = T, label.y_size = 12, label.x_size = 12, angle.x = 30, vjust.x = 0.5){

  library(viridis)
    #ord <- 1:(ncol(meta_categorical) + ncol(meta_continuous))

if (type == "All"){  
  Rotation<-as.data.frame(unclass(PCA_Object$loadings))
  Importance<-(PCA_Object$sdev^2)/sum(PCA_Object$sdev^2)
  adjust<-1-Importance[1]
  pca_adjusted<-Importance[2:length(Importance)]/adjust
  pca_df<-data.frame(adjusted_variance=pca_adjusted, PC=seq(1:length(pca_adjusted)))
  pca_df <- within(pca_df, cumulative_variance <- cumsum(adjusted_variance))
}else{
  Rotation<-as.data.frame(PCA_Object$loadings)
  Center<-as.data.frame(unclass(PCA_Object$center))
  rownames(Rotation) <- rownames(Center)
  Importance<-(PCA_Object$sdev^2)/sum(PCA_Object$sdev^2)
  adjust<-1-Importance[1]
  pca_adjusted<-Importance[2:length(Importance)]/adjust
  pca_df<-data.frame(adjusted_variance=pca_adjusted, PC=seq(1:length(pca_adjusted)))
  pca_df <- within(pca_df, cumulative_variance <- cumsum(adjusted_variance))
}

  # Allows code to work if you only have categorical or continuous data alone
  if (exists("meta_categorical") & exists("meta_continuous")){
    
    stopifnot(identical(rownames(meta_categorical), rownames(Rotation)))
    stopifnot(identical(rownames(meta_continuous), rownames(Rotation)))
    aov_PC_meta <- lapply(1:ncol(meta_categorical), function(covar) sapply(1:ncol(Rotation), 
                                                                           function(PC) summary(aov(Rotation[, PC] ~ meta_categorical[, covar]))[[1]]$"Pr(>F)"[1]))
    
    cor_PC_meta <- lapply(1:ncol(meta_continuous), function(covar) sapply(1:ncol(Rotation), 
                                                                          function(PC) (cor.test(Rotation[, PC], as.numeric(meta_continuous[, 
                                                                                                                                            covar]), alternative = "two.sided", method = "pearson", na.action = na.omit)$p.value)))
    
    names(aov_PC_meta) <- colnames(meta_categorical)
    names(cor_PC_meta) <- colnames(meta_continuous)
    aov_PC_meta <- do.call(rbind, aov_PC_meta)
    cor_PC_meta <- do.call(rbind, cor_PC_meta)
    aov_PC_meta <- rbind(aov_PC_meta, cor_PC_meta)
    aov_PC_meta <- as.data.frame(aov_PC_meta) 

    }else{
      
      if (exists("meta_categorical") & !exists("meta_continuous")){
        
        stopifnot(identical(rownames(meta_categorical), rownames(Rotation)))
        aov_PC_meta <- lapply(1:ncol(meta_categorical), function(covar) sapply(1:ncol(Rotation), 
                                                                               function(PC) summary(aov(Rotation[, PC] ~ meta_categorical[, covar]))[[1]]$"Pr(>F)"[1]))
        
        names(aov_PC_meta) <- colnames(meta_categorical)
        aov_PC_meta <- do.call(rbind, aov_PC_meta)
        aov_PC_meta <- as.data.frame(aov_PC_meta)}else{
          
          if (!exists("meta_categorical") & exists("meta_continuous")){
            stopifnot(identical(rownames(meta_continuous), rownames(Rotation)))
            cor_PC_meta <- lapply(1:ncol(meta_continuous), function(covar) sapply(1:ncol(Rotation), 
                                                                                  function(PC) (cor.test(Rotation[, PC], as.numeric(meta_continuous[, 
                                                                                                                                                    covar]), alternative = "two.sided", method = "pearson", na.action = na.omit)$p.value)))
            
            names(cor_PC_meta) <- colnames(meta_continuous)
            aov_PC_meta <- do.call(rbind, cor_PC_meta)
            aov_PC_meta <- as.data.frame(aov_PC_meta)}else{
              print("Missing both meta_continuous and meta_categorical objects")
            }}}

  aov_PC_meta_MTC <- as.data.frame(t(aov_PC_meta))
  aov_PC_meta_MTC$PCs <- rownames(aov_PC_meta_MTC)
  aov_PC_meta_MTCshort <- aov_PC_meta_MTC[1:(nPCs + 1),]
  aov_PC_meta_MTCm <- melt(aov_PC_meta_MTCshort, id.vars = "PCs")
  aov_PC_meta_MTCm$qvalues <- p.adjust(aov_PC_meta_MTCm$value, "BH")
  aov_PC_meta_MTCm$value <- NULL
  aov_PC_meta_MTC <- dcast(aov_PC_meta_MTCm, formula = factor(aov_PC_meta_MTCm$PCs, levels = unique(aov_PC_meta_MTCm$PCs)) ~ variable)
  colnames(aov_PC_meta_MTC)[1] <- "PCs"
  rownames(aov_PC_meta_MTC) <- as.character(aov_PC_meta_MTC$PCs)
  aov_PC_meta_MTC$PCs <- NULL
  aov_PC_meta_MTC <- as.data.frame(t(aov_PC_meta_MTC))
  
  
  if (MTC == T){
    aov_PC_meta_adjust<-aov_PC_meta_MTC[,2:ncol(aov_PC_meta_MTC)]
  }else{
    aov_PC_meta_adjust<-aov_PC_meta[,2:ncol(aov_PC_meta)]
  }
  
  #reshape
  colnames(aov_PC_meta_adjust) <- paste("PC", 1:ncol(aov_PC_meta_adjust), sep = "")
  avo_heat_num<-apply(aov_PC_meta_adjust,2, as.numeric)
  avo_heat<-as.data.frame(avo_heat_num)
  avo_heat <- avo_heat[,1:nPCs]
  avo_heat$meta<-rownames(aov_PC_meta_adjust)
  avo_heat_melt<-melt(avo_heat, id=c("meta"))
  
  # cluster meta data
  meta_var_order<-unique(avo_heat_melt$meta)[rev(ord)]
  avo_heat_melt$meta <- factor(avo_heat_melt$meta, levels = meta_var_order)
  
  if (MTC == T & Discrete == T){
    avo_heat_melt$Pvalue<-sapply(1:nrow(avo_heat_melt), function(x) if(avo_heat_melt$value[x]<=0.005){"<=0.005"}else{
      if(avo_heat_melt$value[x]<=0.05){"<=0.05"}else{
        if(avo_heat_melt$value[x]<=0.1){"<=0.1"}else{
          if(avo_heat_melt$value[x]<=0.2){"<=0.2"}else{">0.2"}}}})
    
    heat_legend<-ggplot(avo_heat_melt, aes(variable,meta, fill = Pvalue)) +
      geom_tile(color = "black",size=0.5) +
      theme_gray(8)+scale_fill_manual(values=c("<=0.005" = "#440154FF", "<=0.05" = "#443A83FF", "<=0.1" = "#31688EFF", "<=0.2" = "#1F9E89FF", ">0.2" = "#C1D4DD"), name = "q-value")+
      theme(axis.text = element_text(color="black"),
            axis.text.x = element_text(size = label.x_size, angle = angle.x, vjust = vjust.x),
            axis.text.y = element_text(size = label.y_size),
            axis.title = element_text(size =20),
            legend.text = element_text(size =14),
            legend.title = element_text(size =12),
            legend.position = c(1, 0.4), legend.justification = c(1,1),
            plot.margin=unit(c(0,2.25,1,1),"cm"))+
      xlab("Adjusted Principal Component")+ylab(NULL)
    
    heat<-ggplot(avo_heat_melt, aes(variable,meta, fill = Pvalue)) +
      geom_tile(color = "black",size=0.5) +
      theme_gray(8)+scale_fill_manual(values=c("<=0.005" = "#440154FF", "<=0.05" = "#443A83FF", "<=0.1" = "#31688EFF", "<=0.2" = "#1F9E89FF", ">0.2" = "#C1D4DD"))+
      theme(axis.text = element_text(color="black"),
            axis.text.x = element_text(size = label.x_size, angle = angle.x,vjust = vjust.x),
            axis.text.y = element_text(size = label.y_size),
            axis.title = element_text(size =20),
            legend.position = "none",
            plot.margin=unit(c(0,2.25,1,1),"cm"))+
      xlab("Adjusted Principal Component")+ylab(NULL)
    
    scree<-ggplot(pca_df[which(pca_df$PC<(nPCs+1)),],aes(PC,adjusted_variance))+geom_bar(stat = "identity",color="black",fill="darkgrey")+theme_bw()+
      theme(axis.text = element_text(size =18),
            axis.title = element_text(size =24),
            plot.margin=unit(c(1.25,1.6,0.2,3),"cm"))+ylab("Adjusted Variance")+
      scale_x_continuous(breaks = seq(1,nPCs,1))
  }else{
    
    if (MTC == F & Discrete == T){
      
      # color if sig
      avo_heat_melt$Pvalue<-sapply(1:nrow(avo_heat_melt), function(x) if(avo_heat_melt$value[x]<=0.00005){"<=0.00005"}else{
        if(avo_heat_melt$value[x]<=0.0005){"<=0.0005"}else{
          if(avo_heat_melt$value[x]<=0.005){"<=0.005"}else{
            if(avo_heat_melt$value[x]<=0.05){"<=0.05"}else{">0.05"}}}})

      heat_legend <-ggplot(avo_heat_melt, aes(variable,meta, fill = Pvalue)) +
        geom_tile(color = "black",size=0.5) +
        theme_gray(8)+scale_fill_manual(values=c("<=0.00005" = "#440154FF", "<=0.0005" = "#443A83FF", "<=0.005" = "#31688EFF", "<=0.05" = "#1F9E89FF", ">0.05" = "#C1D4DD"), name = "p-value")+
        theme(axis.text = element_text(color="black"),
              axis.text.x = element_text(size = label.x_size, angle = angle.x, vjust = vjust.x),
              axis.text.y = element_text(size = label.y_size),
              axis.title = element_text(size =20),
              legend.text = element_text(size =14),
              legend.title = element_text(size =12),
              legend.position = c(1, 0.4), legend.justification = c(1,1),
              plot.margin=unit(c(0,2.25,1,1),"cm"))+
        xlab("Adjusted Principal Component")+ylab(NULL)
      
      heat<-ggplot(avo_heat_melt, aes(variable,meta, fill = Pvalue)) +
        geom_tile(color = "black",size=0.5) +
        theme_gray(8)+scale_fill_manual(values=c("<=0.00005" = "#440154FF", "<=0.0005" = "#443A83FF", "<=0.005" = "#31688EFF", "<=0.05" = "#1F9E89FF", ">0.05" = "#C1D4DD"))+
        theme(axis.text = element_text(color="black"),
              axis.text.x = element_text(size = label.x_size, angle = angle.x, vjust = vjust.x),
              axis.text.y = element_text(size = label.y_size),
              axis.title = element_text(size =20), 
              legend.position = "none",
              plot.margin=unit(c(0,2.25,1,1),"cm"))+
        xlab("Adjusted Principal Component")+ylab(NULL)
      
      ## Scree plot
      scree<-ggplot(pca_df[which(pca_df$PC<(nPCs+1)),],aes(PC,adjusted_variance))+geom_bar(stat = "identity",color="black",fill="darkgrey")+theme_bw()+
        theme(axis.text = element_text(size =18),
              axis.title = element_text(size =24),
              plot.margin=unit(c(1.25,1.6,0.2,3),"cm"))+ylab("Adjusted Variance")+
        scale_x_continuous(breaks = seq(1,nPCs,1))
    }else{
      if (MTC == T & Discrete == F){
 
        heat_legend<-ggplot(avo_heat_melt, aes(variable, meta, fill = value)) +
          geom_tile(color = "black",size=0.5) +
          theme_gray(8)+scale_fill_viridis(option="magma", name = "q-value") +
          theme(axis.text = element_text(color="black"),
                axis.text.x = element_text(size = label.x_size, angle = angle.x,vjust = vjust.x),
                axis.text.y = element_text(size = label.y_size),
                axis.title = element_text(size =20),
                legend.text = element_text(size =14),
                legend.title = element_text(size =12),
                legend.position = c(1, 0.4), legend.justification = c(1,1),
                plot.margin=unit(c(0,2.25,1,1),"cm"))+
          xlab("Adjusted Principal Component")+ylab(NULL)
        
        heat<-ggplot(avo_heat_melt, aes(variable,meta, fill = value)) +
          geom_tile(color = "black",size=0.5) +
          theme_gray(8)+
          scale_fill_viridis(option="magma", name = "p-value") +
          theme(axis.text = element_text(color="black"),
                axis.text.x = element_text(size = label.x_size, angle = angle.x, vjust = vjust.x),
                axis.text.y = element_text(size = label.y_size),
                axis.title = element_text(size =20),
                legend.position = "none",
                plot.margin=unit(c(0,2.25,1,1),"cm"))+
          xlab("Adjusted Principal Component")+ylab(NULL)
        
        scree<-ggplot(pca_df[which(pca_df$PC<(nPCs+1)),],aes(PC,adjusted_variance))+geom_bar(stat = "identity",color="black",fill="darkgrey")+theme_bw()+
          theme(axis.text = element_text(size =18),
                axis.title = element_text(size =24),
                plot.margin=unit(c(1.25,1.6,0.2,3),"cm"))+ylab("Adjusted Variance")+
          scale_x_continuous(breaks = seq(1,nPCs,1))
      }else{
        if (MTC == F & Discrete == F){
   
          heat_legend<-ggplot(avo_heat_melt, aes(variable, meta, fill = value)) +
            geom_tile(color = "black",size=0.5) +
            theme_gray(8)+scale_fill_viridis(option="magma", name = "p-value") +
            theme(axis.text = element_text(color="black"),
                  axis.text.x = element_text(size = label.x_size, angle = angle.x,vjust = vjust.x),
                  axis.text.y = element_text(size = label.y_size),
                  axis.title = element_text(size =20),
                  legend.text = element_text(size =14),
                  legend.title = element_text(size =12),
                  legend.position = c(1, 0.4), legend.justification = c(1,1),
                  plot.margin=unit(c(0,2.25,1,1),"cm"))+
            xlab("Adjusted Principal Component")+ylab(NULL)
          
          heat<-ggplot(avo_heat_melt, aes(variable,meta, fill = value)) +
            geom_tile(color = "black",size=0.5) +
            theme_gray(8)+
            scale_fill_viridis(option="magma", name = "p-value") +
            theme(axis.text = element_text(color="black"),
                  axis.text.x = element_text(size = label.x_size, angle = angle.x,vjust = vjust.x),
                  axis.text.y = element_text(size = label.y_size),
                  axis.title = element_text(size =20),
                  legend.position = "none",
                  plot.margin=unit(c(0,2.25,1,1),"cm"))+
            xlab("Adjusted Principal Component")+ylab(NULL)
          
          scree<-ggplot(pca_df[which(pca_df$PC<(nPCs+1)),],aes(PC,adjusted_variance))+geom_bar(stat = "identity",color="black",fill="darkgrey")+theme_bw()+
            theme(axis.text = element_text(size =18),
                  axis.title = element_text(size =24),
                  plot.margin=unit(c(1.25,1.6,0.2,3),"cm"))+ylab("Adjusted Variance")+
            scale_x_continuous(breaks = seq(1,nPCs,1))
        }}}}
  
  list(heat_legend, heat, scree)
}