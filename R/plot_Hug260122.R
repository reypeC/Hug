########################################################################################
########################################################################################
########################################################################################

#                     Import library

########################################################################################
########################################################################################
########################################################################################

rm(list=ls())

library(fields)
library(scatterplot3d)
library(xtable)

palette("R4")

########################################################################################
########################################################################################
########################################################################################

#                     Core function

########################################################################################
########################################################################################
########################################################################################

# Create gradient of colors
colfunc<-paste0(colorRampPalette(c("white","royalblue","yellow","red",alpha=0.5))(10),"75")[-c(9,10)]

setwd("C:/Users/creype/Nextcloud/MyDrive/Hug/simHug/Hug250721")



save_plot_statistics=function(typedata,savename=typedata,annealing=T)
{
  stats=read.table(paste0("statistics","_",typedata,".txt"),header=T)
  statsFinal=read.table(paste0("statistics_final","_",typedata,".txt"),header=F)
  colnames(statsFinal)=colnames(stats)
  
  statsComplete=rbind(stats,statsFinal)
  
  statsByPlane=split(statsComplete,factor(paste0(statsComplete[,5],statsComplete[,6])))
  
  numberPatternCooling=nrow(stats)/length(statsByPlane)
  
  energy=rep(0,nrow(statsByPlane[[1]]))
  temperature=rep(NA,nrow(statsByPlane[[1]]))
  t=1000
  c=0.9999
  iteration=0
  save=0
  while(t>0.0001)
  {
    iteration=iteration+1
    t=c*t
    if(iteration%%100==0)
    {
      save=save+1
      temperature[save+1]=t
    }
    
    
  }
  temperature[is.na(temperature)]=t
  
  for (plan in 1:length(statsByPlane)) {
    
    
    png(paste0("stats_",savename,"_plane_",plan,".png"),width = 1000,height = 200)
    par(mfrow=c(1,2))
    par(mar=c(5,6,4,1)+.1)
    plot(statsByPlane[[plan]][,1],type = "l",main = "",ylab = "n_e",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
    abline(v=numberPatternCooling,col="red",lwd=2,lty=2)
    plot(statsByPlane[[plan]][,2],type = "l",main = "",ylab = "n_a",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
    abline(v=numberPatternCooling,col="red",lwd=2,lty=2)
    dev.off()
    
    energy=energy+statsByPlane[[plan]][,1]+1000*statsByPlane[[plan]][,2]
  }
  energy=energy/length(statsByPlane)
  
  statisticsNumber=statsByPlane[[1]][,4]
  tableNumber=table(statisticsNumber)
  meanstat4=mean(statisticsNumber)
  medianstat4=median(statisticsNumber)
  modestat4=as.numeric(names(tableNumber)[which.max(tableNumber)])
  
  png(paste0("stats_",savename,"_number.png"),width = 1000,height = 300)
  par(mfrow = c(1, 2))
  par(mar = c(5, 6, 4, 1) + 0.1)
  # Left Plot: Line Plot with Annotations
  plot(statisticsNumber,
       type = "s",
       main = "",
       ylab = "n",
       col = palette.colors()[6],
       cex.lab = 2.5,
       cex.axis = 2)
  
  abline(h = meanstat4, col = palette.colors()[8], lwd = 2)
  abline(h = modestat4, col = palette.colors()[3], lwd = 2)
  abline(v = numberPatternCooling, col = "red", lwd = 2, lty = 2)
  abline(h = medianstat4, col = palette.colors()[2], lwd = 2)
  
  # Right Plot: Histogram or Barplot based on unique levels
  
  if (length(tableNumber) == 1) {
    hist(tableNumber / length(statisticsNumber),
         col = palette.colors()[6],
         xlab = "n",
         main = "",
         ylab = "Probability",
         cex.lab = 2,
         cex.axis = 2)
    
    legend("topright",
           legend = paste0(c("median = ", "mean = ", "mode = "),
                           round(c(medianstat4,
                                   meanstat4,
                                   modestat4), 2)),
           cex = 1.5)
    
  } else {
    hist(statisticsNumber,
         col = palette.colors()[6],
         prob = TRUE,
         nclass = 20,
         xlab = "n",
         main = "",
         cex.lab = 2,
         cex.axis = 2)
    
    abline(v = meanstat4, col = palette.colors()[8], lwd = 2)
    abline(v = modestat4, col = palette.colors()[3], lwd = 2)
    abline(v = medianstat4, col = palette.colors()[2], lwd = 2)
    
    legend("topright",
           legend = paste0(c("median = ", "mean = ", "mode = "),
                           round(c(medianstat4,
                                   meanstat4,
                                   modestat4), 2)),
           cex = 1.5)
  }
  
  
  dev.off()
  
  png(paste0("stats_",savename,"_strauss.png"),width = 1000,height = 200)
  par(mfrow=c(1,1))
  par(mar=c(5,6,4,1)+.1)
  plot(statsByPlane[[1]][,3],type = "l",main = "",ylab = "n_r",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
  abline(v=numberPatternCooling,col="red",lwd=2,lty=2)
  dev.off()
  
  
  
  png(paste0("energy_",savename,".png"),width = 1000,height = 400)
  par(mfrow=c(1,1))
  par(mar=c(5,6,4,1)+.1)
  plot(energy,type = "l",main = "",ylab = "Energy",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
  abline(v=numberPatternCooling,col="red",lwd=2,lty=2)
  dev.off()
  
  
  png(paste0("energyTemperature_",savename,".png"),width = 1000,height = 400)
  par(mfrow=c(1,1))
  par(mar=c(5,6,4,1)+.1)
  plot(energy/temperature,type = "l",main = "",ylab = "Energy / Temperature",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
  abline(v=numberPatternCooling,col="red",lwd=2,lty=2)
  dev.off()
  
  png(paste0("density_",savename,".png"),width = 1000,height = 400)
  par(mfrow=c(1,1))
  par(mar=c(5,6,4,1)+.1)
  plot(exp(-energy/temperature),type = "l",main = "",ylab = "density",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
  abline(v=numberPatternCooling,col="red",lwd=2,lty=2)
  dev.off()
  
  
  return(modestat4)
}

save_dendrogram=function(typedata,savename=typedata,numberSourcesDetected)
{
  sources=read.table(paste0("sources_final","_",typedata,".txt"),header=T)
  
  clustering=hclust(dist(sources),method = "ward.D2")
  clusters=cutree(clustering,k=20)
  
  # Find the height that corresponds to this cut (i.e., where to "zoom" in)
  cut_height <- sort(clustering$height, decreasing = TRUE)[20]
  
  png(paste0("stats_",savename,"_dendrogram.png"),width = 1000,height = 1000)
  plot(cut(as.dendrogram(clustering),h=cut_height)$upper,ylim=c(cut_height,max(clustering$height)),main="Dendrogram with a maximum of 20 clusters", leaflab = "none")
  dev.off()
  
  return(cutree(clustering,k=numberSourcesDetected))
}


save_detection=function(typedata,typedata2,clustersSources,nameSourcesDetected,n=50,annealing=T)
{
  if(annealing)
  {
    
    data=read.table(paste0("../DATA/",substring(typedata,1,nchar(typedata)-2),"_norm.txt"),header = T)
    dataSources=read.table(paste0("../DATA/",substring(typedata,1,nchar(typedata)-2),"_sources2d.txt"),header = T)
  }else{
    data=read.table(paste0("../DATA/",typedata,"_norm.txt"),header = T)
    dataSources=read.table(paste0("../DATA/",typedata,"_sources2d.txt"),header = T)
  }
  
  sources=read.table(paste0("sources_final","_",typedata2,".txt"),header=T)
  
  statsFinal=read.table(paste0("statistics_final","_",typedata2,".txt"),header=F)
  statsFinalByPlane=split(statsFinal,factor(paste0(statsFinal[,5],statsFinal[,6])))
  nSourcesPattern=statsFinalByPlane[[1]][,4]
  nbsim=nrow(statsFinal)/length(statsFinalByPlane)
  
  
  if((typedata=="mexico")|(typedata=="mexicoSA"))
  {
    dataTrue=read.table("../DATA/mexico.txt",header = T)
    
    res=reverse(dataTrue,data,sources,dataSources,positif=F)
    sources=res[["sources"]]
    dataSources=res[["dataSources"]]
    borne=res[["borne"]]
    data=dataTrue
  }else if((typedata=="athabasca")|(typedata=="athabascaSA"))
  {
    dataTrue=read.table("../DATA/athabasca_8d.txt",header = T)
    
    
    res=reverse(dataTrue,data,sources,dataSources,positif=T)
    sources=res[["sources"]]
    dataSources=res[["dataSources"]]
    borne=res[["borne"]]
    data=dataTrue
  }else{
    
    alldim=as.numeric(unique(c(levels(factor(dataSources[["dim1"]])),levels(factor(dataSources[["dim2"]])))))
    borne=data.frame(minCoord=rep(0,length(alldim)),maxCoord=rep(1,length(alldim)),dim=alldim)
  }
  
  
  
  saveplt(data,sources,dataSources,clustersSources,borne,nameSourcesDetected,nbsim,nSourcesPattern,n,savename=typedata)#n1,n2=> grid row,col
  
}

save_detection_article=function(typedata,typedata2,clustersSources,nameSourcesDetected,n=50,annealing=T)
{
  if(annealing)
  {
    
    data=read.table(paste0("../DATA/",substring(typedata,1,nchar(typedata)-2),"_norm.txt"),header = T)
    dataSources=read.table(paste0("../DATA/",substring(typedata,1,nchar(typedata)-2),"_sources2d.txt"),header = T)
  }else{
    data=read.table(paste0("../DATA/",typedata,"_norm.txt"),header = T)
    dataSources=read.table(paste0("../DATA/",typedata,"_sources2d.txt"),header = T)
  }
  
  sources=read.table(paste0("sources_final","_",typedata2,".txt"),header=T)
  
  statsFinal=read.table(paste0("statistics_final","_",typedata2,".txt"),header=F)
  statsFinalByPlane=split(statsFinal,factor(paste0(statsFinal[,5],statsFinal[,6])))
  nSourcesPattern=statsFinalByPlane[[1]][,4]
  nbsim=nrow(statsFinal)/length(statsFinalByPlane)
  
  
  if((typedata=="mexico")|(typedata=="mexicoSA"))
  {
    dataTrue=read.table("../DATA/mexico.txt",header = T)
    
    res=reverse(dataTrue,data,sources,dataSources,positif=F)
    sources=res[["sources"]]
    dataSources=res[["dataSources"]]
    borne=res[["borne"]]
    data=dataTrue
  }else if((typedata=="athabasca")|(typedata=="athabascaSA"))
  {
    dataTrue=read.table("../DATA/athabasca_8d.txt",header = T)
    
    
    res=reverse(dataTrue,data,sources,dataSources,positif=T)
    sources=res[["sources"]]
    dataSources=res[["dataSources"]]
    borne=res[["borne"]]
    data=dataTrue
  }else{
    
    alldim=as.numeric(unique(c(levels(factor(dataSources[["dim1"]])),levels(factor(dataSources[["dim2"]])))))
    borne=data.frame(minCoord=rep(0,length(alldim)),maxCoord=rep(1,length(alldim)),dim=alldim)
  }
  
  
  
  saveplt_article(data,sources,dataSources,clustersSources,borne,nameSourcesDetected,nbsim,nSourcesPattern,n,savename=typedata)#n1,n2=> grid row,col
  
}


saveplt_article=function(data,sources,real_sources,clustersSources,borne,nameSourcesDetected,nbsim,nSourcesPattern,n=50,savename=typedata)#n1,n2=> grid row,col
{
  
  # The sources detected are the mean of each cluster
  estimatedSources=aggregate(sources, by = list(cluster = clustersSources), FUN = mean)[,-1]
  write.table( round(estimatedSources,4), paste0(nameSourcesDetected,".txt"), row.names=F, col.names=T, append=F )
  estimatedSourcesStats=aggregate(sources, by = list(cluster = clustersSources), FUN = function(x) c(median=median(x),mean=mean(x),sd=sd(x)))
  write.table( round(estimatedSourcesStats,2), paste0(nameSourcesDetected,"_statistics.txt"), row.names=F, col.names=T, append=F )
  
  estimatedSourcesDist=rep(NA,length(levels(factor(clustersSources))))
  
  if(typedata=="syntheticSA")
  {
    real_s=read.table("../DATA/synthetic_sources3d.txt",header = T)
    estimatedSourcesDistanceToSources=data.frame(source=rep(NA,nrow(real_s)),dist=rep(NA,nrow(real_s)),
                                                 distDim1=rep(NA,nrow(real_s)),distDim2=rep(NA,nrow(real_s)),
                                                 distDim3=rep(NA,nrow(real_s)))
  }else{
    real_s=data.frame()
  }
  
  distanceCluster=data.frame(cluster=rep(NA,length(levels(factor(clustersSources)))),dist=rep(NA,length(levels(factor(clustersSources)))),
                             distDim1=rep(NA,length(levels(factor(clustersSources)))),distDim2=rep(NA,length(levels(factor(clustersSources)))),
                             distDim3=rep(NA,length(levels(factor(clustersSources)))))
  
  for(i in 1:length(levels(factor(clustersSources))))
  {
    estimatedSourcesDist[i]=max(dist(sources[clustersSources==i,]))
    clus=sources[clustersSources==i,]
    distanceCluster[i,1]=i
    distanceCluster[i,2]=max(dist(sources[clustersSources==i,]))
    distanceCluster[i,3]=max(clus[,1])-min(clus[,1])
    distanceCluster[i,4]=max(clus[,2])-min(clus[,2])
    distanceCluster[i,5]=max(clus[,3])-min(clus[,3])
    if(i<=nrow(real_s))
    {
      estimatedSourcesDistanceToSources[i,1]=which.min(as.matrix(dist(rbind(real_s[i,],estimatedSources)))[1,-1])
      estimatedSourcesDistanceToSources[i,2]=min(as.matrix(dist(rbind(real_s[i,],estimatedSources)))[1,-1])
      estimatedSourcesDistanceToSources[i,3]=abs(real_s[i,1]-estimatedSources[which.min(as.matrix(dist(rbind(real_s[i,],estimatedSources)))[1,-1]),1])
      estimatedSourcesDistanceToSources[i,4]=abs(real_s[i,2]-estimatedSources[which.min(as.matrix(dist(rbind(real_s[i,],estimatedSources)))[1,-1]),2])
      estimatedSourcesDistanceToSources[i,5]=abs(real_s[i,3]-estimatedSources[which.min(as.matrix(dist(rbind(real_s[i,],estimatedSources)))[1,-1]),3])
      
    }
    
    
  }
  write.table( round(estimatedSourcesDist,4), paste0(nameSourcesDetected,"_distMax.txt"), row.names=F, col.names=T, append=F )
  write.table( round(distanceCluster,4), paste0(nameSourcesDetected,"_distInCluster.txt"), row.names=F, col.names=T, append=F )
  
  print(xtable(round(distanceCluster,4), type = "latex"), file = paste0(nameSourcesDetected,"_distInCluster.tex"))
  
  if(typedata=="syntheticSA")
  {
    write.table( round(estimatedSourcesDistanceToSources,4), paste0(nameSourcesDetected,"_distToSources.txt"), row.names=F, col.names=T, append=F )
    print(xtable(round(estimatedSourcesDistanceToSources,4), type = "latex"), file = paste0(nameSourcesDetected,"_distToSources.tex"))
    
  }
  print("Detected sources saved!")
  
  if(typedata=="athabascaSA")
  {
    sourceRic=read.table("../DATA/athabasca_sourcesRic16.txt",header = T)
    
  }else{
    sourceRic=NULL
    
  }
  nbcluster=nrow(estimatedSources)
  medianpattern=NULL
  
  plan=0
  for (abscisse in 1:(ncol(sources)-1))
  {
    for(ordonnee in (abscisse+1):ncol(sources))
    {
      plan=plan+1
      maxord=borne[borne[["dim"]]==ordonnee,2]
      minord=borne[borne[["dim"]]==ordonnee,1]
      minabs=borne[borne[["dim"]]==abscisse,1]
      maxabs=borne[borne[["dim"]]==abscisse,2]
      
      real_sources_plane=real_sources[(real_sources$dim1==abscisse)&(real_sources$dim2==ordonnee),]
      
      
      maxord=borne[borne[["dim"]]==ordonnee,2]
      minord=borne[borne[["dim"]]==ordonnee,1]
      minabs=borne[borne[["dim"]]==abscisse,1]
      maxabs=borne[borne[["dim"]]==abscisse,2]
      #
      # # Define the grid
      x_breaks <- seq(minabs, maxabs, length.out = n)  # X-axis grid breaks
      y_breaks <- seq(minord, maxord, length.out = n)  # Y-axis grid breaks
      #
      # # Assign points to grid cells
      x_bins_iteration <- cut(sources[1:nSourcesPattern[1], abscisse], breaks = x_breaks, include.lowest = TRUE)
      y_bins_iteration <- cut(sources[1:nSourcesPattern[1], ordonnee], breaks = y_breaks, include.lowest = TRUE)
      notDuplicated = (!duplicated(x_bins_iteration)) & (!duplicated(y_bins_iteration))
      
      x_bins <- x_bins_iteration[notDuplicated]
      y_bins <- y_bins_iteration[notDuplicated]
      
      for(indicePattern in 2:(nbsim))
      {
        sourcesPatternIteration=sources[(sum(nSourcesPattern[1:(indicePattern-1)])+1):(sum(nSourcesPattern[1:(indicePattern-1)])+nSourcesPattern[indicePattern]),]
        x_bins_iteration = cut(sourcesPatternIteration[, abscisse], breaks = x_breaks, include.lowest = TRUE)
        y_bins_iteration = cut(sourcesPatternIteration[, ordonnee], breaks = y_breaks, include.lowest = TRUE)
        notDuplicated = (!duplicated(x_bins_iteration)) & (!duplicated(y_bins_iteration))
        x_bins =c(x_bins, x_bins_iteration[notDuplicated])
        y_bins =c(y_bins, y_bins_iteration[notDuplicated])
      }
      #
      # # Create a 2D table of counts
      grid_counts <- table(x_bins, y_bins)
      #
      # # Convert table to matrix
      z <- as.matrix(grid_counts)
      #
      # # Transpose to align with x/y axes
      z <- z/nbsim
      #
      png(file = paste0("lvl_",savename,"_plane_",plan,".png"), width = 800, height = 600)
      par(mar = c(8, 8, 4, 1) + 0.1,mgp = c(5, 1.7, 0) )
      
      # Plot using `image()` for blocky coloring
      if(typedata=="mexicoSA")
      {
        if(plan==1)
        {
          image.plot(
            x = x_breaks, y = y_breaks, z = z,  # Transpose z for proper orientation
            col = colfunc,                     # Use fewer distinct colors
            cex.axis=3,cex.lab=3,
            xlab = colnames(sources)[[abscisse]],
            ylab = expression(paste(delta^37,"Cl"))
          )
        }else if(plan==2)
        {
          image.plot(
            x = x_breaks, y = y_breaks, z = z,  # Transpose z for proper orientation
            col = colfunc,                     # Use fewer distinct colors
            cex.axis=3,cex.lab=3,
            xlab = colnames(sources)[[abscisse]],
            ylab = expression(paste(delta^81,"Br"))
          )
        }else{
          image.plot(
            x = x_breaks, y = y_breaks, z = z,  # Transpose z for proper orientation
            col = colfunc,                     # Use fewer distinct colors
            cex.axis=3,cex.lab=3,
            xlab = expression(paste(delta^37,"Cl")),
            ylab = expression(paste(delta^81,"Br"))
          )
        }
      }else{
        image.plot(
          x = x_breaks, y = y_breaks, z = z,  # Transpose z for proper orientation
          col = colfunc,                     # Use fewer distinct colors
          cex.axis=3,cex.lab=3,
          xlab = colnames(sources)[[abscisse]],
          ylab = colnames(sources)[[ordonnee]]# expression(paste(delta^81,"Br"))
        )
      }
      
      
      # Add grid lines
      abline(h = seq(minord, maxord, length.out = n), col = "lightgray",lwd=0.1)
      abline(v = seq(minabs, maxabs, length.out = n), col = "lightgray",lwd=0.1)
      points(data[, c(abscisse, ordonnee)], pch = 19, col = "darkgray", cex = 1)
      # points(sourceRic[c(1,2),c(abscisse, ordonnee)],pch=17,cex=3)
      # lines(sourceRic[c(1:2,1),c(abscisse, ordonnee)],lwd=10,col="blue")
      lines(x=c(sourceRic[1,c(abscisse)],sourceRic[2,c( abscisse)]),y=c(sourceRic[1,c(ordonnee)],sourceRic[1,c( ordonnee)]),lwd=5,col="blue")
      lines(x=c(sourceRic[1,c(abscisse)],sourceRic[2,c( abscisse)]),y=c(sourceRic[2,c(ordonnee)],sourceRic[2,c( ordonnee)]),lwd=5,col="blue")
      lines(x=c(sourceRic[1,c(abscisse)],sourceRic[1,c( abscisse)]),y=c(sourceRic[1,c(ordonnee)],sourceRic[2,c( ordonnee)]),lwd=5,col="blue")
      lines(x=c(sourceRic[2,c(abscisse)],sourceRic[2,c( abscisse)]),y=c(sourceRic[1,c(ordonnee)],sourceRic[2,c( ordonnee)]),lwd=5,col="blue")
      
      lines(x=c(sourceRic[3,c(abscisse)],sourceRic[4,c( abscisse)]),y=c(sourceRic[3,c(ordonnee)],sourceRic[3,c( ordonnee)]),lwd=5,col="red")
      lines(x=c(sourceRic[3,c(abscisse)],sourceRic[4,c( abscisse)]),y=c(sourceRic[4,c(ordonnee)],sourceRic[4,c( ordonnee)]),lwd=5,col="red")
      lines(x=c(sourceRic[3,c(abscisse)],sourceRic[3,c( abscisse)]),y=c(sourceRic[3,c(ordonnee)],sourceRic[4,c( ordonnee)]),lwd=5,col="red")
      lines(x=c(sourceRic[4,c(abscisse)],sourceRic[4,c( abscisse)]),y=c(sourceRic[3,c(ordonnee)],sourceRic[4,c( ordonnee)]),lwd=5,col="red")
      
      
      # Add points
      points(real_sources_plane[, c(1, 2)], pch = 17, col = "red", cex = 3)
      
      points(estimatedSources[, c(abscisse, ordonnee)], pch = 19, col = 'darkgreen', cex = 3)
      
      if(typedata=="syntheticSA")
      {
        text(
          estimatedSources[1, abscisse]+0.02,
          estimatedSources[1, ordonnee],
          labels = 1,
          pos = 4,
          col = "darkgreen",
          cex = 3
        )
        text(
          estimatedSources[2, abscisse],
          estimatedSources[2, ordonnee]-0.02,
          labels = 2,
          pos = 1,
          col = "darkgreen",
          cex = 3
        )
        if(plan==1)
        {
          text(
            estimatedSources[3, abscisse],
            estimatedSources[3, ordonnee]+0.02,
            labels = 3,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else{
          text(
            estimatedSources[3, abscisse],
            estimatedSources[3, ordonnee]-0.02,
            labels = 3,
            pos = 1,
            col = "darkgreen",
            cex = 3
          )
        }
        
        text(
          estimatedSources[4, abscisse],
          estimatedSources[4, ordonnee]+0.02,
          labels = 4,
          pos = 3,
          col = "darkgreen",
          cex = 3
        )
      }
      
      if(typedata=="mexicoSA")
      {
        text(
          estimatedSources[1, abscisse]+0.2,
          estimatedSources[1, ordonnee],
          labels = 1,
          pos = 4,
          col = "darkgreen",
          cex = 3
        )
        text(
          estimatedSources[2, abscisse]+0.2,
          estimatedSources[2, ordonnee],
          labels = 2,
          pos = 4,
          col = "darkgreen",
          cex = 3
        )
        text(
          estimatedSources[3, abscisse]-0.2,
          estimatedSources[3, ordonnee],
          labels = 3,
          pos = 2,
          col = "darkgreen",
          cex = 3
        )
        

      }
      
      if(typedata=="athabascaSA")
      {
        stepabs=(maxabs-minabs)/n
        stepord=(maxord-minord)/n
        
        if(plan<5)
        {
          text(
            estimatedSources[1, abscisse],
            estimatedSources[1, ordonnee]+stepord,
            labels = 1,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else{
          text(
            estimatedSources[1, abscisse]+stepabs,
            estimatedSources[1, ordonnee]+stepord,
            labels = 1,
            pos = 4,
            col = "darkgreen",
            cex = 3
          )
        }
        if(plan==1)
        {
          text(
            estimatedSources[2, abscisse]+stepabs,
            estimatedSources[2, ordonnee],
            labels = 2,
            pos = 4,
            col = "darkgreen",
            cex = 3
          )
        }else if((plan==2)|(plan==3)|(plan==4))
        {
          text(
            estimatedSources[2, abscisse]+stepabs,
            estimatedSources[2, ordonnee]+stepord,
            labels = 2,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else if((plan==5)|(plan==6)|(plan==7))
        {
          text(
            estimatedSources[2, abscisse],
            estimatedSources[2, ordonnee]+stepord,
            labels = 2,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else if((plan==8)|(plan==10))
        {
          text(
            estimatedSources[2, abscisse]+stepabs,
            estimatedSources[2, ordonnee]+2*stepord,
            labels = 2,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else if((plan==9))
        {
          text(
            estimatedSources[2, abscisse]+3*stepabs,
            estimatedSources[2, ordonnee]+3*stepord,
            labels = 2,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else if((plan==10))
        {
          text(
            estimatedSources[2, abscisse]+3*stepabs,
            estimatedSources[2, ordonnee]+3*stepord,
            labels = 2,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }
        
        
        if((plan==1)){
          text(
            estimatedSources[3, abscisse]+stepabs,
            estimatedSources[3, ordonnee]+stepord,
            labels = 3,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else if((plan==4)|(plan==7)|(plan==9)|(plan==10))
        {
          text(
            estimatedSources[3, abscisse]+stepabs,
            estimatedSources[3, ordonnee],
            labels = 3,
            pos = 4,
            col = "darkgreen",
            cex = 3
          )
        }else if((plan==5)|(plan==6))
        {
          text(
            estimatedSources[3, abscisse]+stepabs,
            estimatedSources[3, ordonnee]+stepord,
            labels = 3,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else if ((plan==2)|(plan==3)|(plan==8)){
          text(
            estimatedSources[3, abscisse]+3*stepabs,
            estimatedSources[3, ordonnee]+3*stepord,
            labels = 3,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }
        
        if((plan==3)|(plan==6)|(plan==8))
        {
          text(
            estimatedSources[4, abscisse]+stepabs,
            estimatedSources[4, ordonnee],
            labels = 4,
            pos = 4,
            col = "darkgreen",
            cex = 3
          )
        }else if ((plan==1)|(plan==2)|(plan==4)|(plan==7))
        {
          text(
            estimatedSources[4, abscisse]+stepabs,
            estimatedSources[4, ordonnee]+stepord,
            labels = 4,
            pos = 4,
            col = "darkgreen",
            cex = 3
          )
        }else if ((plan==9)|(plan==10))
        {
          text(
            estimatedSources[4, abscisse],
            estimatedSources[4, ordonnee]+stepord,
            labels = 4,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else if ((plan==5))
        {
          text(
            estimatedSources[4, abscisse]+1*stepabs,
            estimatedSources[4, ordonnee]+3*stepord,
            labels = 4,
            pos = 4,
            col = "darkgreen",
            cex = 3
          )
        }else if ((plan==7))
        {
          text(
            estimatedSources[4, abscisse]+stepabs,
            estimatedSources[4, ordonnee]+2*stepord,
            labels = 4,
            pos = 4,
            col = "darkgreen",
            cex = 3
          )
        }
        
        if((plan==2)|(plan==5))
        {
          text(
            estimatedSources[5, abscisse]+stepabs,
            estimatedSources[5, ordonnee],
            labels = 5,
            pos = 4,
            col = "darkgreen",
            cex = 3
          )
        }else if ((plan==1)|(plan==4)|(plan==6)|(plan==10))
        {
          text(
            estimatedSources[5, abscisse]+3*stepabs,
            estimatedSources[5, ordonnee]+3*stepord,
            labels = 5,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else if ((plan==9)|(plan==8))
        {
          text(
            estimatedSources[5, abscisse],
            estimatedSources[5, ordonnee]+stepord,
            labels = 5,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }else if ((plan==3))
        {
          text(
            estimatedSources[5, abscisse]+3*stepabs,
            estimatedSources[5, ordonnee]+1*stepord,
            labels = 5,
            pos = 4,
            col = "darkgreen",
            cex = 3
          )
        }else if ((plan==7))
        {
          text(
            estimatedSources[5, abscisse]+stepabs,
            estimatedSources[5, ordonnee]+stepord,
            labels = 5,
            pos = 3,
            col = "darkgreen",
            cex = 3
          )
        }

      }
      # Add labels to points
      

      
      
      # Close the PNG device
      dev.off()
      
      
      
      
      print(paste("Levels sets saved for the plane number: ",plan,sep=""))
    }
  }
  
  
}

saveplt=function(data,sources,real_sources,clustersSources,borne,nameSourcesDetected,nbsim,nSourcesPattern,n=50,savename=typedata)#n1,n2=> grid row,col
{
  
  # The sources detected are the mean of each cluster
  estimatedSources=aggregate(sources, by = list(cluster = clustersSources), FUN = mean)[,-1]
  write.table( round(estimatedSources,4), paste0(nameSourcesDetected,".txt"), row.names=F, col.names=T, append=F )
  estimatedSourcesStats=aggregate(sources, by = list(cluster = clustersSources), FUN = function(x) c(median=median(x),mean=mean(x),sd=sd(x)))
  write.table( round(estimatedSourcesStats,2), paste0(nameSourcesDetected,"_statistics.txt"), row.names=F, col.names=T, append=F )
  
  estimatedSourcesDist=rep(NA,length(levels(factor(clustersSources))))
  
  
  distanceCluster=data.frame(cluster=rep(NA,length(levels(factor(clustersSources)))),dist=rep(NA,length(levels(factor(clustersSources)))),
                             distDim1=rep(NA,length(levels(factor(clustersSources)))),distDim2=rep(NA,length(levels(factor(clustersSources)))),
                             distDim3=rep(NA,length(levels(factor(clustersSources)))))
  
  for(i in 1:length(levels(factor(clustersSources))))
  {
    estimatedSourcesDist[i]=max(dist(sources[clustersSources==i,]))
    clus=sources[clustersSources==i,]
    distanceCluster[i,1]=i
    distanceCluster[i,2]=max(dist(sources[clustersSources==i,]))
    distanceCluster[i,3]=max(clus[,1])-min(clus[,1])
    distanceCluster[i,4]=max(clus[,2])-min(clus[,2])
    distanceCluster[i,5]=max(clus[,3])-min(clus[,3])
    
    
    
  }
  write.table( round(estimatedSourcesDist,4), paste0(nameSourcesDetected,"_distMax.txt"), row.names=F, col.names=T, append=F )
  write.table( round(distanceCluster,4), paste0(nameSourcesDetected,"_distInCluster.txt"), row.names=F, col.names=T, append=F )
  
  print(xtable(round(distanceCluster,4), type = "latex"), file = paste0(nameSourcesDetected,"_distInCluster.tex"))
  print("Detected sources saved!")
  
  
  nbcluster=nrow(estimatedSources)
  medianpattern=NULL
  
  plan=0
  for (abscisse in 1:(ncol(sources)-1))
  {
    for(ordonnee in (abscisse+1):ncol(sources))
    {
      plan=plan+1
      maxord=borne[borne[["dim"]]==ordonnee,2]
      minord=borne[borne[["dim"]]==ordonnee,1]
      minabs=borne[borne[["dim"]]==abscisse,1]
      maxabs=borne[borne[["dim"]]==abscisse,2]
      
      real_sources_plane=real_sources[(real_sources$dim1==abscisse)&(real_sources$dim2==ordonnee),]
      
      
      maxord=borne[borne[["dim"]]==ordonnee,2]
      minord=borne[borne[["dim"]]==ordonnee,1]
      minabs=borne[borne[["dim"]]==abscisse,1]
      maxabs=borne[borne[["dim"]]==abscisse,2]
      #
      # # Define the grid
      x_breaks <- seq(minabs, maxabs, length.out = n)  # X-axis grid breaks
      y_breaks <- seq(minord, maxord, length.out = n)  # Y-axis grid breaks
      #
      # # Assign points to grid cells
      x_bins_iteration <- cut(sources[1:nSourcesPattern[1], abscisse], breaks = x_breaks, include.lowest = TRUE)
      y_bins_iteration <- cut(sources[1:nSourcesPattern[1], ordonnee], breaks = y_breaks, include.lowest = TRUE)
      notDuplicated = (!duplicated(x_bins_iteration)) & (!duplicated(y_bins_iteration))
      
      x_bins <- x_bins_iteration[notDuplicated]
      y_bins <- y_bins_iteration[notDuplicated]
      
      for(indicePattern in 2:(nbsim))
      {
        sourcesPatternIteration=sources[(sum(nSourcesPattern[1:(indicePattern-1)])+1):(sum(nSourcesPattern[1:(indicePattern-1)])+nSourcesPattern[indicePattern]),]
        x_bins_iteration = cut(sourcesPatternIteration[, abscisse], breaks = x_breaks, include.lowest = TRUE)
        y_bins_iteration = cut(sourcesPatternIteration[, ordonnee], breaks = y_breaks, include.lowest = TRUE)
        notDuplicated = (!duplicated(x_bins_iteration)) & (!duplicated(y_bins_iteration))
        x_bins =c(x_bins, x_bins_iteration[notDuplicated])
        y_bins =c(y_bins, y_bins_iteration[notDuplicated])
      }
      
      
      #
      # # Create a 2D table of counts
      grid_counts <- table(x_bins, y_bins)
      #
      # # Convert table to matrix
      z <- as.matrix(grid_counts)
      #
      # # Transpose to align with x/y axes
      z <- z/nbsim
      #
      png(file = paste0("lvl_",savename,"_plane_",plan,".png"), width = 800, height = 600)
      par(mar = c(8, 8, 4, 1) + 0.1,mgp = c(5, 1.7, 0) )
      
      # Plot using `image()` for blocky coloring
      image.plot(
        x = x_breaks, y = y_breaks, z = z,  # Transpose z for proper orientation
        col = colfunc,                     # Use fewer distinct colors
        cex.axis=3,cex.lab=3,
        xlab = colnames(sources)[[abscisse]],
        ylab = colnames(sources)[[ordonnee]]
      )
      
      # Add grid lines
      abline(h = seq(minord, maxord, length.out = n), col = "lightgray",lwd=0.1)
      abline(v = seq(minabs, maxabs, length.out = n), col = "lightgray",lwd=0.1)
      points(data[, c(abscisse, ordonnee)], pch = 19, col = "darkgray", cex = 1)
      
      # Add points
      points(real_sources_plane[, c(1, 2)], pch = 17, col = "red", cex = 3)
      
      points(estimatedSources[, c(abscisse, ordonnee)], pch = 19, col = 'darkgreen', cex = 3)
      
      # Add labels to points
      text(
        estimatedSources[, abscisse],
        estimatedSources[, ordonnee],
        labels = 1:nrow(estimatedSources),
        pos = 4,
        col = "darkgreen",
        cex = 3
      )
      
      
      
      # Close the PNG device
      dev.off()
      
      
      
      
      print(paste("Levels sets saved for the plane number: ",plan,sep=""))
    }
  }
  
  
}



reverse=function(dataTrue,data,sources,dataSources,positif=F)
{
  
  alldim=as.numeric(unique(c(levels(factor(dataSources[["dim1"]])),levels(factor(dataSources[["dim2"]])))))
  borne=data.frame(minCoord=rep(0,length(alldim)),maxCoord=rep(1,length(alldim)),dim=alldim)
  for(i in 1:ncol(sources))
  {
    etendu=max(dataTrue[,i],na.rm=T)-min(dataTrue[,i],na.rm=T)
    if(positif)
    {
      binf=max(min(dataTrue[,i],na.rm=T)-etendu,0)
      minWindow=0
    }else{
      binf=min(dataTrue[,i],na.rm=T)-etendu
      minWindow=binf
    }
    bsup=max(dataTrue[,i],na.rm=T)+etendu
    minimu=min(data[,i],na.rm=T)
    sources[,i]=sources[,i]*(bsup-binf)+minWindow
    dataSources[dataSources$dim1==i,1]=dataSources[dataSources$dim1==i,1]*(bsup-binf)+minWindow
    dataSources[dataSources$dim2==i,2]=dataSources[dataSources$dim2==i,2]*(bsup-binf)+minWindow
    borne[borne[["dim"]]==i,c(1,2)]=borne[borne[["dim"]]==i,c(1,2)]*(bsup-binf)+minWindow
    
    
  }
  
  return(list(sources=sources,dataSources=dataSources,borne=borne))
  
}



########################################################################################
########################################################################################
########################################################################################

#                     Plot results article

########################################################################################
########################################################################################
########################################################################################


setwd("C++/RESULTS")


typedata=paste0("syntheticSA")

numberSourcesDetected=save_plot_statistics(typedata,typedata)
clustersSources=save_dendrogram(typedata,typedata,numberSourcesDetected)
nameSourcesDetected=paste0("sourcesDetected_",typedata)
save_detection_article("syntheticSA",typedata,clustersSources,paste0("sourcesDetected_",typedata),n=50,annealing=T)


typedata="mexicoSA"

numberSourcesDetected=save_plot_statistics(typedata,typedata)
clustersSources=save_dendrogram(typedata,typedata,numberSourcesDetected)
nameSourcesDetected=paste0("sourcesDetected_",typedata)
save_detection_article(typedata,typedata,clustersSources,paste0("sourcesDetected_",typedata),n=50,annealing=T)

typedata="athabascaSA"

numberSourcesDetected=save_plot_statistics(typedata,typedata)
clustersSources=save_dendrogram(typedata,typedata,numberSourcesDetected)
nameSourcesDetected=paste0("sourcesDetected_",typedata)
save_detection_article(typedata,typedata,clustersSources,paste0("sourcesDetected_",typedata),n=50,annealing=T)

########################################################################################
########################################################################################
########################################################################################

#                     Plot results

########################################################################################
########################################################################################
########################################################################################


setwd("C++/RESULTS")


for(i in 1:35)
{
  typedata=paste0("syntheticSA",i)
  typedata=paste0("syntheticSA")
  
  numberSourcesDetected=save_plot_statistics(typedata,typedata)
  clustersSources=save_dendrogram(typedata,typedata,numberSourcesDetected)
  nameSourcesDetected=paste0("sourcesDetected_",typedata)
  save_detection("syntheticSA",typedata,clustersSources,paste0("sourcesDetected_",typedata),n=50,annealing=T)
}



typedata="mexicoSA"

numberSourcesDetected=save_plot_statistics(typedata,typedata)
clustersSources=save_dendrogram(typedata,typedata,numberSourcesDetected)
nameSourcesDetected=paste0("sourcesDetected_",typedata)
save_detection(typedata,typedata,clustersSources,paste0("sourcesDetected_",typedata),n=50,annealing=T)

typedata="athabascaSA"

numberSourcesDetected=save_plot_statistics(typedata,typedata)
clustersSources=save_dendrogram(typedata,typedata,numberSourcesDetected)
nameSourcesDetected=paste0("sourcesDetected_",typedata)
save_detection(typedata,typedata,clustersSources,paste0("sourcesDetected_",typedata),n=50,annealing=T)


########################################################################################
########################################################################################
########################################################################################

#                     Useful plots

########################################################################################
########################################################################################
########################################################################################

data=read.table("C++/DATA/synthetic_norm.txt",header = T)
real_sources=read.table("C++/DATA/synthetic_sources3d.txt",header = T)
real_sources2d=read.table("C++/DATA/synthetic_sources2d.txt",header = T)

colnames(real_sources)=c("solute 1", "solute 2", "solute 3")
colnames(data)=colnames(real_sources)
real_sources=real_sources[c(3,2,4,1),]

#Plot of the synthetic data set
png(file="plot3dsynthetic.png",width = 800,height = 600)
pl=scatterplot3d(rbind(real_sources,data),cex.lab=2.5,pch=c(rep(9,nrow(real_sources)),
                                                            rep(20,nrow(data))),color=c(rep(4,nrow(real_sources)),rep(1,nrow(data))),
                 tick.marks=F,cex.symbols = c(rep(3,nrow(real_sources)),rep(1.5,nrow(data))),
                 xlim=c(0,1),ylim=c(0,1),zlim=c(0,1),mar = c(7, 6, 5, 5),xpd=NA)
lgd=cbind(real_sources[,1]+0.05,real_sources[,2]+0.05,real_sources[,3]-0.1)
text(pl$xyz.convert(lgd), labels = 1:4,
     cex= 4, col = 4)
dev.off()


#Plot the projection on the first plane of the synthetic data set
png(file="synthetic1.png",width = 800,height = 600)
par(mar=c(5,6,4,1)+.1)
plot(real_sources[,c(1,2)],xlim=c(0,1),ylim=c(0,1),pch=9,col=4,cex=2.5,cex.lab=2.5,cex.axis=2.5,xaxt="n",yaxt="n")
points(data[,c(1,2)],pch=20)
text(
  real_sources[3:4, 1],
  real_sources[3:4, 2]+0.02,
  labels = 3:4,  # Number each red point
  pos = 3,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
text(
  real_sources[1:2, 1]+0.02,
  real_sources[1:2, 2],
  labels = 1:2,  # Number each red point
  pos = 4,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
lines(unique(real_sources[,c(1,2)])[c(1:(nrow(unique(real_sources[,c(1,2)]))),1),],lty=2,lwd=3,col=2)
dev.off()

#Plot the projection on the second plane of the synthetic data set
png(file="synthetic2.png",width = 800,height = 600)
par(mar=c(5,6,4,1)+.1)
plot(real_sources[,c(1,3)],xlim=c(0,1),ylim=c(0,1),pch=9,col=4,cex=2.5,cex.lab=2.5,cex.axis=2.5,xaxt="n",yaxt="n")
points(data[,c(1,3)],pch=20)
text(
  real_sources[3:4, 1],
  real_sources[3:4, 3]+0.02,
  labels = 3:4,  # Number each red point
  pos = 3,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
text(
  real_sources[1, 1]+0.02,
  real_sources[1, 3],
  labels = 1,  # Number each red point
  pos = 4,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
text(
  real_sources[2, 1],
  real_sources[2, 3]-0.02,
  labels = 2,  # Number each red point
  pos = 1,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
lines(unique(real_sources[,c(1,3)])[c(1:(nrow(unique(real_sources[,c(1,3)]))),1),],lty=2,lwd=3,col=2)
dev.off()

#Plot the projection on the third plane of the synthetic data set
png(file="synthetic3.png",width = 800,height = 600)
par(mar=c(5,6,4,1)+.1)
plot(real_sources[,c(2,3)],xlim=c(0,1),ylim=c(0,1),pch=9,col=4,cex=2.5,cex.lab=2.5,cex.axis=2.5,xaxt="n",yaxt="n")
points(data[,c(2,3)],pch=20)
text(
  real_sources[2:3, 2],
  real_sources[2:3, 3]+0.02,
  labels = 2:3,  # Number each red point
  pos = 3,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
text(
  real_sources[4, 2]-0.02,
  real_sources[4, 3],
  labels = 4,  # Number each red point
  pos = 2,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
text(
  real_sources[1, 2]+0.02,
  real_sources[1, 3],
  labels = 1,  # Number each red point
  pos = 4,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
lines(unique(real_sources[,c(2,3)])[c(1:(nrow(unique(real_sources[,c(2,3)]))),1),],lty=2,lwd=3,col=2)
dev.off()


#Plot the region of interest of the first plane of the synthetic data set
png(file="synthetic1bis.png",width = 800,height = 600)
par(mar=c(5,6,4,1)+.1)
plot(real_sources[,c(1,2)],xlim=c(0,1),ylim=c(0,1),pch = 17, col = "red",cex=3,cex.lab=2.5,cex.axis=2.5)
points(data[,c(1,2)],pch=20)
for(i in 2:4)
{
  symbols(real_sources[i,c(1,2)], circles =0.1376, add = TRUE, fg = "red",bg = rgb(1, 0, 0, alpha = 0.1),inches = FALSE)
}
dev.off()

#Plot the region of interest of the second plane of the synthetic data set
png(file="synthetic2bis.png",width = 800,height = 600)
par(mar=c(5,6,4,1)+.1)
plot(real_sources[,c(1,3)],xlim=c(0,1),ylim=c(0,1),pch = 17, col = "red",cex=3,cex.lab=2.5,cex.axis=2.5)
points(data[,c(1,3)],pch=20)
for(i in 2:4)
{
  symbols(real_sources[i,c(1,3)], circles =0.1357, add = TRUE, fg = "red",bg = rgb(1, 0, 0, alpha = 0.1),inches = FALSE)
}
dev.off()

#Plot the region of interest of the third plane of the synthetic data set
png(file="synthetic3bis.png",width = 800,height = 600)
par(mar=c(5,6,4,1)+.1)
plot(real_sources[,c(2,3)],xlim=c(0,1),ylim=c(0,1),pch = 17, col = "red",cex=3,cex.lab=2.5,cex.axis=2.5)
points(data[,c(2,3)],pch=20)
for(i in 1:3)
{
  symbols(real_sources[i,c(2,3)], circles =0.137, add = TRUE, fg = "red",bg = rgb(1, 0, 0, alpha = 0.1),inches = FALSE)
}
dev.off()
