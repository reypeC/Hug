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

setwd("C:/Users/creype/Documents/Hug/simHug/Hug250617")



save_plot_statistics=function(typedata,savename=typedata,annealing=T)
{
  stats=read.table(paste0("statistics","_",typedata,".txt"),header=T)
  statsFinal=read.table(paste0("statistics_final","_",typedata,".txt"),header=F)
  colnames(statsFinal)=colnames(stats)
  
  statsComplete=rbind(stats,statsFinal)
  
  statsByPlane=split(statsComplete,factor(paste0(statsComplete[,5],statsComplete[,6])))
  
  numberPatternCooling=nrow(stats)/length(statsByPlane)
  
  for (plan in 1:length(statsByPlane)) {
    
    
    png(paste0("stats_",savename,"_plane_",plan,".png"),width = 1000,height = 200)
    par(mfrow=c(1,2))
    par(mar=c(5,6,4,1)+.1)
    plot(statsByPlane[[plan]][,1],type = "l",main = "",ylab = "n_e",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
    abline(v=numberPatternCooling,col="red",lwd=2,lty=2)
    plot(statsByPlane[[plan]][,2],type = "l",main = "",ylab = "n_a",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
    abline(v=numberPatternCooling,col="red",lwd=2,lty=2)
    dev.off()
    
  }
  
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
       type = "l",
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

save_detection=function(typedata,clustersSources,nameSourcesDetected,n=50,annealing=T)
{
  if(annealing)
  {
    
    data=read.table(paste0("../DATA/",substring(typedata,1,nchar(typedata)-2),"_norm.txt"),header = T)
    dataSources=read.table(paste0("../DATA/",substring(typedata,1,nchar(typedata)-2),"_sources2d.txt"),header = T)
  }else{
    data=read.table(paste0("../DATA/",typedata,"_norm.txt"),header = T)
    dataSources=read.table(paste0("../DATA/",typedata,"_sources2d.txt"),header = T)
  }
  
  sources=read.table(paste0("sources_final","_",typedata,".txt"),header=T)
  
  statsFinal=read.table(paste0("statistics_final","_",typedata,".txt"),header=F)
  nbsim=nrow(statsFinal)/length(split(statsFinal,factor(paste0(statsFinal[,5],statsFinal[,6]))))
  
  
  if((typedata=="pinti")|(typedata=="pintiSA"))
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
  

  
  saveplt(data,sources,dataSources,clustersSources,borne,nameSourcesDetected,nbsim,n)#n1,n2=> grid row,col
 
}



saveplt=function(data,sources,real_sources,clustersSources,borne,nameSourcesDetected,nbsim,n=50,savename=typedata)#n1,n2=> grid row,col
{
  
  # The sources detected are the mean of each cluster
  estimatedSources=aggregate(sources, by = list(cluster = clustersSources), FUN = mean)[,-1]
  write.table( estimatedSources, paste0(nameSourcesDetected,".txt"), row.names=F, col.names=T, append=F )
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
      x_bins <- cut(sources[, abscisse], breaks = x_breaks, include.lowest = TRUE)
      y_bins <- cut(sources[, ordonnee], breaks = y_breaks, include.lowest = TRUE)
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
      par(mar = c(5, 6, 4, 1) + 0.1)
      
      # Plot using `image()` for blocky coloring
      image.plot(
        x = x_breaks, y = y_breaks, z = z,  # Transpose z for proper orientation
        col = colfunc,                     # Use fewer distinct colors
        cex.axis=3,cex.lab=3,
        xlab = colnames(sources)[[abscisse]],
        ylab = colnames(sources)[[ordonnee]]
      )

      # Add grid lines
      abline(h = seq(minord, maxord, length.out = n), col = "gray")
      abline(v = seq(minabs, maxabs, length.out = n), col = "gray")
      
      
      # Add points
      points(estimatedSources[, c(abscisse, ordonnee)], pch = 19, col = 'darkgreen', cex = 3)
      points(real_sources_plane[, c(1, 2)], pch = 17, col = "red", cex = 3)
      points(data[, c(abscisse, ordonnee)], pch = 19, col = "black", cex = 1)
      
      # Add labels to points
      text(
        estimatedSources[, abscisse],
        estimatedSources[, ordonnee],
        labels = seq_len(nrow(estimatedSources)),
        pos = 3,
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

#                     Plot results

########################################################################################
########################################################################################
########################################################################################


setwd("C++/RESULTS")



typedata="syntheticSA"

numberSourcesDetected=save_plot_statistics(typedata,typedata)
clustersSources=save_dendrogram(typedata,typedata,numberSourcesDetected)
nameSourcesDetected=paste0("sourcesDetected_",typedata)
save_detection(typedata,clustersSources,paste0("sourcesDetected_",typedata),n=50,annealing=T)


typedata="mexicoSA"

numberSourcesDetected=save_plot_statistics(typedata,typedata)
clustersSources=save_dendrogram(typedata,typedata,numberSourcesDetected)
nameSourcesDetected=paste0("sourcesDetected_",typedata)
save_detection(typedata,clustersSources,paste0("sourcesDetected_",typedata),n=50,annealing=T)

typedata="athabascaSA"

numberSourcesDetected=save_plot_statistics(typedata,typedata)
clustersSources=save_dendrogram(typedata,typedata,numberSourcesDetected)
nameSourcesDetected=paste0("sourcesDetected_",typedata)
save_detection(typedata,clustersSources,paste0("sourcesDetected_",typedata),n=50,annealing=T)



########################################################################################
########################################################################################
########################################################################################

#                     Useful plots

########################################################################################
########################################################################################
########################################################################################

data=read.table("../DATA/synthetic_norm.txt",header = T)
real_sources=read.table("../DATA/synthetic_sources2d.txt",header = T)

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
plot(real_sources[,c(1,2)],xlim=c(0,1),ylim=c(0,1),pch=9,col=4,cex=2.5,cex.lab=2.5,cex.axis=2.5)
points(data[,c(1,2)],pch=20)
text(
  real_sources[2:4, 1],
  real_sources[2:4, 2],
  labels = 2:4,  # Number each red point
  pos = 3,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
text(
  real_sources[1, 1],
  real_sources[1, 2],
  labels = 1,  # Number each red point
  pos = 4,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
dev.off()

#Plot the projection on the second plane of the synthetic data set
png(file="synthetic2.png",width = 800,height = 600)
par(mar=c(5,6,4,1)+.1)
plot(real_sources[,c(1,3)],xlim=c(0,1),ylim=c(0,1),pch=9,col=4,cex=2.5,cex.lab=2.5,cex.axis=2.5)
points(data[,c(1,3)],pch=20)
text(
  real_sources[2:4, 1],
  real_sources[2:4, 3],
  labels = 2:4,  # Number each red point
  pos = 3,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
text(
  real_sources[1, 1],
  real_sources[1, 3],
  labels = 1,  # Number each red point
  pos = 4,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
dev.off()

#Plot the projection on the third plane of the synthetic data set
png(file="synthetic3.png",width = 800,height = 600)
par(mar=c(5,6,4,1)+.1)
plot(real_sources[,c(2,3)],xlim=c(0,1),ylim=c(0,1),pch=9,col=4,cex=2.5,cex.lab=2.5,cex.axis=2.5)
points(data[,c(2,3)],pch=20)
text(
  real_sources[1:3, 2],
  real_sources[1:3, 3],
  labels = 1:3,  # Number each red point
  pos = 3,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
text(
  real_sources[4, 2],
  real_sources[4, 3],
  labels = 4,  # Number each red point
  pos = 2,  # Position to the right of the points
  col = 4,
  cex = 4  # Size of text
)
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
