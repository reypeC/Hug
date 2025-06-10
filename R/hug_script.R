rm(list=ls())



# Package names
packages <- c("shiny", "shinyFiles", "dendextend", "fields")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# library(pracma)
# library(igraph)
# library(knitr)
# library(fields)
# library(xtable)
# library(MASS)
# library(fields)
# library(scatterplot3d)
# library(cluster) # pour clusGap
# library(factoextra) # pour fviz_gap_stat
# library(corrplot)
# library(dendextend)#pour afficher que certains clusters du dendrogamme
# palette("R4")

setwd("C:/Users/creype/Documents/Hug/simHug/Hug_finalfixe_FINAL")


########################################################################################
########################################################################################
########################################################################################

#                     Core functions

########################################################################################
########################################################################################
########################################################################################

colfunc<-paste0(colorRampPalette(c("white","royalblue","yellow","red",alpha=0.5))(10),"75")[-c(9,10)]


# Add alpha (transparency) to the colors by modifying the RGB values
colfunc_with_alpha <- rgb(
  col2rgb(colfunc)[1,],
  col2rgb(colfunc)[2,],
  col2rgb(colfunc)[3,],
  alpha = 191,  # 191 corresponds to 75% opacity (255 = fully opaque)
  maxColorValue = 255
)

# Remove the 9th and 10th color from the palette
colfunc_final <- colfunc_with_alpha


normalisation_data=function(dataTrue,bounds=NULL)
{
  numberColumn=ncol(dataTrue)
  normalisedData=dataTrue
  
  if(is.null(bounds))
  {
    for(i in 1:numberColumn)
    {
      scale=max(dataTrue[,i],na.rm=T)-min(dataTrue[,i],na.rm=T)
      lowerBound=min(dataTrue[,i],na.rm=T)-scale
      upperBound=max(dataTrue[,i],na.rm=T)+scale
      bounds=rbind(bounds,c(lowerBound,upperBound))
      
    }
    
  }
  
  # Normalise to [0,1]
  for(i in 1:numberColumn)
  {
    
    normalisedData[,i]=(normalisedData[,i]-bounds[i,1])/(bounds[i,2]-bounds[i,1])
    
  }
  
  
  return(normalisedData)
  
}


reverse_normalisation_data=function(dataTrue,sources,regionOfInterest,bounds=NULL)
{
  numberColumn=ncol(dataTrue)
  normalisedData=dataTrue
  
  if(is.null(bounds))
  {
    for(i in 1:numberColumn)
    {
      scale=max(dataTrue[,i],na.rm=T)-min(dataTrue[,i],na.rm=T)
      lowerBound=min(dataTrue[,i],na.rm=T)-scale
      upperBound=max(dataTrue[,i],na.rm=T)+scale
      bounds=rbind(bounds,c(lowerBound,upperBound))
      
    }
    
  }
  
  reverseSources=sources
  reverseRegionOfInterest=regionOfInterest
  
  for(i in 1:numberColumn)
  {
    
    reverseSources[,i]=reverseSources[,i]*(bounds[i,2]-bounds[i,1])+bounds[i,1]
    reverseRegionOfInterest[dataSources$dim1==i,1]=reverseRegionOfInterest[dataSources$dim1==i,1]*(bounds[i,2]-bounds[i,1])+bounds[i,1]
    reverseRegionOfInterest[dataSources$dim2==i,2]=reverseRegionOfInterest[dataSources$dim2==i,2]*(bounds[i,2]-bounds[i,1])+bounds[i,1]
    
  }
  
  return(list(reverseSources=reverseSources,reverseRegionOfInterest=reverseRegionOfInterest))
  
}

plot_statistic=function(statisticsByPlane,name)
{
  statname=c("n_e","n_a","n_s","n")
  for (plan in 1:length(statisticsByPlane)) {
    
    png(file=paste("statistics_",name,"_",plan,".png",sep=""),width = 1000,height = 200)
    par(mfrow=c(1,2))
    par(mar=c(5,6,4,1)+.1)
    plot(statisticsByPlane[[plan]][,1],type = "l",main = "",ylab = statname[1],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
    abline(v=nrow(statisticsByPlane[[plan]])-1000,col="red",lwd=2,lty=2)
    plot(statisticsByPlane[[plan]][,2],type = "l",main = "",ylab = statname[2],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
    abline(v=nrow(statisticsByPlane[[plan]])-1000,col="red",lwd=2,lty=2)
    dev.off()
    
  }
  
  png(file=paste("statistics_",name,"_strauss.png",sep=""),width = 1000,height = 200)
  par(mfrow=c(1,1))
  par(mar=c(5,6,4,1)+.1)
  plot(statisticsByPlane[[1]][,3],type = "l",main = "",ylab = statname[3],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
  abline(v=nrow(statisticsByPlane[[1]])-1000,col="red",lwd=2,lty=2)
  dev.off()
  
  
  meanstat4=mean(statisticsByPlane[[1]][,4])
  medianstat4=median(statisticsByPlane[[1]][,4])
  
  tableNumber=table(statisticsByPlan[[1]][,4])
  modestat4=as.numeric(names(tableNumber)[which.max(tableNumber)])
  
  
  png(file=paste("statistics_",name,"_number.png",sep=""),width = 1000,height = 300)
  par(mfrow=c(1,2))
  par(mar=c(5,6,4,1)+.1)
  plot(statisticsByPlane[[1]][,4],type = "l",main = "",ylab = statname[4],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
  abline(h=meanstat4,col=palette.colors()[8],lwd=2)
  abline(h=modestat4,col=palette.colors()[3],lwd=2)
  abline(v=nrow(statisticsByPlane[[1]])-1000,col="red",lwd=2,lty=2)
  abline(h=medianstat3,col=palette.colors()[2],lwd=2)
  
  if(length(levels(as.factor(statisticsByPlane[[1]][,4])))==1)
  {
    barplot(tableNumber / nrow(statisticsByPlane[[1]]),col=palette.colors()[6],xlab=statname[4],main="",ylab="Density",cex.lab=2,cex.axis=2)
    legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat4,meanstat4,modestat4),2)),cex=1.5)
  }else
  {
    hist(tableNumber / nrow(statisticsByPlane[[1]]),col=palette.colors()[6],prob=TRUE,nclass=20,xlab=statname[4],main="",cex.lab=2,cex.axis=2)
    abline(v=meanstat4,col=palette.colors()[8],lwd=2)
    abline(v=modestat4,col=palette.colors()[3],lwd=2)
    abline(v=medianstat4,col=palette.colors()[2],lwd=2)
    legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat4,meanstat4,modestat4),2)),cex=1.5)
  }
  
  dev.off()
}

plot_dendrogram=function(sources,name)
{
  clustering = hclust(dist(sources), method = "ward.D2")
  clusters=cutree(clustering, k = 20)
  cutHeight=sort(clustering$height, decreasing = TRUE)[20]

    
    png(paste0("dendrogram_",name, ".png"),width = 1000,height = 1000)
  plot(cut(as.dendrogram(clustering),h=cutHeight)$upper,ylim=c(cutHeight,max(clustering$height)),main="Dendrogram with a maximum of 20 clusters", leaflab = "none")
  dev.off()
  
  return(clustering)
}

compute_proposed_sources=function(sources,numberCluster,clustering)
{
  clusters=cutree(clustering,k=numberCluster)
  proposedSources=aggregate(sources, by = list(Num = clusters), FUN = mean)
  return(proposedSources)
}

saveplt=function(data,sources,detail,real_sources,nbcluster,borne,simsources,MeanDetectedSources,indice,n=50,nbsim)#n1,n2=> grid row,col
{
  
  #clustering
  clustering=kmeans(sources,nbcluster,nstart = 1000)
  estimatedSources=clustering[[2]]
  colnames(estimatedSources)=colnames(sources)
  
  # clustering=hclust(dist(sources))
  # plot(hcSources,labels=F,ylim = c(height_cut, max(clustering$height)),)
  # clusters=cutree(hcSources,k=nbcluster)
  
  sourcesSplit=split(sources,clustering$cluster)
  estimatedSourcesToSave=matrix(NA,ncol=2*ncol(sources),nrow=nbcluster)
  rownames(estimatedSourcesToSave)=names(sourcesSplit)
  colnames(estimatedSourcesToSave)=rep(NA,2*ncol(sources))
  
  for(numSource in 1:length(sourcesSplit))
  {
    for (i in 1:ncol(sources)) {
      estimatedSourcesToSave[numSource,2*(i)]=mean(abs(sourcesSplit[[numSource]][,i]-estimatedSources[numSource,i]))
      estimatedSourcesToSave[numSource,2*(i-1)+1]=estimatedSources[numSource,i]
      colnames(estimatedSourcesToSave)[c(2*(i-1)+1,2*i)]=paste0(colnames(sources)[i],c(""," mean"))
    }
    
    
  }
  
  write.table( estimatedSourcesToSave, paste0(MeanDetectedSources,".txt"), row.names=T, col.names=T, append=F )
  print(xtable(estimatedSourcesToSave, type = "latex"), file = paste0(MeanDetectedSources,".tex"))
  
  namecol=c()
  for (j in colnames(sources)) {
    namecol=c(namecol,paste(j,".median",sep = ""),paste(j,".mean",sep = ""),paste(j,".sd",sep = ""))
  }
  write.table(t(namecol) , paste0(MeanDetectedSources,"_resume",".txt"), row.names=F, col.names=F, append=F )
  
  
  medianpattern=NULL
  for (numcluster in 1:nbcluster) {
    sourcemedian=NULL
    rowmedianpattern=c()
    for (dim in 1:ncol(sources)) {
      sourcemedian=c(sourcemedian,round(c(quantile(sources[clustering$cluster==numcluster,dim])[3],mean(sources[clustering$cluster==numcluster,dim]),sd(sources[clustering$cluster==numcluster,dim])),2))
      rowmedianpattern=c(rowmedianpattern,quantile(sources[clustering$cluster==numcluster,dim])[3])
    }
    write.table(t(sourcemedian) , paste0(MeanDetectedSources,"_resume",".txt"), row.names=F, col.names=F, append=T )
    
    medianpattern=rbind(medianpattern,rowmedianpattern)
  }
  # print(xtable(rowmedianpattern, type = "latex"), file = paste0(MeanDetectedSources,"_resume",".tex"))
  
  print("Detected sources saved!")
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
      
      kde <- kde2d(sources[,abscisse], sources[,ordonnee], n = n,lims = c(c(minabs, maxabs),  c(minord, maxord)))
      #
      kde$z <- (kde$z)/nbsim
      
      png(file=paste("lvl_planeSA_",plan,'_',indice,'.png',sep=''),width = 800,height = 600)
      par(mar=c(5,6,4,1)+.1)
      filled.contour(
        kde,
        xlim = c(minabs, maxabs), ylim = c(minord, maxord),
        color.palette = function(n) colfunc,
        levels = seq(min(kde$z), max(kde$z), length.out = 9),
        plot.axes = {
          abline(h = seq(minord, maxord, length.out = n), col = "gray")  # Horizontal grid lines
          abline(v = seq(minabs, maxabs, length.out = n), col = "gray")
          points(estimatedSources[,c(abscisse,ordonnee)],pch=19,col='darkgreen',cex=3)
          points(real_sources_plane[, c(1, 2)], pch = 17, col = "red",cex=3)  # Add points for real sources
          points(data[, c(abscisse, ordonnee)], pch = 19, col = "black",cex=1)  # Add points for data
          # Add numbers for red points
          text(
            estimatedSources[, abscisse],
            estimatedSources[, ordonnee],
            labels = seq_len(nrow(estimatedSources)),  # Number each red point
            pos = 3,  # Position to the right of the points
            col = "darkgreen",
            cex = 3  # Size of text
          )
          axis(1)  # Add x-axis
          axis(2)  # Add y-axis
        },
        xlab = colnames(sources)[[abscisse]],
        ylab = colnames(sources)[[ordonnee]]
      )
      dev.off()
      
      # ####
      #
      # plan=plan+1
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
      # # Plot using filled.contour
      
      png(file=paste("lvlb_planeSA_",plan,'_',indice,'.png',sep=''),width = 800,height = 600)
      par(mar=c(5,6,4,1)+.1)
      filled.contour(
        x = x_breaks[-50],
        y = y_breaks[-50],
        z = z,
        level=seq(min(z), max(z), length.out = 9),
        xlim = c(minabs, maxabs), ylim = c(minord, maxord),
        color.palette = function(n) colfunc,
        plot.axes = {
          abline(h = seq(minord, maxord, length.out = 50), col = "gray")  # Horizontal grid lines
          abline(v = seq(minabs, maxabs, length.out = 50), col = "gray")
          points(estimatedSources[,c(abscisse,ordonnee)],pch=19,col='darkgreen',cex=3)
          points(real_sources_plane[, c(1, 2)], pch = 17, col = "red",cex=3)  # Add points for real sources
          points(data[, c(abscisse, ordonnee)], pch = 19, col = "black",cex=1)  # Add points for data
          # Add numbers for red points
          text(
            estimatedSources[, abscisse],
            estimatedSources[, ordonnee],
            labels = seq_len(nrow(estimatedSources)),  # Number each red point
            pos = 3,  # Position to the right of the points
            col = "darkgreen",
            cex = 3  # Size of text
          )
          axis(1)  # Add x-axis
          axis(2)  # Add y-axis
        },
        xlab = colnames(sources)[[abscisse]],
        ylab = colnames(sources)[[ordonnee]]
      )
      dev.off()
      #
      # ###
      # dat=data.frame(Li=c(900,3000,520,6000),Na=c(80000,100000,15000,22000),Mg=c( 4000,  9000, 22000  ,   40000  ),K=c(1700,5200,8000,17000),Ca=c(11000,32000,27000,60000 ))
      
      png(file = paste("lvlc_planeSA_", plan, '_', indice, '.png', sep = ''), width = 800, height = 600)
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
      
      # segments(dat[1,abscisse],dat[1,ordonnee],dat[2,abscisse],dat[1,ordonnee],col = "blue",lwd=2)
      # segments(dat[1,abscisse],dat[2,ordonnee],dat[2,abscisse],dat[2,ordonnee],col = "blue",lwd=2)
      #
      # segments(dat[1,abscisse],dat[1,ordonnee],dat[1,abscisse],dat[2,ordonnee],col = "blue",lwd=2)
      # segments(dat[2,abscisse],dat[2,ordonnee],dat[2,abscisse],dat[1,ordonnee],col = "blue",lwd=2)
      #
      # segments(dat[3,abscisse],dat[3,ordonnee],dat[4,abscisse],dat[3,ordonnee],col = "red",lwd=2)
      # segments(dat[3,abscisse],dat[4,ordonnee],dat[4,abscisse],dat[4,ordonnee],col = "red",lwd=2)
      #
      # segments(dat[3,abscisse],dat[3,ordonnee],dat[3,abscisse],dat[4,ordonnee],col = "red",lwd=2)
      # segments(dat[4,abscisse],dat[4,ordonnee],dat[4,abscisse],dat[3,ordonnee],col = "red",lwd=2)
      
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
      
      
      print(paste("Clusters and level sets saved for the plane number: ",plan,sep=""))
    }
  }
  
  
}


########################################################################################
########################################################################################
########################################################################################

#                     Applications

########################################################################################
########################################################################################
########################################################################################





######



############################### Tetra

############################### Pinti

############################### Athabasca



###################################

# setwd("C:/Users/creype/Documents/Hug/simHug/Hug_finalfixe")




setwd("RESULTS")
setwd("..")
setwd("RESULTStrue")

statname=c("n_e","n_a","n_r","n")
listetype="mean"

listetype="simtetra"
listetype=paste0(listetype,"pinti")
listetype=paste0(listetype,"athab")
listetype=paste0(listetype,"SA")
# indice=0

typedata="tetra"
typedata="pinti"
typedata="athab"
# type=listetype
# numero="5"
for (type in listetype) {
  
  
  for(numero in c("5"))
  {
    # numero=""
    # indice=indice+1
    
    sources=read.table(paste0("sources_final","_",type,numero,".txt"),header=T)
    
    
    
    stats=read.table(paste0("statistics","_",type,numero,".txt"),header=T)
    statsFinal=read.table(paste0("statistics_final","_",type,numero,".txt"),header=F)
    colnames(statsFinal)=colnames(stats)
    
    statsSplit=split(stats,factor(paste0(stats[,5],stats[,6])))
    statsFinalSplit=split(statsFinal,factor(paste0(statsFinal[,5],statsFinal[,6])))
    statsplan=list()
    for(i in 1:length(statsSplit))
    {
      statsplan[[names(statsSplit)[i]]]=rbind(statsSplit[[i]],statsFinalSplit[[i]][seq(from=1,to=100000,by=100),])
      
    }
    
    cumSumN=cumsum(statsFinalSplit[[1]][,4])
    
    numberSourcesToSave=NULL
    for(i in seq(from=1,to=100000,by=100))
    {
      numberSourcesToSave=c(numberSourcesToSave,(cumSumN[i]-statsFinalSplit[[1]][i,4]+1):cumSumN[i])
    }
    sources=sources[numberSourcesToSave,]
    
    clustering=hclust(dist(sources),method = "ward.D2")
    clusters=cutree(clustering,k=20)
    
    # Find the height that corresponds to this cut (i.e., where to "zoom" in)
    cut_height <- sort(clustering$height, decreasing = TRUE)[20]
    
    plot(cut(as.dendrogram(clustering),h=0)$upper,ylim=c(cut_height,max(clustering$height)),main="Dendrogram with a maximum of 20 clusters")
    
    png(file=paste("dendro_",type,"_",numero,".png",sep=""),width = 1000,height = 1000)
    plot(cut(as.dendrogram(clustering),h=cut_height)$upper,ylim=c(cut_height,max(clustering$height)),main="Dendrogram with a maximum of 20 clusters", leaflab = "none")
    dev.off()
    
    # gap.stat <- clusGap(sources, FUNcluster = kmeans, K.max = 20)
    #
    # png(file=paste("gapStat_",type,"_",numero,".png",sep=""),width = 1000,height = 1000)
    # fviz_gap_stat(gap.stat)
    # dev.off()
    #
    # nclusGap <- maxSE(f         = gap.stat$Tab[,"gap"],
    #                   SE.f      = gap.stat$Tab[,"SE.sim"],
    #                   method    = "firstSEmax",
    #                   SE.factor = 1)
    
    if(typedata=="pinti")
    {
      dataTrue=read.table("../../DATA/pinti.txt",header = T)
      
      data=read.table("../../DATA/pintinorm.txt",header = T)
      dataSources=read.table("../../DATA/pinti_sources2d.txt",header = T)
      
      res=reverse(dataTrue,data,sources,dataSources,positif=F)
      sources=res[["sources"]]
      dataSources=res[["dataSources"]]
      borne=res[["borne"]]
      data=dataTrue
    }else if(typedata=="athab")
    {
      dataTrue=read.table("../../DATA/athabasca_8d.txt",header = T)
      
      data=read.table("../../DATA/athabasca_8dnorm2.txt",header = T)
      dataSources=read.table("../../DATA/athabasca_sources2d.txt",header = T)
      
      res=reverse(dataTrue,data,sources,dataSources,positif=T)
      sources=res[["sources"]]
      dataSources=res[["dataSources"]]
      borne=res[["borne"]]
      data=dataTrue
    }else{
      data=read.table("../../DATA/tetraedre_datanorm.txt",header = T)
      dataSources=read.table("../../DATA/tetraedre_sources2d.txt",header = T)
      
      alldim=as.numeric(unique(c(levels(factor(dataSources[["dim1"]])),levels(factor(dataSources[["dim2"]])))))
      borne=data.frame(minCoord=rep(0,length(alldim)),maxCoord=rep(1,length(alldim)),dim=alldim)
    }
    
    
    for (plan in 1:length(statsplan)) {
      dstat1=density(statsplan[[plan]][,1],kernel="epanechnikov")
      pstat1=which.max(dstat1$y)
      meanstat1=mean(statsplan[[plan]][,1])
      medianstat1=median(statsplan[[plan]][,1])
      modestat1=dstat1$x[pstat1]
      
      dstat2=density(statsplan[[plan]][,2],kernel="epanechnikov")
      pstat2=which.max(dstat2$y)
      meanstat2=mean(statsplan[[plan]][,2])
      medianstat2=median(statsplan[[plan]][,2])
      modestat2=dstat2$x[pstat2]
      
      dstat3=density(statsplan[[plan]][,3],kernel="epanechnikov")
      pstat3=which.max(dstat3$y)
      meanstat3=mean(statsplan[[plan]][,3])
      medianstat3=median(statsplan[[plan]][,3])
      modestat3=dstat3$x[pstat3]
      
      dstat4=density(statsplan[[plan]][,4],kernel="epanechnikov")
      pstat4=which.max(dstat4$y)
      meanstat4=mean(statsplan[[plan]][,4])
      medianstat4=median(statsplan[[plan]][,4])
      modestat4=dstat4$x[pstat4]
      
      png(file=paste("truestat_",type,"_",plan,"_",numero,".png",sep=""),width = 1000,height = 1000)
      par(mfrow=c(4,2))
      
      plot(statsplan[[plan]][,1],type = "l",main = "",ylab = statname[1],col=palette.colors()[6])
      abline(h=meanstat1,col=palette.colors()[8],lwd=2)
      abline(h=modestat1,col=palette.colors()[3],lwd=2)
      abline(h=medianstat1,col=palette.colors()[2],lwd=2)
      
      hist(statsplan[[plan]][,1],col=palette.colors()[6],prob=TRUE,ylim=c(0,1.25*max(dstat1$y)),nclass=20,xlab=statname[1],main="")
      lines(dstat1,col="red")
      points(dstat1$x[pstat1],dstat1$y[pstat1],col=palette.colors()[3],pch="+",cex=3)
      abline(v=meanstat1,col=palette.colors()[8],lwd=2)
      abline(v=modestat1,col=palette.colors()[3],lwd=2)
      abline(v=medianstat1,col=palette.colors()[2],lwd=2)
      legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat1,meanstat1,dstat1$x[pstat1]),4)),cex=1.5)
      
      
      
      plot(statsplan[[plan]][,2],type = "l",main = "",ylab = statname[2],col=palette.colors()[6])
      abline(h=meanstat2,col=palette.colors()[8],lwd=2)
      abline(h=modestat2,col=palette.colors()[3],lwd=2)
      abline(h=medianstat2,col=palette.colors()[2],lwd=2)
      
      hist(statsplan[[plan]][,2],col=palette.colors()[6],prob=TRUE,ylim=c(0,1.25*max(dstat2$y)),nclass=20,xlab=statname[2],main="")
      lines(dstat2,col="red")
      points(dstat2$x[pstat2],dstat2$y[pstat2],col=palette.colors()[3],pch="+",cex=3)
      abline(v=meanstat2,col=palette.colors()[8],lwd=2)
      abline(v=modestat2,col=palette.colors()[3],lwd=2)
      abline(v=medianstat2,col=palette.colors()[2],lwd=2)
      legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat2,meanstat2,dstat2$x[pstat2]),4)),cex=1.5)
      
      plot(statsplan[[plan]][,3],type = "l",main = "",ylab = statname[3],col=palette.colors()[6])
      abline(h=meanstat3,col=palette.colors()[8],lwd=2)
      abline(h=modestat3,col=palette.colors()[3],lwd=2)
      abline(h=medianstat3,col=palette.colors()[2],lwd=2)
      
      hist(statsplan[[plan]][,3],col=palette.colors()[6],prob=TRUE,ylim=c(0,1.25*max(dstat3$y)),nclass=20,xlab=statname[3],main="")
      lines(dstat3,col="red")
      points(dstat3$x[pstat3],dstat3$y[pstat3],col=palette.colors()[3],pch="+",cex=3)
      abline(v=meanstat3,col=palette.colors()[8],lwd=2)
      abline(v=modestat3,col=palette.colors()[3],lwd=2)
      abline(v=medianstat3,col=palette.colors()[2],lwd=2)
      legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat3,meanstat3,dstat3$x[pstat3]),4)),cex=1.5)
      
      plot(statsplan[[plan]][,4],type = "l",main = "",ylab = statname[4],col=palette.colors()[6])
      abline(h=meanstat3,col=palette.colors()[8],lwd=2)
      abline(h=modestat3,col=palette.colors()[3],lwd=2)
      abline(h=medianstat3,col=palette.colors()[2],lwd=2)
      
      hist(statsplan[[plan]][,4],col=palette.colors()[6],prob=TRUE,ylim=c(0,1.25*max(dstat4$y)),nclass=20,xlab=statname[4],main="")
      lines(dstat4,col="red")
      points(dstat4$x[pstat4],dstat4$y[pstat4],col=palette.colors()[4],pch="+",cex=3)
      abline(v=meanstat4,col=palette.colors()[8],lwd=2)
      abline(v=modestat4,col=palette.colors()[3],lwd=2)
      abline(v=medianstat4,col=palette.colors()[2],lwd=2)
      legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat4,meanstat4,dstat4$x[pstat4]),4)),cex=1.5)
      
      dev.off()
      
      png(file=paste("truestat_byPlan_",type,"_",plan,"_",numero,".png",sep=""),width = 1000,height = 200)
      par(mfrow=c(1,2))
      par(mar=c(5,6,4,1)+.1)
      plot(statsplan[[plan]][,1],type = "l",main = "",ylab = statname[1],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
      abline(v=nrow(statsplan[[plan]])-1000,col="red",lwd=2,lty=2)
      plot(statsplan[[plan]][,2],type = "l",main = "",ylab = statname[2],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
      abline(v=nrow(statsplan[[plan]])-1000,col="red",lwd=2,lty=2)
      # plot(statsplan[[plan]][,3],type = "l",main = "",ylab = statname[3],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
      dev.off()
      
      png(file=paste("truestat_byPlanStraussGlobal_",type,"_",plan,"_",numero,".png",sep=""),width = 1000,height = 200)
      par(mfrow=c(1,1))
      par(mar=c(5,6,4,1)+.1)
      plot(statsplan[[plan]][,3],type = "l",main = "",ylab = statname[3],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
      abline(v=nrow(statsplan[[plan]])-1000,col="red",lwd=2,lty=2)
      dev.off()
      
      png(file=paste("truestat_byPlanStrauss_",type,"_",plan,"_",numero,".png",sep=""),width = 1000,height = 300)
      par(mfrow=c(1,2))
      par(mar=c(5,6,4,1)+.1)
      plot(statsplan[[plan]][,3],type = "l",main = "",ylab = statname[3],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
      abline(h=meanstat3,col=palette.colors()[8],lwd=2)
      abline(h=modestat3,col=palette.colors()[3],lwd=2)
      abline(h=medianstat3,col=palette.colors()[2],lwd=2)
      abline(v=nrow(statsplan[[plan]])-1000,col="red",lwd=2,lty=2)
      
      if(length(levels(as.factor(statsplan[[plan]][,3])))==1)
      {
        barplot(table(statsplan[[plan]][,3])/table(statsplan[[plan]][,3]),col=palette.colors()[6],xlab=statname[3],main="",ylab="Density",cex.lab=2,cex.axis=2)
        legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat3,meanstat3,dstat3$x[pstat3]),2)),cex=1.5)
      }else
      {
        hist(statsplan[[plan]][,3],col=palette.colors()[6],prob=TRUE,ylim=c(0,1.25*max(dstat3$y)),nclass=20,xlab=statname[3],main="")
        lines(dstat3,col="red")
        points(dstat3$x[pstat3],dstat3$y[pstat3],col=palette.colors()[3],pch="+",cex=3)
        abline(v=meanstat3,col=palette.colors()[8],lwd=2)
        abline(v=modestat3,col=palette.colors()[3],lwd=2)
        abline(v=medianstat3,col=palette.colors()[2],lwd=2)
        legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat3,meanstat3,dstat3$x[pstat3]),2)),cex=1.5)
      }
      
      
      dev.off()
      
      png(file=paste("truestat_byPlanNumber_",type,"_",plan,"_",numero,".png",sep=""),width = 1000,height = 300)
      par(mfrow=c(1,2))
      par(mar=c(5,6,4,1)+.1)
      plot(statsplan[[plan]][,4],type = "l",main = "",ylab = statname[4],col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
      abline(h=meanstat4,col=palette.colors()[8],lwd=2)
      abline(h=modestat4,col=palette.colors()[3],lwd=2)
      abline(v=nrow(statsplan[[plan]])-1000,col="red",lwd=2,lty=2)
      abline(h=medianstat3,col=palette.colors()[2],lwd=2)
      
      if(length(levels(as.factor(statsplan[[plan]][,4])))==1)
      {
        barplot(table(statsplan[[plan]][,4])/table(statsplan[[plan]][,4]),col=palette.colors()[6],xlab=statname[4],main="",ylab="Density",cex.lab=2,cex.axis=2)
        legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat4,meanstat4,dstat4$x[pstat4]),2)),cex=1.5)
      }else
      {
        hist(statsplan[[plan]][,4],col=palette.colors()[6],prob=TRUE,ylim=c(0,1.25*max(dstat4$y)),nclass=20,xlab=statname[4],main="",cex.lab=2,cex.axis=2)
        lines(dstat4,col="red")
        points(dstat4$x[pstat4],dstat4$y[pstat4],col=palette.colors()[4],pch="+",cex=3)
        abline(v=meanstat4,col=palette.colors()[8],lwd=2)
        abline(v=modestat4,col=palette.colors()[3],lwd=2)
        abline(v=medianstat4,col=palette.colors()[2],lwd=2)
        legend("topright",legend = paste0(c("median = ","mean = ","mode = "),round(c(medianstat4,meanstat4,dstat4$x[pstat4]),2)),cex=1.5)
      }
      
      dev.off()
      
      
    }
    MeanDetectedSources=paste0("proposedMeanSources_",paste0(type,numero))
    simsources=paste0("quantileSources_",paste0(type,numero),".txt")
    if((mean(stats[,4])<=10)&&(mean(stats[,4])>2))
    {
      saveplt(data,sources,statsFinal,dataSources,round(mean(statsFinal[,4])),borne,simsources,MeanDetectedSources,paste0(type,numero),n=50,length(seq(from=1,to=100000,by=100)))
      
    }
    
  }
}

