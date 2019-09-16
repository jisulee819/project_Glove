library(CCA)
library(reshape2)
library(ggplot2)
library(HDCI)

multiplot <- function(..., plotlist=NULL, file, cols=2, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#subjectpool <- c(1,2,5,7,8,11:15,18:21,23,26:38,42,44)
#nsubject <- length(subjectpool)
#for (subjectID in subjectpool){
#  
#  if (subjectID < 10){
#    filename <- paste0("GA0",subjectID)
#  } else {filename <- paste0("GA",subjectID)
#  }
#  
#  dat1 <- read.csv(paste0("/Volumes/clmnlab/GA/behavior_data/behavior_csv/",filename,"_sensor_GAp.csv"))
#  dat2 <- read.csv(paste0("/Volumes/clmnlab/GA/behavior_data/behavior_csv/",filename,"_sensor_GBp.csv"))
#}

###

library(R.matlab)


subjectpool <- c(1,2,5,7,8,11:15,18:21,23,26:38,42,44)
nsubject <- length(subjectpool)

for (subjectID in subjectpool){

    if (subjectID < 10){
      filename <- paste0("GA0",subjectID)
    } else {filename <- paste0("GA",subjectID)
    }
    
    dat1 <- readMat(paste0("/Volumes/clmnlab/GA/behavior_data/",filename,"/",filename,"-fmri_preproc.mat"))
    sensors <- dat1$filteredData[,c(1:87300)]
    cursors <- dat1$allXY[,(1:87300)]
    
    m.sensors <- matrix(NA, nrow=14, ncol=2910)
    m.cursors <- matrix(NA, nrow=2, ncol=2910)
    
    for (j in 1:2910){
      sub.sensors <- sensors[, (30*j-29) : (30*j)]
      sub.cursors <- cursors[, (30*j-29) : (30*j)]
      m.sensors[,j] <- apply(sub.sensors, 1, mean) 
      m.cursors[,j] <- apply(sub.cursors, 1, mean) 
    }
    
    dist.cursors <- rep(NA, dim(m.cursors)[2])
    
    
    diff.pool <- c(1:2910)[-seq(10,2910,10)]
    
    for (j in diff.pool){
      if (j == 1){
        d <- m.sensors[,j+1] - m.sensors[,j]
      } else {
        d2 <- m.sensors[,j+1] - m.sensors[,j]
        d <- cbind(d, d2)
      }
    }
    d.sensors <- matrix(NA, nrow=14, ncol=length(diff.pool))
    d.sensors <- d
    dim(d.sensors)
    d.sensors3 <- t(d.sensors)
    
    stat.data <- data.frame("run" = rep(c(1:3), each = 9*97), 
                            "trial" = rep(c(1:97), each = 9, 3),
                            "number" = rep(c(1:9), 97*3),
                            d.sensors3)
    
    dim(d.sensors3)
    
    targetID <- dat1$targetID[1:(97*3)]
    a <- data.frame(targetID)
    a$target_pos1 <- ifelse(targetID == 1 | targetID == 21, -160, 160)
    a$target_pos2 <- ifelse(targetID == 1 | targetID == 5, 160, -160)
    #c(a$target_pos1[a$targetID == 25][1], a$target_pos2[a$targetID == 25][1])
    
    target_pos <- rbind(a$target_pos1, a$target_pos2)
    dim(target_pos)
    
    for (i in 1:dim(m.cursors)[2]){
      d <- rbind(m.cursors[,i], target_pos[,((i-1)%/%10)+1])
      d2 <- dist(d, method="euclidean")
      dist.cursors[i] <- d2
      #print((i%/%10)+1)
    }
    
    d.dist <- matrix(NA, nrow=1, ncol=length(diff.pool))
    for (j in diff.pool){
      if (j == 1){
        d <- dist.cursors[j+1] - dist.cursors[j]
      } else {
        d2 <- dist.cursors[j+1] - dist.cursors[j]
        d <- append(d, d2)
      }
    }
    
    d.dist <- matrix(d)
    #d.dist 확인
    #1           2           3           4           5           6           7           8           9 
    #-47.556094 -120.223473  -64.265624  -31.993393  -11.055892  -14.447085   -2.413967   -3.191865   -5.190311 
    
    
    stat.data <- data.frame(stat.data, "d.dist" = d.dist)
    colnames(stat.data) = c("run", "trial", "number", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                            "10", "11", "12", "13", "14", "d.dist")
  
    
    ### Path indicator
    # the first 1 in "NAN" and the others are 8 times repeatition of 12 paths sequence
    # paths are repeated for 9 times each
    targetID <- dat1$targetID[1:(97*3)]
    pathID <- NA
    for (j in 2:96){
      pathID[j-1] <- paste0(targetID[j],"to",targetID[j+1])  
    }
    path <- rep(pathID[1:12])
    rep.path<- c("Nan",rep(pathID[1:12], 8))
    
    stat.data$path <- rep(rep(rep.path , each = 9), 3)
    stat.data$ID <- subjectID
    
    if (subjectID == 1){
      GA.stat.data <- stat.data 
    } else {
      GA.stat.data <- rbind(GA.stat.data, stat.data)
    }
    
    for (k in 1:3){
      for (l in 1:12){
      
        subdata <- stat.data[stat.data$run == k,]
        subdata2 <- subdata[subdata$path == path[l],]
        
        cca.data <- cc(subdata2[,c(4:17)], matrix(subdata2$d.dist))
        subdata2$`1`
        names(subdata2[,c(4:17)])
        
        
        cordata <- data.frame("run" = rep(k), "path" = rep(path[l]), "ID" = rep(subjectID), "cor" = cca.data$cor)
        if ((subjectID == 1) && (k == 1) && (l == 1)){
          GA.cor.data <- cordata
        } else {
          GA.cor.data <- rbind(GA.cor.data, cordata)
        }
        
        coefdata <- matrix(NA)
        coefdata <- data.frame("run" = rep(k, 14), "path" = rep(path[l], 14), "path2" = rep(path2[l], 14), "ID" = rep(subjectID, 14), "coef" = cca.data$xcoef)
        dataname <- paste0("run",k, "-path",l, "-ID",subjectID)
        
        if ((subjectID == 1) && (k == 1) && (l == 1)){
          GA.coef.data <- coefdata
        } else {
          GA.coef.data <- rbind(GA.coef.data, coefdata)
        }



}}

}






subjectpool <- c(1,2,5,7,8,11:15,18:21,23,26:38,42,44)
nsubject <- length(subjectpool)

for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("GA0",subjectID)
  } else {filename <- paste0("GA",subjectID)
  }
  
  dat1 <- readMat(paste0("/Volumes/clmnlab/GA/behavior_data/",filename,"/",filename,"-refmri_preproc.mat"))
  sensors <- dat1$filteredData[,c(1:87300)]
  cursors <- dat1$allXY[,(1:87300)]
  
  m.sensors <- matrix(NA, nrow=14, ncol=2910)
  m.cursors <- matrix(NA, nrow=2, ncol=2910)
  
  for (j in 1:2910){
    sub.sensors <- sensors[, (30*j-29) : (30*j)]
    sub.cursors <- cursors[, (30*j-29) : (30*j)]
    m.sensors[,j] <- apply(sub.sensors, 1, mean) 
    m.cursors[,j] <- apply(sub.cursors, 1, mean) 
  }
  
  dist.cursors <- rep(NA, dim(m.cursors)[2])
  
  
  diff.pool <- c(1:2910)[-seq(10,2910,10)]
  
  for (j in diff.pool){
    if (j == 1){
      d <- m.sensors[,j+1] - m.sensors[,j]
    } else {
      d2 <- m.sensors[,j+1] - m.sensors[,j]
      d <- cbind(d, d2)
    }
  }
  d.sensors <- matrix(NA, nrow=14, ncol=length(diff.pool))
  d.sensors <- d
  dim(d.sensors)
  d.sensors3 <- t(d.sensors)
  
  stat.data <- data.frame("run" = rep(c(1:3), each = 9*97), 
                          "trial" = rep(c(1:97), each = 9, 3),
                          "number" = rep(c(1:9), 97*3),
                          d.sensors3)
  
  dim(d.sensors3)
  
  
  a <- data.frame(targetID)
  a$target_pos1 <- ifelse(targetID == 1 | targetID == 21, -160, 160)
  a$target_pos2 <- ifelse(targetID == 1 | targetID == 5, 160, -160)
  #c(a$target_pos1[a$targetID == 25][1], a$target_pos2[a$targetID == 25][1])
  
  target_pos <- rbind(a$target_pos1, a$target_pos2)
  dim(target_pos)
  
  for (i in 1:dim(m.cursors)[2]){
    d <- rbind(m.cursors[,i], target_pos[,((i-1)%/%10)+1])
    d2 <- dist(d, method="euclidean")
    dist.cursors[i] <- d2
    #print((i%/%10)+1)
  }
  
  d.dist <- matrix(NA, nrow=1, ncol=length(diff.pool))
  for (j in diff.pool){
    if (j == 1){
      d <- dist.cursors[j+1] - dist.cursors[j]
    } else {
      d2 <- dist.cursors[j+1] - dist.cursors[j]
      d <- append(d, d2)
    }
  }
  
  d.dist <- matrix(d)
  #d.dist 확인
  #1           2           3           4           5           6           7           8           9 
  #-47.556094 -120.223473  -64.265624  -31.993393  -11.055892  -14.447085   -2.413967   -3.191865   -5.190311 
  
  
  stat.data <- data.frame(stat.data, "d.dist" = d.dist)
  colnames(stat.data) = c("run", "trial", "number", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                          "10", "11", "12", "13", "14", "d.dist")
  
  
  ### Path indicator
  # the first 1 in "NAN" and the others are 8 times repeatition of 12 paths sequence
  # paths are repeated for 9 times each
  targetID <- dat1$targetID[1:(97*3)]
  pathID <- NA
  for (j in 2:96){
    pathID[j-1] <- paste0(targetID[j],"to",targetID[j+1])  
  }
  path <- rep(pathID[1:12])
  
  path2 <- ifelse(path == "1to5" | path == "5to1" | path == "21to25" | path == "25to21", "hor",
                  ifelse(path == "1to21" | path == "21to1" | path == "5to25" | path == "25to5", "ver", 
                         ifelse(path == "1to25" | path == "25to1" | path == "5to21" | path == "21to5", "dia", "NA")))
  
  
  rep.path<- c("Nan",rep(pathID[1:12], 8))
  
  stat.data$path <- rep(rep(rep.path , each = 9), 3)
  stat.data$ID <- subjectID
  
  if (subjectID == 1){
    GB.stat.data <- stat.data 
  } else {
    GB.stat.data <- rbind(GB.stat.data, stat.data)
  }
  
  for (k in 1:3){
    for (l in 1:12){
      
      subdata <- stat.data[stat.data$run == k,]
      subdata2 <- subdata[subdata$path == path[l],]
      
      cca.data <- cc(subdata2[,c(4:17)], matrix(subdata2$d.dist))
      subdata2$`1`
      names(subdata2[,c(4:17)])
      
      
      cordata <- data.frame("run" = rep(k), "path" = rep(l), "ID" = rep(subjectID), "cor" = cca.data$cor)
      if ((subjectID == 1) && (k == 1) && (l == 1)){
        GB.cor.data <- cordata
      } else {
        GB.cor.data <- rbind(GB.cor.data, cordata)
      }
      
      coefdata <- matrix(NA)
      coefdata <- data.frame("run" = rep(k, 14), "path" = rep(l, 14), "ID" = rep(subjectID, 14), "coef" = cca.data$xcoef)
      dataname <- paste0("run",k, "-path",l, "-ID",subjectID)
      
      if ((subjectID == 1) && (k == 1) && (l == 1)){
        GB.coef.data <- coefdata
      } else {
        GB.coef.data <- rbind(GB.coef.data, coefdata)
      }
      
      
      
    }}
  
}

#### GA - GB compile
table(GB.coef.data$run)
dim(GB.coef.data) #3*12*14*30
GB.coef.data$run <- ifelse( GB.coef.data$run == 1, 4, 
                            ifelse( GB.coef.data$run == 2, 5, 
                                    ifelse(GB.coef.data$run == 3, 6, NA)))
GB.cor.data$run <- ifelse( GB.cor.data$run == 1, 4, 
                            ifelse( GB.cor.data$run == 2, 5, 
                                    ifelse(GB.cor.data$run == 3, 6, NA)))

GB.stat.data$run <- ifelse( GB.stat.data$run == 1, 4, 
                              ifelse( GB.stat.data$run == 2, 5, 
                                      ifelse(GB.stat.data$run == 3, 6, NA)))

coef.data <- rbind(GA.coef.data, GB.coef.data)
cor.data <- rbind(GA.cor.data, GB.cor.data)
stat.data <- rbind(GA.stat.data, GB.stat.data)

#### ploting the correlation
plot.cor.data <- tapply(cor.data$cor, list(cor.data$run, cor.data$path), mean)
plot.cor.data.ID <- tapply(cor.data$cor, list(cor.data$run, cor.data$path, cor.data$ID), mean)

cor.sd <- matrix(NA, 6, 12)
for (i in 1:6){
  for (j in 1:12){
    cor.sd[i,j] <- print(sd(plot.cor.data.ID[i,j,]))
  }}
cor.sd

plot.cor.data <- reshape2::melt(plot.cor.data)
cor.sd <- reshape2::melt(cor.sd)
plot.cor.data$sd <- cor.sd$value


for (j in 1:12){
  p <- ggplot(plot.cor.data[plot.cor.data$Var2 == j,], aes(x=Var1, y=value)) +
    geom_point() +
    geom_line() +
    scale_y_continuous("Correlation", limit = c(0.85, 1.03), breaks=seq(0.85, 1.05, 0.01)) +
    scale_x_continuous(paste0("Path ",j),breaks=c(1:6)) +
    theme(legend.position="none") +
    coord_cartesian(ylim = c(0.87, 1)) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05))

  assign(paste0("p",j), p)
}

png(paste0("/Users/jisu/Documents/Glove/cca_cor.png"), width = 4000, height = 2500, res=220)
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, cols=4)
dev.off()


colnames(plot.cor.data) = c("Run", "Path", "Value")
plot.cor.data$Path2 <- as.character(plot.cor.data$Path)
p <- ggplot(plot.cor.data, aes(x=Run, y=Value, group=Path2)) +
    geom_point(aes(color=Path2)) +
    geom_line(aes(color=Path2)) +
    scale_y_continuous("Correlation", limit = c(0.89, 1), breaks=seq(0.90, 1.00, 0.01)) +
    scale_x_continuous(paste0("Run"),breaks=c(1:6)) +
    scale_color_discrete(name = "Path",
                         breaks = c(1:12),
                         label = c("1","2","3","4","5","6","7","8","9","10","11","12"))
  

png(paste0("/Users/jisu/Documents/Glove/cca_cor_group.png"), width = 2000, height = 2000, res=220)
p
dev.off()


#### cosine distance from PCAs
#### take out PCAs! 

for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("GA0",subjectID)
  } else {filename <- paste0("GA",subjectID)
  }
  
  dat1 <- readMat(paste0("/Volumes/clmnlab/GA/behavior_data/",filename,"/",filename,"-fmri_preproc.mat"))
  dat1$A
  hor <- matrix(NA, nrow=14, ncol=2)
  hor[1:14,1] <- subjectID
  hor[1:14,2] <- dat1$A[,1]
  
  ver <- matrix(NA, nrow=14, ncol=2)
  ver[1:14,1] <- subjectID
  ver[1:14,2] <- dat1$A[,2]
  
  dia <- matrix(NA, nrow=14, ncol=2)
  dia[1:14,1] <- subjectID
  dia[1:14,2] <- 0.5*dat1$A[,1] + 0.5*dat1$A[,2]
  
  if (subjectID == 1){
    PCA1 <- hor
    PCA2 <- ver
    PCA3 <- dia
  } else { 
    PCA1 <- rbind(PCA1, hor)
    PCA2 <- rbind(PCA2, ver)
    PCA3 <- rbind(PCA3, dia)
  }
}
  
coef.data

for (subjectID in subjectpool){
  
  sub.data <- coef.data[coef.data$ID == subjectID,]
  sPCA1 <- PCA1[PCA1[,1] == subjectID,]
  sPCA2 <- PCA2[PCA1[,1] == subjectID,]
  sPCA3 <- PCA3[PCA1[,1] == subjectID,]
  
  sub.data$PCA <- ifelse(sub.data$path2 == "hor", sPCA1[,2], 
                        ifelse(sub.data$path2 == "ver", sPCA2[,2],
                               ifelse(sub.data$path2 == "dia", sPCA3[,2], "NA")))
  
  if (subjectID == 1){
    cos.data <- sub.data
  } else {
    cos.data <- rbind(cos.data, sub.data)
  }
}

### get cosine distance
cos.d.data <- rep(NA, 4)

for (ID in subjectpool){
  for (run in 1:6){
    for (path in 1:12){
      
      sub.data <- cos.data[cos.data$run == run & cos.data$path == path & cos.data$ID == ID,]
      vec1 <- sub.data$coef
      vec2 <- as.numeric(sub.data$PCA)
      cos <- cosine(vec1, vec2)
      
      cos.d.data <- rbind(cos.d.data, c(ID, run, path, cos))
      
    }}
}

cos.d.data <- data.frame(cos.d.data[-1,])
names(cos.d.data) <- c("ID", "run", "path", "cosine dist")

### plot cosine distance
plot.cos.data <- tapply(cos.d.data$`cosine dist`, list(cos.d.data$run, cos.d.data$path), mean)
plot.cos.data.ID <- tapply(cos.d.data$`cosine dist`, list(cos.d.data$run, cos.d.data$path, cos.d.data$ID), mean)


cos.sd <- matrix(NA, 6, 12)
for (i in 1:6){
  for (j in 1:12){
    cos.sd[i,j] <- (sd(plot.cos.data.ID[i,j,]))
  }}
cos.sd

plot.cos.data <- reshape2::melt(plot.cos.data)
cos.sd <- reshape2::melt(cos.sd)
plot.cos.data$sd <- cos.sd$value


for (j in 1:12){
  p <- ggplot(plot.cos.data[plot.cos.data$Var2 == j,], aes(x=Var1, y=value)) +
    geom_point() +
    geom_line() +
    scale_y_continuous("Cosine Distance")+#, limit = c(0.85, 1.03), breaks=seq(0.85, 1.05, 0.01)) +
    scale_x_continuous(paste0("Path ",j),breaks=c(1:6)) +
    theme(legend.position="none") +
    #coord_cartesian(ylim = c(0.87, 1)) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05))
  
  assign(paste0("p",j), p)
}

png(paste0("/Users/jisu/Documents/Glove/cca_cos.png"), width = 4000, height = 2500, res=220)
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, cols=4)
dev.off()

plot.cos.data$Path.label <- path2[plot.cos.data$Var2]
colnames(plot.cos.data) = c("Run", "Path", "Value", "sd", "Path.label")
plot.cos.data$Path2 <- as.character(plot.cos.data$Path)


p <- ggplot(plot.cos.data, aes(x=Run, y=Value, group=Path2)) +
  geom_point(aes(color=Path2)) +
  geom_line(aes(color=Path2)) +
  scale_y_continuous("Cosine Distance")+#, limit = c(0.89, 1), breaks=seq(0.90, 1.00, 0.01)) +
  scale_x_continuous(paste0("Run"),breaks=c(1:6)) +
  scale_color_discrete(name = "Path",
                       breaks = c(1:12),
                       label = c("1","2","3","4","5","6","7","8","9","10","11","12"))


png(paste0("/Users/jisu/Documents/Glove/cca_cos_group.png"), width = 2000, height = 2000, res=220)
p
dev.off()


names(stat.data) <- c(names(stat.data[1:3]), paste0("s",1:14), names(stat.data[18:19]), "ID")
names(stat.data)
stat.data <- data.frame(stat.data)
### plotting d-sensor and d-dist
### wtf 

for (j in 1:14){
  p <- ggplot(data = stat.data) +
    geom_point(aes(x = stat.data[,j+3], y=d.dist), size=0.01) +
    scale_x_continuous(paste0("Sensor", j)) +
    scale_y_continuous("Differences of Distance")
  
    #assign(paste0("p",j), p)
    png(paste0("/Users/jisu/Documents/Glove/cca_sensor_dist_sensor",j,".png"), width = 1000, height = 1000, res=220)
    p
    dev.off()
}

#png(paste0("/Users/jisu/Documents/Glove/cca_sensor_dist.png"), width = 4000, height = 1200, res=220)
#multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, cols=7)
#dev.off()




  sub.data <- stat.data[(stat.data$ID == i & stat.data$path == path[j] & stat.data$run == k),]
  result.lasso <- Lasso(sub.data[,c(4:17)], sub.data$d.dist, fix.lambda = FALSE)
  lm(d.dist ~ s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14, sub.data)
  

  p1 <- ggplot(data = sub.data)+
    geom_point( aes(x=s1, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[1], intercept = result.lasso$beta0, color="red")
  
  p2 <- ggplot(data = sub.data)+
    geom_point( aes(x=s2, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[2], intercept = result.lasso$beta0, color="red")
  
  p3 <- ggplot(data = sub.data)+
    geom_point( aes(x=s3, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[3], intercept = result.lasso$beta0, color="red")
  
  p4 <- ggplot(data = sub.data)+
    geom_point( aes(x=s4, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[4], intercept = result.lasso$beta0, color="red")
  
  p5 <- ggplot(data = sub.data)+
    geom_point( aes(x=s5, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[5], intercept = result.lasso$beta0, color="red")
  
  p6 <- ggplot(data = sub.data)+
    geom_point( aes(x=s6, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[6], intercept = result.lasso$beta0, color="red")
  
  p7 <- ggplot(data = sub.data)+
    geom_point( aes(x=s7, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[7], intercept = result.lasso$beta0, color="red")
  
  p8 <- ggplot(data = sub.data)+
    geom_point( aes(x=s8, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[8], intercept = result.lasso$beta0, color="red")
  
  p9 <- ggplot(data = sub.data)+
    geom_point( aes(x=s9, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[9], intercept = result.lasso$beta0, color="red")
  
  p10 <- ggplot(data = sub.data)+
    geom_point( aes(x=s10, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[10], intercept = result.lasso$beta0, color="red")
  
  p11 <- ggplot(data = sub.data)+
    geom_point( aes(x=s11, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[11], intercept = result.lasso$beta0, color="red")
  
  p12 <- ggplot(data = sub.data)+
    geom_point( aes(x=s12, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[12], intercept = result.lasso$beta0, color="red")
  
  p13 <- ggplot(data = sub.data)+
    geom_point( aes(x=s13, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[13], intercept = result.lasso$beta0, color="red")
  
  p14 <- ggplot(data = sub.data)+
    geom_point( aes(x=s14, y=d.dist), size=0.01)+
    scale_x_continuous(limits=c(-1, 1))+
    geom_abline(slope = result.lasso$beta[14], intercept = result.lasso$beta0, color="red")
    
  
  png(paste0("/Users/jisu/Documents/Glove/ID26path1run6.png"), width = 4000, height = 1000, res=220)
  multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, cols=7)
  dev.off()
  
  
  p <- ggplot(plot.cos.data[plot.cos.data$Var2 == j,], aes(x=Var1, y=value)) +
    geom_point() +
    geom_line() +
    scale_y_continuous("Cosine Distance")+#, limit = c(0.85, 1.03), breaks=seq(0.85, 1.05, 0.01)) +
    scale_x_continuous(paste0("Path ",j),breaks=c(1:6)) +
    theme(legend.position="none") +
    #coord_cartesian(ylim = c(0.87, 1)) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05))
  
  assign(paste0("p",j), p)


png(paste0("/Users/jisu/Documents/Glove/cca_cos.png"), width = 4000, height = 2500, res=220)
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, cols=4)
dev.off()


# dsi - dr의 scatter plot 참가자별로 / path별로 
# dsensor - dr의 기울기

head(stat.data)
dim(stat.data)
write.csv(stat.data, "/Users/jisu/Documents/Glove/sensor_dist_data.csv")
  