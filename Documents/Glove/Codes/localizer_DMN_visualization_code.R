################
### 190906 DMN visualization ###
############


library(ggplot2)
library(reshape2)
library(FSA)

lda <- read.csv("/Users/jisu/Documents//Glove/20190906_DMN_lda_localizer.csv")
#svc <- read.csv("/Users/jisu/Documents//Glove/20190906_DMN_svc_all.csv")


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
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


data <- lda
dim(data)
data$period <- ifelse(substr(data$subject, 1, 2) == "GA", "early", 
                      ifelse(substr(data$subject, 1, 2) == "GB", "late", "NA"))

subdata = data[data$target == "pos",]    

# 1. Early - Late 2. Prac - Unprac 3.roi
visudata.m = tapply(subdata$accu, list(subdata$practiced, subdata$roi, subdata$period), mean)
visudata.se = tapply(subdata$accu, list(subdata$practiced, subdata$roi, subdata$period), se)

visu.m <- melt(visudata.m)
visu.se <- melt(visudata.se) 
visudata <- cbind(visu.m, visu.se$value)
names(visudata) <- c("mapping", "roi", "period", "mean", "se")
visudata <- data.frame(visudata)


avdata.m <- tapply(visudata$mean, list(visudata$mapping, visudata$period), mean)
avdata.se <- tapply(visudata$mean, list(visudata$mapping, visudata$period), se)

av.m <- melt(avdata.m)
av.se <- melt(avdata.se)
avdata <- data.frame(cbind(c("practiced", "unpracticed", "practiced", "unpracticed"), 10, 
c("early", "early", "late", "late"), 
av.m$value, av.se$value))
names(avdata) <- c("mapping", "roi", "period", "mean", "se")

visudata <- rbind(visudata, avdata)
visudata$mean <- as.numeric(visudata$mean)
visudata$se <- as.numeric(visudata$se)


ROInames = c(
   'SPL',
   'MFG',
   'MTG',
   'IFG',
   'SPL',
   'MFG',
   'MFG2',
   'IFG',
   'MTG',
   'Average'
)

for (i in 1:length(ROInames)){
  p <- ggplot(data = visudata[visudata$roi == i,], aes(x = period, y = mean)) +
    scale_y_continuous(name="mean_acc", limits = c(0.26, 0.42)) + 
    scale_x_discrete(name = paste0(ROInames[i]," (LDA)")) +
    geom_point(aes(col = mapping)) +
    geom_line(aes(group=mapping, col=mapping)) +
    geom_errorbar(aes(x = period, ymin = mean - se, ymax = mean + se, col=mapping, width = .2)) +
    theme(legend.title = element_text(size=8, color = "salmon", face="bold"),
          legend.justification=c(1,0), 
          legend.position=c(0.3, 0.05),  
          legend.background = element_blank(),
          legend.key = element_blank())
  
  assign(paste0("p", i), p)
}


png(paste0("/Users/jisu/Documents/Glove/20190911_localizer_DMN_LDA.png"), width = 1000, height = 9000, res=220)
multiplot(
  p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
)
dev.off()




#data <- svc
#dim(data)
#data$period <- ifelse(substr(data$subject, 1, 2) == "GA", "early", 
#                      ifelse(substr(data$subject, 1, 2) == "GB", "late", "NA"))
#                     
#subdata = data[data$target == "pos",]     
## 1. Early - Late 2. Prac - Unprac 3.roi
#visudata.m = tapply(subdata$accu, list(subdata$practiced, subdata$roi, subdata$period), mean)
#visudata.se = tapply(subdata$accu, list(subdata$practiced, subdata$roi, subdata$period), se)
#
#visu.m <- melt(visudata.m)
#visu.se <- melt(visudata.se) 
#visudata <- cbind(visu.m, visu.se$value)
#visudata <- data.frame(visudata)
#names(visudata) <- c("mapping", "roi", "period", "mean", "se")
#
#ROInames = c("Core", "DMsub", "MTLsub")
#
#for (i in 1:3){
#p0 <- ggplot(data = visudata[visudata$roi == 1,], aes(x = period, y = mean)) +
#        scale_y_continuous(name="mean_acc") + 
#        scale_x_discrete(name = paste0(ROInames[i]," (SVC)")) +
#        geom_point(aes(col = mapping)) +
#        geom_line(aes(group=mapping, col=mapping)) +
#        geom_errorbar(aes(x = period, ymin = mean - se, ymax = mean + se, col=mapping, width = .2)) +
#           theme(legend.title = element_text(size=8, color = "salmon", face="bold"),
#           legend.justification=c(1,0), 
#           legend.position=c(1, 0.05),  
#           legend.background = element_blank(),
#           legend.key = element_blank())
#
#assign(paste0("p0", i), p0)
#}


#png(paste0("/Users/jisulee/Documents/CLMN LAB/Glove Project/20190906_DMN_LDA_SVC.png"), width = 3000, height = 3000, res=220)
#multiplot(p1, p2, p3, p01, p02, p03)
#dev.off()


### Corrlation w/ learning amount
learning.amount <- c(0.417418981481482,
                     0.360081018518519,
                     0.346666666666667,
                     0.412719907407407,
                     0.212870370370370,
                     0.475243055555556,
                     0.345543981481481,
                     0.329375000000000,
                     0.545196759259259,
                     0.419328703703704,
                     0.219363425925926,
                     0.317615740740741,
                     0.264768518518519,
                     0.539016203703704,
                     0.164085648148148,
                     0.526053240740741,
                     0.292233796296296,
                     0.297141203703704,
                     0.474976851851852,
                     0.150625000000000,
                     0.221759259259259,
                     0.165706018518519,
                     0.687083333333333,
                     0.402210648148148,
                     0.373993055555556,
                     0.411678240740741,
                     0.529895833333333,
                     0.473599537037037,
                     0.413530092592593,
                     0.442511574074074)

subdata <- data[data$practiced == "practiced" & data$target == "pos",]
GAsub <- subdata[substr(subdata$subject, 1, 2) == "GA",]
GBsub <- subdata[substr(subdata$subject, 1, 2) == "GB",]

diffsub <- cbind(GAsub, GBsub$accu, "accu.diff" = (GBsub$accu - GAsub$accu), "learn" = rep(learning.amount, each=length(ROInames)))
#diffsub.Core <- diffsub[diffsub$roi == 1,]
#diffsub.DMsub <- diffsub[diffsub$roi == 2,]
#diffsub.MTLsub <- diffsub[diffsub$roi == 3,]


for (i in 1:length(ROInames)){
  cor.value = cor(diffsub$learn[diffsub$roi == i], diffsub$accu.diff[diffsub$roi == i])
  grob = grobTree(textGrob(paste0("r = ",round(cor.value, 2)), x=0.05,  y=0.90, hjust=0), gp=gpar(fontsize=20, fonface="Bold"))
  
  c <- ggplot(data = diffsub[diffsub$roi == i,], aes(x=accu.diff, y=learn)) +
    scale_y_continuous(name = "reward difference (GB-GA)") + 
    scale_x_continuous(name = paste0(ROInames[i]," (prac_late-early)")) +
    #geom_text(x=0.05, y=0.7, label = paste0("r = ",round(cor.value, 2))) +
    geom_point() +    # Use hollow circles
    geom_smooth(method=lm) + # Add linear regression line 
    annotation_custom(grob)
  # Plot
  
  assign(paste0("c", i), c)
}



#png(paste0("/Users/jisu/Documents//Glove/20190906_DMN_LDA_corr_clusters.png"), width = 2000, height = 2000, res=220)
#multiplot(p1, p6, p16,
#          c1, c6, c16)
#dev.off()




for (i in 1:54){
  print(paste0("c",i))
}


