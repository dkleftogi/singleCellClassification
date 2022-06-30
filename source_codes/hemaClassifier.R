#BEGIN COPYRIGHT NOTICE

#hemaClassifier -- (c) 2022 Dimitrios Kleftogiannis -- UiB and CCBIO
#Copyright 2022 University of Bergen (UiB), Norway.
#This Program is free software licensed under the MIT License.
#You may only use the source code in this repository in compliance with the license provided in this repository. 
#For more details, please refer to the file named "LICENSE.md".
#This Program is distributed as a service to the research community and is experimental in nature and may have hazardous properties. 
#The Program is distributed WITHOUT ANY WARRANTY, express or implied. In particular all warranties as to SATISFACTORY QUALITY or FITNESS FOR A PARTICULAR PURPOSE are excluded. 
#Published reports of research using this code (or a modified version) should cite the relevant article of this code
#Comments and bug reports are welcome.
#Email to dimitrios.kleftogiannis@uib.no
#I would also appreciate hearing about how you used this code, improvements that you have made to it.
#You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 

#END COPYRIGHT NOTICE

#UTILITY
#This code visualises the annotated healthy data and trains a hematopoietic classifier
#It can be used to make predictions in leukemia data and estimate the cellular abundance 
#The data used in this tutorial are only for illustration purposes.

#RUNNING
#We provide a step-by-step execution example, this is useful for illustration and to reproduce parts of the analysis. 
#Remember that this is a minimal example. To obtain the input dataplease contact Dimitrios Kleftogiannis.

#DEPENDENCIES
#To run the tutorial without problem, is required to install and load all packages listed in implementGatingStrategy example

#load a "toy" dataset
load('HD_annotated.Rdata')

#merge individual cell types
HD_annotated <- rbind(Bcells,NK,TcellsCyto,TcellsHelper,
                      MPP,CMP,Myeloid_a1,Myeloid_a2)

#estimate summary abundances
cellularAbundance <- as.data.frame(table(HD_annotated$Type))
cellularAbundance$Ratio <- 100*cellularAbundance$Freq/nrow(HD_annotated)
colnames(cellularAbundance)[1]<- 'CellType'

#visualise using t-SNE the myeloid and lymphoid compartments

#since MPP and CMP contain very few cellsbwe subsample Granulocytes and Monocytes and project it all together 
#10% of granulocytes
a <- sample(nrow(Myeloid_a1))
a_select <- a[1:round(0.1*nrow(Myeloid_a1))]
myel1 <- Myeloid_a1[a_select,]
#10% of monocytes
a_select <- a[1:round(0.1*nrow(Myeloid_a2))]
myel2 <- Myeloid_a2[a_select,]
myeloid_sampled <- rbind(myel1,myel2,MPP,CMP)
myeloidLabels <- myeloid_sampled$Type

#select specific specific cell identity markers
selectedMarkers <- c('CD45','CD33','CD38','CD34','CD123','CD20','CD3','CD8a','CD4','CD7','CD56','CD66b','CD64','HLA-DR','CD14','CD16')

myeloid_sampled <- myeloid_sampled[,c(selectedMarkers)]
#perform tsne for myeloid; run Rtsne with the following setting
d <- myeloid_sampled
duplicates <- duplicated(d) | duplicated(d, fromLast = TRUE) 
tsne_out<-Rtsne(d[!duplicates,1:ncol(d)],max_iter=3000,perplexity=30,seed=1234,num_threads=4)
tsne_data <- myeloid_sampled[!duplicates,]
tsne_data$TSNE1 <- tsne_out$Y[,1]
tsne_data$TSNE2 <- tsne_out$Y[,2]
tsne_data$CellType <- myeloidLabels[!duplicates]
tsne_data_myeloid <- tsne_data

#visualise  tsne
tsne_myeloid <- ggplot(tsne_data_myeloid,  aes(x = TSNE1, y = TSNE2, color = CellType )) +
  geom_point(size = 0.58,alpha=0.6) +
  guides(size=FALSE)+
  scale_color_brewer(palette = "Set1",name='Cell type')+
  #scale_color_manual(values = rep(1:8, len = 8),name='Cell type')+
  theme_minimal()+
  theme(axis.text.y = element_text( size = 14 ),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 14),
        legend.position = "bottom",
        legend.text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),aspect.ratio = 1)+
  guides(color=guide_legend(nrow=2,override.aes = list(size=2.8)))
#save the plot
myfile<-paste('plot8.pdf')
pdf(myfile,onefile = TRUE)
print(tsne_myeloid)
dev.off()

#for lymphoid cell types we visualise 20% of all data points to have ~10K cells similar to myeloid
lymphoid_sampled <- rbind(Bcells,NK,TcellsCyto,TcellsHelper)
a <- sample(nrow(lymphoid_sampled))
a_select <- a[1:round(0.2*nrow(lymphoid_sampled))]
lymphoid_sampled <- lymphoid_sampled[a_select,]
lymphoidLabels <- lymphoid_sampled$Type

#perform tsne for myeloid; run Rtsne with the following setting
d <- lymphoid_sampled
duplicates <- duplicated(d) | duplicated(d, fromLast = TRUE) 
tsne_out<-Rtsne(d[!duplicates,1:ncol(d)],max_iter=3000,perplexity=30,seed=1234,num_threads=4)
tsne_data <- lymphoid_sampled[!duplicates,]
tsne_data$TSNE1 <- tsne_out$Y[,1]
tsne_data$TSNE2 <- tsne_out$Y[,2]
tsne_data$CellType <- lymphoidLabels[!duplicates]
tsne_data_lymphoid <- tsne_data

#visualise tsne
tsne_lymphoid <- ggplot(tsne_data_lymphoid,  aes(x = TSNE1, y = TSNE2, color = CellType )) +
  geom_point(size = 0.58,alpha=0.6) +
  guides(size=FALSE)+
  scale_color_brewer(palette = "Set2",name='Cell type')+
  #scale_color_manual(values = rep(1:8, len = 8),name='Cell type')+
  theme_minimal()+
  theme(axis.text.y = element_text( size = 14 ),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 14),
        legend.position = "bottom",
        legend.text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),aspect.ratio = 1)+
  guides(color=guide_legend(nrow=2,override.aes = list(size=2.8)))
#save the plot
myfile<-paste('plot9.pdf')
pdf(myfile,onefile = TRUE)
print(tsne_myeloid)
dev.off()

#perform PCA using all annotated healthy data

#perform PCA using all data - d contains only markers
d <- HD_annotated[,selectedMarkers]
#save separately the labels
allLabels <- HD_annotated$Type
#run pca --> inherited from stats package
pca <- prcomp(d) 
#check the important components
summary(pca)

#use ggbiplot package to make a "ready-made" PCA plot
#load the package, install it first if not ready
library(ggbiplot)
o <- ggbiplot(pca, groups=HD_annotated$Type,ellipse=TRUE,
              circle = FALSE,pc.biplot = TRUE,varname.size=1.2)+
  theme_classic()+
  scale_colour_brewer(palette = "Set2",name='Cell types')+
  theme(legend.position = "bottom")

#save the plot
myfile<-paste('plot10.tiff')
tiff(myfile,width=6, height=6, units="in", res=360)
print(o)
dev.off()

#in the previous plot we use circles to visualise the cell populations -> clusters 
#we can estimate the centroid of each cluster using euclidean distance
cellTypes <- unique(HD_annotated$Type)
centroids <- data.frame()
for(idx in cellTypes){
  
  a1 <- which(HD_annotated$Type==idx)
  #estimate the "average" coordinate
  dt1 <- pca$x[a1,1]
  ct1 <- mean(dt1)
  dt2 <- pca$x[a1,2]
  ct2 <- mean(dt2)
  #estimate the average distance of all points in the cell type from the centroid
  myDist <- data.frame(X=dt1,CT1=rep(ct1,len=length(dt1)),Y=dt2,CT2=rep(ct2,len=length(dt1)))
  myD <- matrix(-1,nrow = nrow(myDist),ncol = 1)
  for(i in 1:nrow(myDist)){
    myD[i] <- sqrt( (myDist[i,1]-myDist[i,2])^2 + (myDist[i,3]-myDist[i,4])^2)
  }
  df <- data.frame(PC1=ct1,PC2=ct2,Label=idx, Radius=median(myD))
  centroids <- rbind(centroids,df)
}

#to visualise the PCA we sample 30K cells randomly
#first shuffle the annotated healthy cells
a <- sample(nrow(HD_annotated))
#select them
a_select <- a[1:30000]
plot_data <- HD_annotated[a_select,]
d <- plot_data[,selectedMarkers]
myLabels <- plot_data$Type
#generate a new data frame with the selected cells in the PCA space and the expression of  specific markers
#modify this code to visualise specific markers
plot_data <- data.frame(PC1=pca$x[a_select,1],PC2=pca$x[a_select,2],Label=myLabels, 
                        CD45=HD_annotated[a_select,'CD45'],
                        CD34=HD_annotated[a_select,'CD34'],
                        CD38=HD_annotated[a_select,'CD38'],
                        CD33=HD_annotated[a_select,'CD33'],
                        `HLA-DR`=HD_annotated[a_select,'HLA-DR'],
                        CD66b=HD_annotated[a_select,'CD66b'])

#run min-max normalisation just to rescale in the same range the values - it is just for visualisation not for "real" analyses
plot_data$CD45 <- (plot_data$CD45-min(plot_data$CD45))/(max(plot_data$CD45)-min(plot_data$CD45))
plot_data$CD34 <- (plot_data$CD34-min(plot_data$CD34))/(max(plot_data$CD34)-min(plot_data$CD34))
plot_data$CD38 <- (plot_data$CD38-min(plot_data$CD38))/(max(plot_data$CD38)-min(plot_data$CD38))
plot_data$CD33 <- (plot_data$CD33-min(plot_data$CD33))/(max(plot_data$CD33)-min(plot_data$CD33))
plot_data$HLA.DR <- (plot_data$HLA.DR-min(plot_data$HLA.DR))/(max(plot_data$HLA.DR)-min(plot_data$HLA.DR))
plot_data$CD66b <- (plot_data$CD66b-min(plot_data$CD66b))/(max(plot_data$CD66b)-min(plot_data$CD66b))

#generate the plot with centroids
myPCA <- ggplot() +
  geom_point(data=plot_data,  aes(x = PC1, y = PC2, color = CD45),size = 0.2,alpha=0.6)+
  scale_color_gradientn("CD45",colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50),
                        limits=c(0,1),breaks = c(0,0.5,1),values=c('0','0.5','1'))+
  geom_point(data=centroids,  aes(x = PC1, y = PC2, shape = Label,group=1),size = 2.8)+
  geom_point(data=centroids,  aes(x=PC1,y=PC2,size=Radius),alpha = 0.2)+
  guides(size=FALSE)+
  #geom_circle(data=a2,  aes(x=PC1,y=PC2,radius=Radius))+
  #geom_circle(data=a3,  aes(x=PC1,y=PC2,radius=Radius))
  #geom_polygon(data = hulls, alpha = 0.5)+
  #stat_ellipse(aes(x = PC1, y = PC2, group = Label ),type = "norm")+
  theme_minimal()+
  ylab('PC2')+
  xlab('PC1')+
  scale_shape_manual(values = rep(1:8, len = 8),name='Cell type')+
  theme(axis.text.y = element_text( size = 14 ),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 14),
        legend.position = "bottom",
        legend.text=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),aspect.ratio = 1)+
  guides(shape=guide_legend(ncol=3))

myfile<-paste('plot11.pdf')
pdf(myfile,onefile = TRUE)
print(myPCA)
dev.off()

#redo the PCA plot but color according to marker expression
p1 <- ggplot() +
  geom_point(data=plot_data,  aes(x = PC1, y = PC2, color = CD38),size = 0.2,alpha=0.6)+
  scale_color_gradientn("CD38",colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50),
                        limits=c(0,1),breaks = c(0,0.5,1),values=c('0','0.5','1'))+
  geom_point(data=centroids,  aes(x = PC1, y = PC2, shape = Label,group=1),size = 2.8)+
  geom_point(data=centroids,  aes(x=PC1,y=PC2,size=Radius),alpha = 0.2)+
  guides(size=FALSE)+
  #geom_circle(data=a2,  aes(x=PC1,y=PC2,radius=Radius))+
  #geom_circle(data=a3,  aes(x=PC1,y=PC2,radius=Radius))
  #geom_polygon(data = hulls, alpha = 0.5)+
  #stat_ellipse(aes(x = PC1, y = PC2, group = Label ),type = "norm")+
  theme_minimal()+
  scale_shape_manual(values = rep(1:8, len = 8),name='Cell type')+
  theme(axis.text.y = element_text( size = 14 ),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 14),
        legend.position = "bottom",
        legend.text=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),aspect.ratio = 1)+
  guides(shape=guide_legend(ncol=3))

myfile<-paste('plot12.pdf')
pdf(myfile,onefile = TRUE)
print(p1)
dev.off()

#in the same way we can visualise any marker available in our antibody panel
#repeat accordindly or wrap around a for-loop
#here just show for CD66b the typical granulocytic marker
p2 <- ggplot() +
  geom_point(data=plot_data,  aes(x = PC1, y = PC2, color = CD66b),size = 0.2,alpha=0.6)+
  scale_color_gradientn("CD66b",colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50),
                        limits=c(0,1),breaks = c(0,0.5,1),values=c('0','0.5','1'))+
  geom_point(data=centroids,  aes(x = PC1, y = PC2, shape = Label,group=1),size = 2.8)+
  geom_point(data=centroids,  aes(x=PC1,y=PC2,size=Radius),alpha = 0.2)+
  guides(size=FALSE)+
  #geom_circle(data=a2,  aes(x=PC1,y=PC2,radius=Radius))+
  #geom_circle(data=a3,  aes(x=PC1,y=PC2,radius=Radius))
  #geom_polygon(data = hulls, alpha = 0.5)+
  #stat_ellipse(aes(x = PC1, y = PC2, group = Label ),type = "norm")+
  theme_minimal()+
  scale_shape_manual(values = rep(1:8, len = 8),name='Cell type')+
  theme(axis.text.y = element_text( size = 14 ),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 14),
        legend.position = "bottom",
        legend.text=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),aspect.ratio = 1)+
  guides(shape=guide_legend(ncol=3))

myfile<-paste('plot13.pdf')
pdf(myfile,onefile = TRUE)
print(p2)
dev.off()

#for sanity reasons we can visualise the expression profiles of selected markers using density plots
#we need to make some numeric tranformations of the data as follows
data <- HD_annotated[,1:(ncol(HD_annotated)-1)]
data <- as.matrix(data)
rng <- colQuantiles(data, probs = c(0.01, 0.99))
data01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01[data01 < 0] <- 0
data01[data01 > 1] <- 1
#more transformations and name changing due to R 
test <- data.frame(data01)

selectedMarkers <- c('CD45','CD33','CD38','CD34','CD123','CD20','CD3','CD8a','CD4','CD7','CD56','CD66b','CD64','HLA.DR','CD14','CD16')
test <- test[,selectedMarkers]
test$Class <- HD_annotated$Type

test$Class <- factor(test$Class,levels=c('B-cells','NK',
                                         'T-cytotoxic','T-helper',
                                         'Granulocytes','Monocytes',
                                         'CMP','MPP'))
test <- melt(test)
HV_profiles <- ggplot(test, aes(x = value, y = Class, fill = variable)) + 
  geom_density_ridges(scale=1)+
  facet_wrap(~ variable, scales = "free", nrow = 4)+
  theme_ridges() + 
  ylab('Density')+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1,size = 7),
        axis.text.y = element_text(hjust=1,size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),aspect.ratio = 1.58)+
  scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))

myfile<-paste('plot14.pdf')
pdf(myfile,onefile = TRUE)
print(HV_profiles)
dev.off()

########################################################################################


#train a simple LDA classifier and assess the performance using 2-fold cross validation
#install and use MASS package for that
library(MASS)
#how many cells we have, split into half
myCells <- nrow(HD_annotated)
myprop <- round(0.5*myCells)
#shuffle the cells
a <- sample(nrow(HD_annotated))
#take the first half of it for training, the rest remaining for testing
a_select <- a[1:myprop]
training <- HD_annotated[a_select,]
a_select <-a[(myprop+1):length(a)]
testing <- HD_annotated[a_select,]

#develop a model using the high quality data of this iteration

#FOLD 1: train using the training data
model <- lda(Type~.,data=na.omit(training))
myTest1 <- testing[,1:(ncol(testing)-1)]
myLabels <- testing[,ncol(testing)]
a<-predict(model,myTest1)
#store the posterior probabilities
prob_matrix<-matrix(0L, nrow = nrow(myTest1), ncol =1)
prob_matrix<-as.matrix(apply(a$posterior,1,max))
DR<-data.frame(myTest1,Previous=as.factor(myLabels),New=a$class,Prob=prob_matrix)

#filter based on probability - here we used 0.85 but this can be changed accordingly. In the paper this optimisation was skipped
#estimate the number of "rejected cells" - it is a good idea to investigate the rejected cells with clustering
rejectedCells <- which(DR$Prob<0.85)
#filter out the rejected cells
DR <- DR[-rejectedCells,]
#generate miss classification plots
clusters <- unique(DR$Previous)
performance_fold1 <- data.frame()
#store missclassification plots for fold1
pf1 <- list()
myC <- 1
for(idx in clusters){
  
  a <- which(DR$New==idx)
  tmp <- DR[a,]
  Nj <- nrow(tmp)
  predicted_res <- table(tmp$Previous)
  plot_data <- data.frame(ClusterSize=nrow(tmp),t(predicted_res))
  colnames(plot_data)[2] <- idx
  colnames(plot_data)[3] <- 'PredictedCellType'
  plot_data$Freq <- 100*plot_data$Freq/nrow(tmp)
  str <- paste('Gated cell type: ',idx,sep='')
  #one subplot per cell type -- plots fold 1
  pf1[[myC]] <- ggplot(plot_data, aes(x=PredictedCellType,y=Freq,size=Freq,color=PredictedCellType))+
    geom_point(alpha=0.5)+
    geom_linerange(aes(x=PredictedCellType, ymin=0, ymax=Freq),size=1)+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle=45, hjust=1,size = 10),
          axis.text.y = element_text(hjust=1,size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=10),aspect.ratio = 0.5)+ ylim(0,100)+
    ggtitle(str)+scale_colour_brewer(palette = "Set2")
  myName <- which.max(table(tmp$Previous))
  s1 <- which(names(table(tmp$Previous))==idx)
  Nij <- table(tmp$Previous)[s1]
  predictedName <- idx
  a <- which(DR$Previous==predictedName)
  tmp <- DR[a,]
  Ni <- nrow(tmp)
  df <- data.frame(CellType=paste(idx,sep=''),PredictedCellType=predictedName,
                   Recall=Nij/Ni,
                   Precision=Nij/Nj,
                   F1=2*(Nij/Ni)*(Nij/Nj)/((Nij/Ni)+(Nij/Nj)))
  #store precision and recall per cell type
  performance_fold1 <- rbind(performance_fold1,df)
  myC <- myC + 1
}

#FOLD 2: train using the previous testing data
model <- lda(Type~.,data=na.omit(testing))
myTest2 <- training[,1:(ncol(training)-1)]
myLabels <- training[,ncol(training)]
a<-predict(model,myTest2)
#store the posterior probabilities
prob_matrix<-matrix(0L, nrow = nrow(myTest2), ncol =1)
prob_matrix<-as.matrix(apply(a$posterior,1,max))
DR<-data.frame(myTest2,Previous=as.factor(myLabels),New=a$class,Prob=prob_matrix)
#filter based on probability - we used 0.85 but this can be changed accordingly
#estimate the number of "rejected cells"
rejectedCells <- which(DR$Prob<0.85)
#filter out the rejected cells
DR <- DR[-rejectedCells,]
#generate miss classification plots
clusters <- unique(DR$Previous)
performance_fold2 <- data.frame()
#store missclassification plots for fold 2
pf2 <- list()
myC <- 1
for(idx in clusters){
  
  a <- which(DR$New==idx)
  tmp <- DR[a,]
  Nj <- nrow(tmp)
  predicted_res <- table(tmp$Previous)
  plot_data <- data.frame(ClusterSize=nrow(tmp),t(predicted_res))
  colnames(plot_data)[2] <- idx
  colnames(plot_data)[3] <- 'PredictedCellType'
  plot_data$Freq <- 100*plot_data$Freq/nrow(tmp)
  str <- paste('Gated cell type: ',idx,sep='')
  #one subplot per cell type -- plots fold 1
  pf2[[myC]] <- ggplot(plot_data, aes(x=PredictedCellType,y=Freq,size=Freq,color=PredictedCellType))+
    geom_point(alpha=0.5)+
    geom_linerange(aes(x=PredictedCellType, ymin=0, ymax=Freq),size=1)+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle=45, hjust=1,size = 10),
          axis.text.y = element_text(hjust=1,size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=10),aspect.ratio = 0.5)+ ylim(0,100)+
    ggtitle(str)+scale_colour_brewer(palette = "Set2")
  myName <- which.max(table(tmp$Previous))
  s1 <- which(names(table(tmp$Previous))==idx)
  Nij <- table(tmp$Previous)[s1]
  predictedName <- idx
  a <- which(DR$Previous==predictedName)
  tmp <- DR[a,]
  Ni <- nrow(tmp)
  df <- data.frame(CellType=paste(idx,sep=''),PredictedCellType=predictedName,
                   Recall=Nij/Ni,
                   Precision=Nij/Nj,
                   F1=2*(Nij/Ni)*(Nij/Nj)/((Nij/Ni)+(Nij/Nj)))
  #store precision and recall per cell type
  performance_fold2 <- rbind(performance_fold2,df)
  myC <- myC + 1
}

#save cell misclassification plots for fold 1 and fold 2 
#cell type classification performance can be found in performance_fold1 and performance_fold2 data frames
combined_plot_f1 <- ggarrange(pf1[[1]],pf1[[2]],pf1[[3]],pf1[[4]],pf1[[5]],pf1[[6]],pf1[[7]],pf1[[8]],
                         nrow=3,ncol=3)
myfile<-paste('plot15_f1.pdf')
pdf(myfile,onefile = TRUE)
print(combined_plot_f1)
dev.off()

combined_plot_f2 <- ggarrange(pf2[[1]],pf2[[2]],pf2[[3]],pf2[[4]],pf2[[5]],pf2[[6]],pf2[[7]],pf2[[8]],
                              nrow=3,ncol=3)
myfile<-paste('plot15_f2.pdf')
pdf(myfile,onefile = TRUE)
print(combined_plot_f2)
dev.off()


#after tuning the hematopoietic classifier we can generate one model using all available annotated cells

#first choose the markers for further analysis - surface for cell type annotation and signalling too
studiedMarkers <- c('CD45','CD33','CD38','CD34','CD123','CD20','CD3','CD8a','CD4','CD7','CD56','CD66b','CD64','HLA-DR','CD14','CD16',
                     'pNFkB','pErk','pSTAT1','pSTAT3','pSTAT5','pP38','pHistone3','pRB','Caspase3_cleaved','pAxl','pAkt','pS6')
#selected surface markers for cell type annotation
selectedMarkers <- c('CD45','CD33','CD38','CD34','CD123','CD20','CD3','CD8a','CD4','CD7','CD56','CD66b','CD64','HLA-DR','CD14','CD16')
#modelling using selected markers from annotated HD data frame
tmp <- HD_annotated[,c('Type',selectedMarkers)]
#this model can be used for transfer-learning using data from leukemia patients
model <- lda(Type~.,data=na.omit(tmp))


#############################################################
#train a simple ΚΝΝ classifier and assess the performance using 2-fold cross validation

#how many cells we have, split into half
myCells <- nrow(HD_annotated)
myprop <- round(0.5*myCells)
#shuffle the cells
a <- sample(nrow(HD_annotated))
#take the first half of it for training, the rest remaining for testing
a_select <- a[1:myprop]
training <- HD_annotated[a_select,]
a_select <-a[(myprop+1):length(a)]
testing <- HD_annotated[a_select,]

#develop a model using the high quality data of this iteration

#projection based on k-nn
data <- training[,1:16]
dataLabels <- training$Type
dataLabels <- as.character(dataLabels)

query <- testing[,1:16]
queryLabels <- testing$Type
queryLabels <- as.character(queryLabels)

res <- FNN::knn(data, query, dataLabels, k=5,algorithm=c("kd_tree"),prob=TRUE)

enseRes <- data.frame(KNNProb=attr(res,"prob"))
enseRes$RealLabel <- queryLabels
enseRes$KNNpred <- as.character(res)

clusters <- unique(enseRes$RealLabel)
performance <- data.frame()
pf1 <- list()
myC <- 1
for(idx in clusters){
  
  a <- which(enseRes$KNNpred==idx)
  tmp <- enseRes[a,]
  Nj <- nrow(tmp)
  predicted_res <- table(tmp$RealLabel)
  plot_data <- data.frame(ClusterSize=nrow(tmp),t(predicted_res))
  colnames(plot_data)[2] <- idx
  colnames(plot_data)[3] <- 'PredictedCellType'
  plot_data$Freq <- 100*plot_data$Freq/nrow(tmp)
  str <- paste('Gated cell type: ',idx,sep='')
  pf1[[myC]] <- ggplot(plot_data, aes(x=PredictedCellType,y=Freq,size=Freq,color=PredictedCellType))+
    geom_point(alpha=0.5)+
    geom_linerange(aes(x=PredictedCellType, ymin=0, ymax=Freq),size=1)+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle=45, hjust=1,size = 10),
          axis.text.y = element_text(hjust=1,size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=10),aspect.ratio = 0.5)+ ylim(0,100)+
    ggtitle(str)+scale_colour_brewer(palette = "Set2")
  myName <- which.max(table(tmp$RealLabel))
  s1 <- which(names(table(tmp$RealLabel))==idx)
  Nij <- table(tmp$RealLabel)[s1]
  predictedName <- idx
  a <- which(enseRes$RealLabel==predictedName)
  tmp <- enseRes[a,]
  Ni <- nrow(tmp)
  df <- data.frame(CellType=paste(idx,sep=''),PredictedCellType=predictedName,
                   Recall=Nij/Ni,
                   Precision=Nij/Nj,
                   F1=2*(Nij/Ni)*(Nij/Nj)/((Nij/Ni)+(Nij/Nj)))
  performance <- rbind(performance,df)
  myC <- myC + 1
}

#FOLD 2: train using the previous testing data
#projection based on k-nn
data <- testing[,1:16]
dataLabels <- testing$Type 
dataLabels <- as.character(dataLabels)

query <-  training[,1:16]
queryLabels <- training$Type
queryLabels <- as.character(queryLabels)

res <- FNN::knn(data, query, dataLabels, k=5,algorithm=c("kd_tree"),prob=TRUE)

enseRes <- data.frame(KNNProb=attr(res,"prob"))
enseRes$RealLabel <- queryLabels
enseRes$KNNpred <- as.character(res)

clusters <- unique(enseRes$RealLabel)

pf2 <- list()
myC <- 1
for(idx in clusters){
  
  a <- which(enseRes$KNNpred==idx)
  tmp <- enseRes[a,]
  Nj <- nrow(tmp)
  predicted_res <- table(tmp$RealLabel)
  plot_data <- data.frame(ClusterSize=nrow(tmp),t(predicted_res))
  colnames(plot_data)[2] <- idx
  colnames(plot_data)[3] <- 'PredictedCellType'
  plot_data$Freq <- 100*plot_data$Freq/nrow(tmp)
  str <- paste('Gated cell type: ',idx,sep='')
  pf2[[myC]] <- ggplot(plot_data, aes(x=PredictedCellType,y=Freq,size=Freq,color=PredictedCellType))+
    geom_point(alpha=0.5)+
    geom_linerange(aes(x=PredictedCellType, ymin=0, ymax=Freq),size=1)+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle=45, hjust=1,size = 10),
          axis.text.y = element_text(hjust=1,size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=10),aspect.ratio = 0.5)+ ylim(0,100)+
    ggtitle(str)+scale_colour_brewer(palette = "Set2")
  myName <- which.max(table(tmp$RealLabel))
  s1 <- which(names(table(tmp$RealLabel))==idx)
  Nij <- table(tmp$RealLabel)[s1]
  predictedName <- idx
  a <- which(enseRes$RealLabel==predictedName)
  tmp <- enseRes[a,]
  Ni <- nrow(tmp)
  df <- data.frame(CellType=paste(idx,sep=''),PredictedCellType=predictedName,
                   Recall=Nij/Ni,
                   Precision=Nij/Nj,
                   F1=2*(Nij/Ni)*(Nij/Nj)/((Nij/Ni)+(Nij/Nj)))
  performance <- rbind(performance,df)
  myC <- myC + 1
}


#save cell misclassification plots for fold 1 and fold 2 
#cell type classification performance can be found in performance_fold1 and performance_fold2 data frames
combined_plot_f1 <- ggarrange(pf1[[1]],pf1[[2]],pf1[[3]],pf1[[4]],pf1[[5]],pf1[[6]],pf1[[7]],pf1[[8]],
                              nrow=3,ncol=3)
myfile<-paste('plot15_f3.pdf')
pdf(myfile,onefile = TRUE)
print(combined_plot_f1)
dev.off()

combined_plot_f2 <- ggarrange(pf2[[1]],pf2[[2]],pf2[[3]],pf2[[4]],pf2[[5]],pf2[[6]],pf2[[7]],pf2[[8]],
                              nrow=3,ncol=3)
myfile<-paste('plot15_f4.pdf')
pdf(myfile,onefile = TRUE)
print(combined_plot_f2)
dev.off()

#load "toy" data from leukemia patients with ~1M cells from 6 patients
load('leukemiaData.Rdata')

#retrieve patient names
patientsList <- unique(leukemiaData$ID)
PCALeukemiaData <- data.frame()
#from each patient retrieve ~5K cells and project them into the healthy PCA space
#use the following sets of markers for annotation - notice only difference in HLA-DR versus HLA.DR transformed from R sometimes....
#TODO here, the problem with HLA-DR needs some fixing
selectedMarkers <- c('CD45','CD33','CD38','CD34','CD123','CD20','CD3','CD8a','CD4','CD7','CD56','CD66b','CD64','HLA.DR','CD14','CD16')
for(eachPatient in patientsList){
  
  a <- which(leukemiaData$ID==eachPatient)
  patientData <- leukemiaData[a,]
  #select random cells
  a <- sample(nrow(patientData))
  a_select <- a[1:5000]
  subPatientData <- patientData[a_select,]
  tmp <- subPatientData[,selectedMarkers]
  #produce PCA coordinates
  s.sc <- scale(tmp, center= pca$center)
  s.pred <- s.sc %*% pca$rotation
  PCAtransformed <- data.frame(PC1=s.pred[,1],PC2=s.pred[,2],tmp,ID=eachPatient,Survival=unique(subPatientData$Survival))
  PCALeukemiaData <- rbind(PCALeukemiaData,PCAtransformed)
  
}

#visualise cells from leukemia patients and color them based on CD45 marker expression
#for illustration purposes only we scale the values to 0-1,for real analysis better to z-score the values globally
PCALeukemiaData$CD45 <- (PCALeukemiaData$CD45-min(PCALeukemiaData$CD45))/(max(PCALeukemiaData$CD45)-min(PCALeukemiaData$CD45))
PCALeukemiaData$CD34 <- (PCALeukemiaData$CD34-min(PCALeukemiaData$CD34))/(max(PCALeukemiaData$CD34)-min(PCALeukemiaData$CD34))

leukemiaProjection1 <- ggplot()+
  geom_point(data=PCALeukemiaData,  aes(x = PC1, y = PC2,colour=CD45),size = 0.2,alpha=0.6)+
  scale_color_gradientn("CD45",colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50),
                        limits=c(0,1),breaks = c(0,0.5,1),values=c('0','0.5','1'))+
  geom_point(data=centroids,  aes(x = PC1, y = PC2, shape = Label,group=1),size = 2.8)+
  geom_point(data=centroids,  aes(x=PC1,y=PC2,size=Radius),alpha = 0.2)+
  guides(shape=FALSE)+
  guides(size=FALSE)+
  facet_wrap(~ID)+
  theme_minimal()+
  xlab('PC1')+
  ylab('PC2')+
  scale_shape_manual(values = rep(1:8, len = 8))+
  theme(axis.text.y = element_text( size = 10 ),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 10),
        legend.position = "bottom",
        legend.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),aspect.ratio = 1)


leukemiaProjection2 <- ggplot()+
  geom_point(data=PCALeukemiaData,  aes(x = PC1, y = PC2,colour=CD34),size = 0.2,alpha=0.6)+
  scale_color_gradientn("CD34",colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50),
                        limits=c(0,1),breaks = c(0,0.5,1),values=c('0','0.5','1'))+
  geom_point(data=centroids,  aes(x = PC1, y = PC2, shape = Label,group=1),size = 2.8)+
  geom_point(data=centroids,  aes(x=PC1,y=PC2,size=Radius),alpha = 0.2)+
  guides(shape=FALSE)+
  guides(size=FALSE)+
  facet_wrap(~ID)+
  theme_minimal()+
  xlab('PC1')+
  ylab('PC2')+
  scale_shape_manual(values = rep(1:8, len = 8))+
  theme(axis.text.y = element_text( size = 10 ),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 10),
        legend.position = "bottom",
        legend.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),aspect.ratio = 1)

#save the plots - in the same way all markers can be visualised
myfile<-paste('plot16.pdf')
pdf(myfile,onefile = TRUE)
print(leukemiaProjection1)
dev.off()

myfile<-paste('plot17.pdf')
pdf(myfile,onefile = TRUE)
print(leukemiaProjection2)
dev.off()

#use the healthy hematopoietic classifierb (lda model variable) to annotate cells from leukemia patients
#due to R changing of column names, marker HLA.DR must be renamed to HLA-DR, trained model uses name 'HLA-DR'
#this is done manually, but this part has to be improved
colnames(leukemiaData)[16] <- 'HLA-DR'
selectedMarkers <- c('CD45','CD33','CD38','CD34','CD123','CD20','CD3','CD8a','CD4','CD7','CD56','CD66b','CD64','HLA-DR','CD14','CD16')
#for the annotation part we use subset of the original data frame
surfaceMarkerExpr <-leukemiaData[,selectedMarkers]

#predict the class, here we show an example using the developed LDA model because it is fast
#to use KNN model un-comment the following chunk of code, takes ~40 min in a commodity computer
#data <- HD_annotated[,1:16]
#dataLabels <- HD_annotated$Type
#dataLabels <- as.character(dataLabels)
#query <- surfaceMarkerExpr
#res <- FNN::knn(data, query, dataLabels, k=5,algorithm=c("kd_tree"),prob=TRUE)
#enseRes <- data.frame(KNNProb=attr(res,"prob"))
#enseRes$KNNpred <- as.character(res)

predicted_class <- predict(model,surfaceMarkerExpr)
#estimate the posterior probabilities
prob_matrix<-matrix(0L, nrow = nrow(tumor_expr), ncol =1)
prob_matrix<-as.matrix(apply(predicted_class$posterior,1,max))
#add new columns in the data frame about the predicted class and the probabilities
df <- leukemiaData
df$PredictedClass <- predicted_class$class
df$Prob <- prob_matrix
#filter the predictions based on the probability threshold, thats after tuning - must be filled carefully
uncertainPredictions <- df %>% filter(df$Prob<0.85)
goodPredictions <- df %>% filter(df$Prob>=0.85)
#LowConfidence predictions are the uncertain predictions that we do not consider for further analysis
#this is TODO work, analysing the uncertain predictions is also very important to study tumour heterogeneity
uncertainPredictions$PredictedClass <- 'LowConfidence'
#summary statistics of annotated cells with good quality
leukemiaData_annotated <- goodPredictions
summaryStats <- as.data.frame(100*table(leukemiaData_annotated$PredictedClass)/nrow(leukemiaData_annotated))
colnames(summaryStats)[1] <- 'CellType' 
#discretize the survival data
leukemiaData_annotated$Survival <- ifelse(leukemiaData_annotated$Survival==60,'Survivors','Deceased')
#and make some visualisations of the summary, similar to Figure 2 of the manuscript but with stacked bar plot instead of boxplot
cell_count_class <- leukemiaData_annotated %>% 
  group_by(Survival,PredictedClass) %>% 
  summarise(cells = n())

#uncomment to remove scientific notation in r
#options(scipen=999)

abundancePlot <- ggplot(cell_count_class, aes(x=PredictedClass,y=cells,fill=Survival)) + 
  geom_bar(stat="identity") +
  theme_bw()+
  theme(axis.text.y = element_text( size = 12 ),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text( size = 12 ),
        strip.text = element_text(size = 14,face='bold',lineheight=1),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        strip.background = element_rect(colour = "black", fill = "white"),
        plot.title = element_text(size = 14, face = "bold",hjust = 0.5),aspect.ratio = 0.4)+
  ylab('Cells in population')+
  scale_fill_manual(values=c("#D73027","#74ADD1","#E6AB02"),name='')+
  theme(panel.spacing = unit(1.8, "lines"))+coord_flip()
  

myfile<-paste('plot18.pdf')
pdf(myfile,onefile = TRUE)
print(abundancePlot)
dev.off()
  
#remove the posterior probability column, it is not needed for further analyses
leukemiaData_annotated <- leukemiaData_annotated[ , !(names(leukemiaData_annotated) %in% c('Prob'))]
#also remove marker p4E_BP1 due to unspecific staining for some cases in the study
leukemiaData_annotated <- leukemiaData_annotated[ , !(names(leukemiaData_annotated) %in% c('p4E_BP1'))]
#rename column PredictedClass to CellType and save the data for further use
colnames(leukemiaData_annotated)[ncol(leukemiaData_annotated)] <- 'CellType'

save(leukemiaData_annotated,file='leukemiaDataAnnot.Rdata')
