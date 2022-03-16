#BEGIN COPYRIGHT NOTICE

#generateFeatures -- (c) 2022 Dimitrios Kleftogiannis -- UiB and CCBIO
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
#This code take as input annotated cells from leukemia patients and generates DREMI features.
#It can be used to make machine learning modelling and to develop predicted models of survival.
#For the latter tasks, different machine learning models can be used, and there is no clear answer what is best.
#In the paper we used glmnet and mRMR feature selection implemented by packages glmnet and mRMRe.
#But i would be happy to see improvements by more advanced approaches.
#The data used in this tutorial are only for illustration purposes.

#RUNNING
#We provide a step-by-step execution example, this is useful for illustration and to reproduce parts of the analysis. 
#Remember that this is a minimal example. To obtain the input dataplease contact Dimitrios Kleftogiannis.

#DEPENDENCIES
#To run the tutorial without problem, is required to install and load all packages listed in implementGatingStrategy example
#additional functions and packages are required for DREMI estimations
source('utilityFuncs.R')
       
library(pracma)
library(pmpp)
library(matlab)
library(ggExtra)
library(minpack.lm)

#load a "toy" dataset
load('leukemiaDataAnnot.Rdata')

#define the set of signaling markers
mySignalingSet <- c("CyclinB1","pRB","pNFkB","pErk","pSTAT1",
           "pSTAT3","pCREB", "pHistone3","CASP3","pS6",
           "pSTAT5","pAkt",'pAxl','pP38')

#generate all possible combinations of the previous markers
#remember here that there is no need to consider the reverse relationships, A->B is the same as B->A
allRelation <- combn(mySignalingSet, 2)
#concatenate
b <- paste(allRelation[1,],'-',allRelation[2,],sep='')
#and make a list for every cell type, we have 8 cell type
#remember here that there is no need to consider the reverse relationships, A->B is the same as B->A
selectFeatList <- list(b,b,b,b,b,b,b,b)
#dewfine the cell types in the study, we will process each of them separately
CellTypes <- unique(leukemiaData_annotated$CellType)
#initiate counter for cell types
myC <- 1
#set an upper limit on the number of cells to use. Just to speedup processing 
Ncells <- 20000
#initiate list to store cell type-specific feature matrices. We have 8 cell types and for each we will compute 91 features with DREMI
featureMatrices <- list()
#iterate over the cell types
for(myClass in CellTypes){

  summaryData <- data.frame()
  myFeat <- selectFeatList[[myC]]
  msg<-paste('Working with cell-type: ',myClass, '\n',sep='')
  cat(msg)
  #iterate over the patient data
  for(myID in unique(leukemiaData$ID)){
    
  
    #for(myFile in c(1,6)){
    msg<-paste('Processing features for patient: ',myID,'\n',sep='')
    cat(msg)
    
    #select the data of each patient for the particular cell type we are studying
    mySel <- which(leukemiaData_annotated$ID==myID & leukemiaData_annotated$CellType==myClass)
    patientData <- leukemiaData_annotated[mySel,]
    #for cases where we have less than 500 cells per cell type, do not provide DREMI estimation and flag the features 
    if(nrow(patientData)<500){
      str <- paste('       Dangerous estimation for:', myID,' cell-type: ',myClass,sep='')
      cat(str)
      #use flag -1 everywhere
      featVector  <- matrix(-1,nrow = 1,ncol=length(myFeat))
      
    }else{
      #select a smaller random subset of cells to speedup, if we have more than 50000 cells per cell type.
      #this can be tuned by user
      if(nrow(patientData)>50000){
        a <- sample(nrow(patientData))
        asel <- a[1:Ncells]
        patientData <- patientData[asel,]
      }
      #store the DREMI scores
      featVector  <- matrix(-1,nrow = 1,ncol=length(myFeat))
      #another internal counter
      myK <- 1
      #parse one feature at a time, considering combinations of 2 markers that are already generated
      for(myRel in myFeat){
        #split the pair into the markers
        a <- strsplit(myRel, "-")
        a <- a[[1]]
        #two markers
        markerX <- a[1]
        markerY <- a[2]
        #generate two vectors of data
        Y <- patientData[,markerY]
        X <- patientData[,markerX]
        #remove outliers
        X <- remove_outliers(X)
        Y <- remove_outliers(Y)
        #combine
        mB<- as.matrix(cbind(X,Y))
        #remove if one of the two values is outlier with NA
        mB <- na.omit(mB)
        #sometimes the staining is not good and the markers have always zero values
        #we can to avoid this case, and we filter as follows:
        test1 <- max(mB[,1], na.rm = TRUE)
        test2 <- max(mB[,2],na.rm=TRUE)
        if(round(test1)==0 | round(test2)==0){
          dremi <- -1 # or impute the value using the distribution of other features of the survival group
        }else{
          myLabel <- myID
          maxy <- 0;
          myCutoffs <- find_data_cutoffs(mB, 50, 255)
          minx1 <- myCutoffs[[1]]
          miny1 <- myCutoffs[[2]]
          maxx1 <- myCutoffs[[3]]
          maxy1 <- myCutoffs[[4]]
          if(maxy1>maxy){
            maxy <- maxy1;
          }
          #adjust the noise threshold accordinly - default is 0.8
          noise_threshold <- 0.8;
          set_maxy <- 0
          #set to 1 to generate DREVI plot - this part is not used in the current analysis
          makePlot <- 0
          #function that returns the DREMI score
          dremi <- compute_dremi(mB,markerX,markerY,noise_threshold,set_maxy,maxy,makePlot,myLabel)
          #the procedure returns a list of output arguments, but for this example we only need the DREMI score which is the first
          dremi <- dremi[[1]]
        }
        #store it in the feature vector
        featVector[1,myK]<- dremi
        myK <- myK + 1
      }
    }#move to the next patient
    featVector <- as.data.frame(featVector)
    colnames(featVector) <- myFeat
    featVector$ID <- myID
    featVector$Class <- myClass
    featVector$Survival <- unique(patientData$Survival)
    summaryData <- rbind(summaryData,featVector)
  }#move to next cell type
  
  featureMatrices[[myC]] <- summaryData
  myC <- myC +1
  rm(summaryData)
}

#make a simple mds 2D plot for visualisation of the feature vector
#visualise the first cell type, but equally works for all other
str <- paste('MDS plot for cell type: ',  CellTypes[1],sep='')
summaryData <- featureMatrices[[1]]
mds <- plotMDS(t(summaryData[,1:91]), plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, Group=summaryData$Survival,ID=summaryData$ID)
out1 <- ggplot(ggdf, aes(x = MDS1, y = MDS2,color=Group,label=ID)) +
  geom_point(size = 2.4, alpha = 1.8)+
  theme_classic()+
  ylab('MDS2')+
  xlab('MDS1')+
  ggtitle(str)+
  geom_label_repel(aes(label = ID ),size = 1.8, max.overlaps = 22,)+
  scale_color_manual(values=c("#D73027","#74ADD1","#E6AB02"),name='')+
  theme(axis.text.y = element_text( size = 12 ),
        axis.text.x = element_text( size = 12 ),
        axis.title.y = element_text( size = 12 ),
        axis.title.x = element_text( size = 12 ),
        legend.position = "bottom",legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size=5),nrow=1))+
  guides(fill = guide_legend(override.aes = list(size=5),ncol=1))

myfile<-paste('plot19.pdf')
pdf(myfile,onefile = TRUE)
print(out1)
dev.off()
