#BEGIN COPYRIGHT NOTICE

#implementGatingStrategy -- (c) 2022 Dimitrios Kleftogiannis -- UiB and CCBIO
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
#This code implements a gating strategy similar to the one presented in Supplementary Figure 1 of the paper.
#It can be used to study and visualise healthy hematopoietic cell types. 
#The data used in this tutorial are only for illustration purposes.

#RUNNING
#We provide a step-by-step execution example, this is useful for illustration and to reproduce parts of the analysis. 
#Remember that this is a minimal example. To obtain the input dataplease contact Dimitrios Kleftogiannis.

#DEPENDENCIES
#To run the tutorial without problem, is required to install the following packages. 

#load packages after installation
pkgs2load <- c('readxl','plyr', 'dplyr', 'magrittr', 'tidyr', 'flowCore', 'FlowSOM', 
               'data.table', 'Rtsne', 'ggplot2', 'gridExtra', 'survival', 'glmnet', 
               'ggpubr', 'survminer', 'readr', 'mixOmics', 'RColorBrewer', 
               'corrr', 'tidygraph', 'ggraph', 'sva','matrixStats','reshape2','limma','ggrepel',
               'ConsensusClusterPlus','Rphenograph','pheatmap','dbscan',
               'flowDensity','ggridges','caret','corrplot','dendextend')
sapply(pkgs2load, require, character.only = TRUE)

#load a "toy" healthy dataset named dataValues
load('HD_input.Rdata')
#generate the expression matrix of the cell dataset
b <- dataValues@exprs

#use density estimation of CD45 to distinguish lymphoid from myeloid cells 
CD45_cut <- deGate(obj = dataValues,channel = c('CD45'),all.cuts=T)

#generate the first plot of the gating strategy
myfile<-paste('plot1.pdf')
pdf(myfile,onefile = TRUE)
plot(density(b[,'CD45']),main='',ylab="Density (healthy cells)",xlab='CD45 expression',
     xlim=c(0,6),lwd=3,cex.axis=1.48,cex.lab = 1.48)
abline(v=c(CD45_cut[2]),col=c(2),lwd=2,lty=3)
dev.off()

#identify lymphoid vs. myeloid compartments given the expression of CD45
part1 <- which(b[,c('CD45')]>=CD45_cut[2])
Lymphoid <- b[part1,]
Myeloid <- b[-part1,]

#take the Lymphoid compartment and gate to identify different immune cell types
Lymphoid <- flowFrame(Lymphoid)
#modify appropriately the percentiles
CD3pCD20p <- flowDensity(Lymphoid,channels = c('CD20','CD3'),position = c(T,T),
                         percentile=c(0.80,0.20),use.percentile=c(T,T))
CD3pCD20n <- flowDensity(Lymphoid,channels = c('CD20','CD3'),position = c(F,T),
                         percentile=c(0.85,0.38),use.percentile=c(T,T))
CD3nCD20p <- flowDensity(Lymphoid,channels = c('CD20','CD3'),position = c(T,F),
                         percentile=c(0.90,0.38),use.percentile=c(T,T))
CD3nCD20n <- flowDensity(Lymphoid,channels = c('CD20','CD3'),position = c(F,F),
                         percentile=c(0.88,0.37),use.percentile=c(T,T))

#generate the gating plot -> 
myfile<-paste('plot2.pdf')
pdf(myfile,onefile = TRUE)
plotDens(Lymphoid,c('CD20','CD3'),ylab='CD3',xlab='CD20',main='',cex.axis=1.48,cex.lab = 1.48)
lines(CD3nCD20n@filter,type="l",lwd=2)
lines(CD3pCD20n@filter,type="l",lwd=2)
lines(CD3nCD20p@filter,type="l",lwd=2)
dev.off()

#gating for NK cells using CD3 negative and CD20 negative cells
CD3nCD20n <- getflowFrame(CD3nCD20n)
NKa2 <- flowDensity(CD3nCD20n,channels = c('CD16','CD56'),position = c(T,T),
                    percentile=c(0.90,0.74),use.percentile=c(T,T)) 
NKa4 <- flowDensity(CD3nCD20n,channels = c('CD16','CD56'),position = c(F,T),
                    percentile=c(0.64,0.74),use.percentile=c(T,T)) 

#produce the third plot
myfile<-paste('plot3.pdf')
pdf(myfile,onefile = TRUE)
plotDens(CD3nCD20n,c('CD16','CD56'),ylab='CD56',xlab='CD16',main='',cex.axis=1.48,cex.lab = 1.48)
lines(NKa2@filter,type="l",lwd=2)
lines(NKa4@filter,type="l",lwd=2)
#plot(CD3nCD20n,NKa3)
dev.off()

#NK cell identification is not perfect and it requires further improvements
#here we just merge 3 different subcategories into one; further work is required
a1 <- getflowFrame(NKa2)
a3 <- getflowFrame(NKa4)
a1 <- a1@exprs
a3 <- a3@exprs
tmp <- rbind(a1,a3)
tmp <- flowFrame(tmp)
NK <- tmp@exprs

#gating for T cells - cell type defining marker is CD3 and absence of CD20
CD3pCD20n <- getflowFrame(CD3pCD20n)

#double positive and double negative T-cells are excluded from the analysis, although can be identified, if so uncomment the following code
#TcellsDP <- flowDensity(CD3pCD20n,channels = c('CD8a','CD4'),position = c(T,T),
#                        percentile=c(0.80,0.50),use.percentile=c(T,T)) 
#TcellsDN <- flowDensity(CD3pCD20n,channels = c('CD8a','CD4'),position = c(F,F),
#                        percentile=c(0.68,0.35),use.percentile=c(T,T)) 

#use CD8 and CD4 for further T-cell gating
CD8posCD4neg <- flowDensity(CD3pCD20n,channels = c('CD8a','CD4'),position = c(T,F),
                            percentile=c(0.74,0.365),use.percentile=c(T,T)) 
CD8negCD4pos <- flowDensity(CD3pCD20n,channels = c('CD8a','CD4'),position = c(F,T),
                            percentile=c(0.70,0.42),use.percentile=c(T,T)) 

#plot 4: about T-cells
myfile<-paste('plot4.pdf')
pdf(myfile,onefile = TRUE)
plotDens(CD3pCD20n,c('CD8a','CD4'),ylab='CD4',xlab='CD8a',main='',cex.axis=1.48,cex.lab = 1.48)
lines(CD8posCD4neg@filter,type="l",lwd=2)
lines(CD8negCD4pos@filter,type="l",lwd=2)
dev.off()

#with the use of CD7 we can have more accurate annotation of T-cytotoxic and T-helper cells
CD8negCD4pos <- getflowFrame(CD8negCD4pos)
CD8posCD4neg <- getflowFrame(CD8posCD4neg)

TcellsCytotoxic <- flowDensity(CD8posCD4neg,channels = c('CD8a','CD7'),position = c(T,T),
                               percentile=c(0.10,0.10),use.percentile=c(T,T)) 

TcellsHelper <- flowDensity(CD8negCD4pos,channels = c('CD4','CD7'),position = c(T,T),
                            percentile=c(0.10,0.10),use.percentile=c(T,T)) 

#generate the correspinding plots
myfile<-paste('plot5.pdf')
pdf(myfile,onefile = TRUE)
#plot(CD8posCD4neg,TcellsCytotoxic)
plotDens(CD8posCD4neg,c('CD8a','CD7'),ylab='CD7',xlab='CD8a',main='',cex.axis=1.48,cex.lab = 1.48)
lines(TcellsCytotoxic@filter,type="l",lwd=2)
dev.off()

myfile<-paste('plot6.pdf')
pdf(myfile,onefile = TRUE)
plotDens(CD8negCD4pos,c('CD4','CD7'),ylab='CD7',xlab='CD4',main='',cex.axis=1.48,cex.lab = 1.48)
lines(TcellsHelper@filter,type="l",lwd=2)
#plot(CD8negCD4pos,TcellsHelper)
dev.off()

#summarize the lympohoid data and make the appropriate conversions for further use from flowFrames
Bcells <- getflowFrame(CD3nCD20p)
TcellsCyto <- getflowFrame(TcellsCytotoxic)
TcellsHelper <- getflowFrame(TcellsHelper)

#we process the myeloid compartment seperately, we use CD34 and CD38 cell type defining markers
#we also note that detecting MPP and CMP cells in healthy donors is quite difficult because they are rare, samples from BM provide better resolution
#please tune the thresholds to collect sufficient number of cells
Myeloid <- flowFrame(Myeloid)
CD34pCD38p <- flowDensity(Myeloid,channels = c('CD34','CD38'),position = c(T,T),
                          percentile=c(0.99865,0.988),use.percentile=c(T,T))
CD34pCD38n <- flowDensity(Myeloid,channels = c('CD34','CD38'),position = c(T,F),
                          percentile=c(0.9975,0.60),use.percentile=c(T,T))
CD34nCD38n <- flowDensity(Myeloid,channels = c('CD34','CD38'),position = c(F,F),
                          percentile=c(0.95,0.50),use.percentile=c(T,T))
CD34nCD38p <- flowDensity(Myeloid,channels = c('CD34','CD38'),position = c(F,T),
                          percentile=c(0.95,0.958),use.percentile=c(T,T))

#generate the last plot
myfile<-paste('plot7.pdf')
pdf(myfile,onefile = TRUE)
plotDens(Myeloid,c('CD34','CD38'),ylab='CD38',xlab='CD34',main='',cex.axis=1.48,cex.lab = 1.48)
lines(CD34pCD38n@filter,type="l",lwd=2)
lines(CD34pCD38p@filter,type="l",lwd=2)
lines(CD34nCD38n@filter,type="l",lwd=2)
lines(CD34nCD38p@filter,type="l",lwd=2)
dev.off()

#summarise the myeloid compartements
CD34pCD38n <- getflowFrame(CD34pCD38n)
MPP <- CD34pCD38n
CD34pCD38p <- getflowFrame(CD34pCD38p)
CMP <- CD34pCD38p
CD34nCD38n <- getflowFrame(CD34nCD38n)
CD34nCD38p <- getflowFrame(CD34nCD38p)

#aggregate all myeloid and lymphoid data,take the expression values and provide indicative names to the data frames
CMP <- CMP@exprs
MPP <- MPP@exprs
Myeloid_CD34NCD38N <- CD34nCD38n@exprs
Myeloid_CD34NCD38P <- CD34nCD38p@exprs
Bcells <- Bcells@exprs
TcellsCyto <- TcellsCyto@exprs
TcellsHelper <- TcellsHelper@exprs
#NK is already processed --> no need to add here

#give the labels to the detected cell types
CMP <- as.data.frame(CMP)
CMP$Type <- 'CMP'
MPP <- as.data.frame(MPP)
MPP$Type <- 'MPP'
#CD34-low/CD38-low
Myeloid_a1 <- as.data.frame(Myeloid_CD34NCD38N)
Myeloid_a1$Type <- 'Granulocytes'
#CD34-low/CD38-high
Myeloid_a2 <- as.data.frame(Myeloid_CD34NCD38P)
Myeloid_a2$Type <- 'Monocytes'
Bcells <- as.data.frame(Bcells)
Bcells$Type <- 'B-cells'
NK <- as.data.frame(NK)
NK$Type <- 'NK'
TcellsCyto <- as.data.frame(TcellsCyto)
TcellsCyto$Type <- 'T-cytotoxic'
TcellsHelper <- as.data.frame(TcellsHelper)
TcellsHelper$Type <- 'T-helper'

save(Bcells,NK,
     TcellsCyto,TcellsHelper,
     MPP,CMP,Myeloid_a1,Myeloid_a2,file='HD_annotated.Rdata')
