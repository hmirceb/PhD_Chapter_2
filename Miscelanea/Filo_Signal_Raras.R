library(dplyr)
library(emmeans)
library(visreg)
library(DHARMa)
library(lme4)
library(reshape2)
library(readxl)
library(ggplot2)
library(PhyloMeasures)
library(popbio)
library(scales)
library(ggtree)
library(ggeffects)
library(ggpubr)
library(directlabels)
library(gridExtra)

setwd('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 1 - Diversidad y vulnerabilidad de plantas del Pirineo/Figuras y tablas')
load('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 1 - Diversidad y vulnerabilidad de plantas del Pirineo/Analisis/Resultados/PD_Pirineo_y_rareza.RData')
load('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 1 - Diversidad y vulnerabilidad de plantas del Pirineo/Analisis/Data/Datos_Div_Vul.RData')
load("/home/hector/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")

library("ape")
source("code.R")

#SOME PARAMETERS... 
lambda0 <- 0.1   #rate parameter of the proposal 
se      <- 0.5   #standard deviation of the proposal
sim     <- 10000 #number of iterations
thin    <- 10    #we kept only each 10th iterate 
burn    <- 100   #100 iterates are burned-in



#CALCULATE DELTA A
deltaA <- delta(trait,tree,lambda0,se,sim,thin,burn)
print(deltaA)

#DRAW THE TREE...
par(mfrow=c(1,2))
tree$tip.label <- rep("",ns)
plot(tree,main="SCENARIO A")
ar <- ace(trait,tree,type="discret",method="ML",model="ARD")$lik.anc
nodelabels(pie = ar, cex = 1,frame="n") 
mtrait <- matrix(0,ncol=3,nrow=ns)
for ( i in 1:ns) {
  mtrait[i,trait[i]] <- 1
}
tiplabels(pie=mtrait,cex=0.5)

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait,tree,lambda0,se,sim,thin,burn)
}
p_value <- sum(random_delta>deltaA)/length(random_delta)
boxplot(random_delta)
abline(h=deltaA,col="red")
p_value