#Script to determine optimum number of clusters using NbClust package

#import modules
library(dplyr)
library(tidyr)
library(NbClust)
library(factoextra)

set.seed(123)

setwd("NIH_Glycopedia")

#read excel file
NIH_all_mono <- readxl::read_excel("NIH_all_mono_newFGv3.xlsx", 
                                   col_types = c("text", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric"))

#setup matrices/DF
##extract list of food_groups
foodGrpList <- NIH_all_mono[,1]
##extract mono values only
monoMat <- NIH_all_mono[, -(1:2)]

#replace missing values with 0
monoMat[is.na(monoMat)] <- 0

#replace 0 values with min.value / 5 (for each column)
monoMat[1:ncol(monoMat)] <- lapply(monoMat[1:ncol(monoMat)],
                                   function(x) replace(x, (x==0)|(x<0), min(x[x>0], na.rm = TRUE)/5))

#determine # of clusters
res.nbclust <- NbClust(monoMat, distance = "euclidean",
                       min.nc = 3, max.nc = 12, 
                       method = "complete", index ="all")

#plot NbClust results
fviz_nbclust(res.nbclust) + ggtitle("NbClust's optimal number of clusters")

#plot WSS plot
fviz_nbclust(monoMat, FUN=hcut, method='wss')

