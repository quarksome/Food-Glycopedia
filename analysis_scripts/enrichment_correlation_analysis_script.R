#import modules
library(RColorBrewer)
library(ComplexHeatmap)
library(gridBase)
library(dplyr)
library(tidyr)
library(corrplot)
library(gplots)

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

#log10 transform mono values for heatmap
monoMat2 <- log10(monoMat[,1:ncol(monoMat)])

#change sample_ID to numerical string '0001'
sampleName <- as.numeric(unlist(NIH_all_mono[,2]))
sampleName <- formatC(sampleName, width = 4, format = "d", flag = "0")
NIH_all_mono <- NIH_all_mono %>%
  mutate(sample_ID=sampleName)

#get cluster membership
hc <- hclust(dist(monoMat), method = 'complete')
k <-  5 #set number of clusters: use determine_optimum_cluster.R script 

clu_df <- data.frame(Cluster = paste('cluster', as.character(cutree(hc, k=k)), sep=""))
clu_df <- cbind(clu_df, sample_ID = NIH_all_mono$sample_ID)
NIH_all_mono2 <- clu_df %>% inner_join(NIH_all_mono)

#get cluster vs food_group frequency matrix for Enrichment Analysis
datalist = list() 
for (cluster in unique(NIH_all_mono2$Cluster)) {
  c_vec <- c()
  for (food_grp in unique(NIH_all_mono2$food_group)) {
    count <- sum(NIH_all_mono2$Cluster == cluster & NIH_all_mono2$food_group == food_grp)
    c_vec <- append(c_vec, count)
  }
  names(c_vec) <- unique(NIH_all_mono2$food_group)
  c_vec <- data.frame(c_vec,
                      row.names = unique(NIH_all_mono2$food_group)
  )
  colnames(c_vec) <- c(cluster)
  datalist <- append(datalist, c_vec)
}
freq_table = do.call(cbind, datalist)
rownames(freq_table) <- unique(NIH_all_mono2$food_group)

freq_table <- rbind(freq_table, colSums(freq_table))
rownames(freq_table)[10] <- "sum"
freq_table <- cbind(freq_table, rowSums(freq_table))
colnames(freq_table)[k+1] <- "sum"

#generating p-vals and Enrichment Factors from hypergeometric test 
A <- function(x,y,z,N) {phyper(x-1,y,N-y,z, lower.tail = FALSE)}

pVal_list = list()
EF_list = list()

for (i in c(1:length(unique(clu_df[["Cluster"]])))) {
  cl <- paste('cluster',as.character(i), sep='')
  # print(cl)
  x <- freq_table[,cl]
  y <- tail(freq_table[,cl], n=1)
  N <- freq_table['sum','sum']
  z <- freq_table[,'sum']
  
  pval <- A(x,y,z,N)
  pval <- data.frame(pval)
  colnames(pval) <- cl
  pVal_list <- append(pVal_list, pval)
  
  EF <- (x/z)/(y/N)
  EF <- data.frame(EF)
  colnames(EF) <- cl
  EF_list <- append(EF_list, EF)
}
pval_df <- do.call(cbind, pVal_list)
rownames(pval_df) <- rownames(freq_table)
pval_df <- pval_df[!rownames(pval_df) %in% c('sum'), ]    #delete row sum

EF_df <- do.call(cbind, EF_list)
rownames(EF_df) <- rownames(freq_table)
EF_df <- EF_df[!rownames(EF_df) %in% c('sum'), ]    #delete row sum

#get EF values where EF>1
EF_label <- ifelse(EF_df > 1, as.character("*"), "")

#adjusting p-value for multiple comparisons
pval_adj <- pval_df %>% 
  as.matrix %>% 
  as.vector %>% 
  p.adjust(method='fdr') %>% 
  matrix(ncol=length(unique(clu_df[["Cluster"]])))
colnames(pval_adj) <- colnames(pval_df)
rownames(pval_adj) <- rownames(pval_df)

#get -log10(p_val_adjusted)
pval_adj2 <- -log10(pval_adj)
pval_adj2[pval_adj2 == "Inf"] <- NA

#get adj p-val < 0.05 
pval_label <- ifelse(pval_adj2 > -log10(0.05), as.character("*"), "")

#color palette for EnrichmentFactor and pval heatmaps
pal <- brewer.pal(8, "Blues")

#generate heatmap object for adjusted p-values
pval_hm <- heatmap.2(pval_adj2, Rowv = FALSE, Colv = FALSE, trace = 'none', density.info = 'none', 
                     dendrogram = 'none', col=pal, colsep=1:nrow(pval_adj2), rowsep=1:nrow(pval_adj2), 
                     sepcolor = 'black', cellnote = pval_label, notecol='red', notecex = 2.5,
                     margins = c(6,22),
                     key.xlab = "-log10(adj_p_val)")

#generate heatmap object for Enrichment Factors
#add significance marker at pval < 0.05
EF_hm <- heatmap.2(EF_df, Rowv = FALSE, Colv = FALSE, trace = 'none', density.info = 'none',
                   dendrogram = 'none', col=pal, colsep=1:nrow(EF_df), rowsep=1:nrow(EF_df),
                   sepcolor = 'black', cellnote = pval_label, notecol='red', notecex = 2.5,
                   cexCol = 1.5, margins = c(6,22),
                   key.xlab = "Enrichment Factor")

#p_value function for CORRELATION MATRIX (WITH FDR ADJUSTMENT)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  #adjusting p-value
  p.mat <- p.mat %>% 
    as.matrix %>% 
    as.vector %>% 
    p.adjust(method='fdr') %>% 
    matrix(ncol=n)
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  
  p.mat
}

#plot correlation matrix heatmap with significance
p.mat <- cor.mtest(monoMat)
cor_mono <- cor(monoMat)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corr_hm <- corrplot(cor_mono, method="color", col=col(200),  
                    type="upper", order="original", 
                    addCoef.col = "black", number.cex = 0.6, # Add coefficient of correlation
                    tl.col="black", tl.srt=45, #Text label color and rotation
                    # Combine with significance
                    p.mat = p.mat, sig.level = 0.05, insig = "blank", 
                    # hide correlation coefficient on the principal diagonal
                    diag=FALSE,
                    addgrid.col = 'black'
)