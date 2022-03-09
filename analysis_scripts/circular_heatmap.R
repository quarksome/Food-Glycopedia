#import modules
library(circlize)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)
library(gridBase)
library(dplyr)
library(tidyr)
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

#set color pallette for circular heatmaps
##for mono values
col_fun1 = colorRamp2(c(min(monoMat2), max(monoMat2)), c("white", "blue"))

##for food groups
col_foodGrp = c("Fruits" = "#386cb0", "Vegetables" = "#7fc97f", "Grain Products" = "#bf5b17",
                "Sugars, Sweets, and Beverages" = "#ffff99", 
                "Beans, Peas, Other Legumes, Nuts, Seeds" = "#F15CF8",
                "Milk and Milk Products" = "#fdc086", "Meat, Poultry, Fish, and Mixtures" = "#beaed4",
                "Fats, Oils, and Salad Dressings" = '#463B38', 'Eggs' = '#FA030E')

#setup legend for circular heatmap
lgd_mono = Legend(title = "log(monosaccharide)", col_fun = col_fun1)
lgd_foodGrp = Legend(title = "Food group", at = names(col_foodGrp),
                     legend_gp = gpar(fill = col_foodGrp))

#get cluster membership
hc <- hclust(dist(monoMat), method = 'complete')
k <-  5 #set number of clusters: use determine_optimum_cluster.R script 

#function for plotting circular heatmap (wrapper function)
circlize_plot = function() {
  
  circos.clear()
  circos.par(gap.after = c(15))
  circos.heatmap.initialize(monoMat, cluster = hc,
                            dend.callback = function(dend, m, si) {
                              color_branches(dend, k = k, col = 1:k)})
  circos.heatmap(foodGrpList, col = col_foodGrp, track.height = 0.03) 
  circos.heatmap(monoMat2, col = col_fun1, dend.side = "inside",
                 dend.track.height = 0.5,
                 bg.border = "black",
                 track.height = 0.4
  )
  
  circos.track(track.index = get.current.track.index()-1, panel.fun = function(x, y) {
    cn = colnames(monoMat)
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                1:n - 0.3, rev(cn),
                cex = 0.5, adj = c(0, 0.3), facing = "inside")
  }, bg.border = NA)
  
}

#plotting
dev.off()
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region

pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circlize_plot()
upViewport()

h = dev.size()[2]
lgd_list = packLegend(lgd_mono, lgd_foodGrp, max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")