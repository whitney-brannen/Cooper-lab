# -----------------------------------------
#
# This script wil perform admixture analysis of a vcf file and 
# display results in a floating pie charts by gps coordinates on a map
#
# start with a new R project and only this script and input files in the folder
# or fatal errors will throw
#
#-----------------------------------------

library(adegenet)
library(LEA)
library(dplyr)
library(purrr)
library(reshape2)
library(ggplot2)
library(rnaturalearth)
library(ggmap)

# inputs
input_file_handle = "WGS_subset"
my.colors = c("tomato", "lightblue", "olivedrab", "gold")
gps = read.csv("gps.csv")
datafilt = read.csv("pcs_withregions.csv")

# get list of all IDs
subsites = gps$ID

# do PCA
pc = pca(paste(input_file_handle,".vcf",sep = ""))
tw = tracy.widom(pc)

# scree plot
plot(tw$percentage, pch = 19, col = "blue", cex = .8)

# calculate cross entropy (this will take a while)
project = NULL
project = snmf(paste(input_file_handle,".geno",sep = ""),
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new")

plot(project, col = "blue", pch = 19, cex = 1.2)

# chose K value
my.K = 3

lowest.ce = which.min(cross.entropy(project, K = my.K))

# plot admixture bar graph
barchart(project, K = my.K, run = lowest.ce,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix K=3") -> bp

axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)


# q matrix and header names
qmatrix = as.data.frame(Q(project, K = my.K, run = lowest.ce))

cluster_names = c()
for (i in 1:ncol(qmatrix)) {
  cluster_names[i] = paste("Cluster", i)
}
colnames(qmatrix) = cluster_names

# average admix proportions
qmatrix$ID = datafilt$Region

clusters = grep("Cluster", names(qmatrix))
avg_admix = aggregate(qmatrix[, clusters], list(qmatrix$ID), mean)

avg_admix = melt(avg_admix, id.vars = "Group.1")

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, id, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, Group.1 == id),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = cols)+
    theme_void()
}

# perform pie chart function on all data 
pies = list()
for (i in subsites) {
  pies[[i]] = pie_charts(admix_df = avg_admix, id = i, cols = my.colors) 
}

# get coordinates for mapping
coord.list = list()
for (i in subsites){
  coord.list[[i]] = c(subset(gps, ID == i)$Lon, subset(gps, ID == i)$Lat)
}

# set radius of pie chart
radius = 1

# add to plot
pies.ac = list()

for (i in gps$ID){
    pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

# create basemap of north america
map <- ne_countries(scale = "medium", returnclass = "sf", country = "United States of America")

basemap <- ggplot() + geom_sf(data = map, fill="light grey", color="black") + xlim(-125,-70) + ylim(25,50) + theme_bw()

# display map and add titles
my.title = "ADMIXTURE Proportions of Culex tarsalis for K=3"

pie.map = basemap + pies.ac
pie.map = pie.map + ggtitle(my.title) + xlab("Latitude") + ylab("Longitude")
pie.map


