#########################################
#MAIN FIGURES
#########################################
library(dartR)
library(adegenet)
library(vcfR)
library(SNPRelate)
library(remotes)
library(gwscaR)
library("readxl")
library(ggplot2)
library(grid)
library(gridExtra)
library("FactoMineR")
library("factoextra")
library(MetBrewer)

#########################################
#FIGURE 1A
#########################################
#SNPRelate
setwd("~/R/barley_collab/parental_selection/")
vcf.fn <- "Supplemental_dataset_1.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "barley.gds", method="copy.num.of.ref")
snpgdsSummary("barley.gds")
genofile <- snpgdsOpen("barley.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) 
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

#For PCA 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#make a data frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
write.csv(tab, "barley_pca.csv")

#Plot Hierarchal Cluster
library("FactoMineR")
library("factoextra")
pca <- read.csv("barley_pca.csv")
head(pca)
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$sample.id
res.pca <- PCA(pca2, graph = FALSE)
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
res.hcpc <- HCPC(res.pca, graph = FALSE)


f<- fviz_dend(res.hcpc, 
              cex = 0.06,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 1.0,     # Augment the room for labels
              labels_cols = TRUE, pca$Type)
#f + scale_color_manual(values=met.brewer("Hokusai1", 8))

f+ scale_color_manual(values = c( "purple",  "red",  "blue", "yellow1", "green", "yellow2", "grey"))


#########################################
#FIGURE 1B
#########################################
#Plot 2D PCA 
library("ggrepel")
#pca <- read.csv("barley_pca.csv")
pca <- read.csv("barley_pca_copy.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p

p<-ggplot(data=pca, aes(x=EV1, y=EV2,color=HCPC_group))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (9.84%)") + 
  ylab("PC 2 (6.96%)")
p + scale_color_manual(values = c( "red",  "yellow",  "green", "blue", "purple"))


#########################################
#FIGURE 1C
#########################################
#MAP
##December 5th
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
setwd("~/R/barley_collab/parental_selection/cross_validation/")
pca <- read.csv("barley_pca_copy.csv")
# Assuming you have a world map dataset named 'world'
world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()

# Set the limits for longitude and latitude to focus on Europe, Asia, and Africa
europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-30, 160), ylim = c(-40, 80))

# Add your points or other layers as needed
final_plot <- europe_asia_africa_plot +
  geom_point(data = pca, 
             aes(x = Longitude, y = Latitude, colour = HCPC_group), 
             pch = 21, size = 1, alpha = I(0.7)) 


final_plot <- europe_asia_africa_plot +
  geom_point(data = pca, 
             aes(x = Longitude, y = Latitude, colour = HCPC_group), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) 
# Display the final plot
final_plot



#########################################
#FIGURE 2
#########################################
#ENV PCA BIPLOT
library("FactoMineR")
library("factoextra")

# Extract the environmental variables
setwd("~/R/barley_collab/parental_selection/cross_validation/")
data <- read.csv("n784_climate_data_n7.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)
# Create a biplot
fviz_pca_biplot(pca_result, geom = "point", col.ind = "black", col.var = "blue", addEllipses = TRUE)
####
# Assuming 'data' is your original data
data <- read.csv("n784_climate_data_n7.csv")
env_variables <- data[, 9:15]

# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)

# Combine PCA scores with grouping variable
pca_data <- data.frame(pca_result$x, HCPC_group = data$HCPC_group)

# Create a biplot with color by 'HCPC_group'
p<- fviz_pca_biplot(pca_result, geom = "point", col.ind = "black", col.var = "blue", addEllipses = FALSE, habillage = pca_data$HCPC_group) +
  ggtitle("PCA Biplot with Color by HCPC_group")

# Adjust x and y axis limits
p1<- p + xlim(-8, 8) + ylim(-8, 8)

#ggsave("pca_biplot_Figure2F.pdf", plot = p1, width = 8, height = 6)


#Group1
data <- read.csv("biodata_HCPC_Group_1.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)
# Create a biplot
biplot<- fviz_pca_biplot(pca_result, geom = "point", col.ind = "black", col.var = "blue", addEllipses = TRUE)
biplot + ggtitle("PCA Biplot - HCPC Group 1")

p2<- biplot + ggtitle("PCA Biplot - HCPC Group 1") + xlim(-8, 8) + ylim(-8, 8)
p2
#ggsave("pca_biplot_Figure2A.pdf", plot = p2, width = 8, height = 6)


#Group2
data <- read.csv("biodata_HCPC_Group_2.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)
# Create a biplot
biplot<- fviz_pca_biplot(pca_result, geom = "point", col.ind = "black", col.var = "blue", addEllipses = TRUE)
biplot + ggtitle("PCA Biplot - HCPC Group 2")


p3<- biplot + ggtitle("PCA Biplot - HCPC Group 2") + xlim(-8, 8) + ylim(-8, 8)
p3
#ggsave("pca_biplot_Figure2B.pdf", plot = p3, width = 8, height = 6)


#Group3
data <- read.csv("biodata_HCPC_Group_3.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)
# Create a biplot
biplot<- fviz_pca_biplot(pca_result, geom = "point", col.ind = "black", col.var = "blue", addEllipses = TRUE)
biplot + ggtitle("PCA Biplot - HCPC Group 3")


p4<- biplot + ggtitle("PCA Biplot - HCPC Group 3") + xlim(-8, 8) + ylim(-8, 8)
p4
#ggsave("pca_biplot_Figure2C.pdf", plot = p4, width = 8, height = 6)


#Group4
data <- read.csv("biodata_HCPC_Group_4.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)
# Create a biplot
biplot<- fviz_pca_biplot(pca_result, geom = "point", col.ind = "black", col.var = "blue", addEllipses = TRUE)
biplot + ggtitle("PCA Biplot - HCPC Group 4")

p5<- biplot + ggtitle("PCA Biplot - HCPC Group 4") + xlim(-8, 8) + ylim(-8, 8)
p5
#ggsave("pca_biplot_Figure2D.pdf", plot = p5, width = 8, height = 6)


#Group5
data <- read.csv("biodata_HCPC_Group_5.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)
# Create a biplot
biplot<- fviz_pca_biplot(pca_result, geom = "point", col.ind = "black", col.var = "blue", addEllipses = TRUE)
biplot + ggtitle("PCA Biplot - HCPC Group 5")

p6<- biplot + ggtitle("PCA Biplot - HCPC Group 5") + xlim(-8, 8) + ylim(-8, 8)
p6
#ggsave("pca_biplot_Figure2E.pdf", plot = p6, width = 8, height = 6)



#########################################
#FIGURE 3
#########################################
#Overlaps of SDM
#Plotting the 5x population SDMs

library(raster)
library(rworldmap)
# Load your rasters
r1 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop1.tif")
r2 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop2.tif")
r3 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop3.tif")
r4 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop4.tif")
r5 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop5.tif")

# Create a mask for values greater than 0.2
s1 <- calc(r1, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )
s2 <- calc(r2, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )
s3 <- calc(r3, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )
s4 <- calc(r4, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )
s5 <- calc(r5, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )

# Prepare a map to plot on
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster Datasets")

# Define colors and add transparency
#colors <- c("#00FF7F44", "#1E90FF44", "#B2222244")  # RGBA values

# Define unique colors for each raster
#colors <- c('red', 'yellow', 'green', 'blue', 'purple')

colors <- c(adjustcolor('red', alpha=0.5), 
            adjustcolor('yellow', alpha=0.5),
            adjustcolor('green', alpha=0.5),
            adjustcolor('blue', alpha=0.2),
            adjustcolor('purple', alpha=0.2))

# Plot the masked rasters on the world map
plot(s1, col=colors[1], legend =F, add=T)
plot(s2, col=colors[2], legend =F, add=T)
plot(s3, col=colors[3], legend =F, add=T)
plot(s4, col=colors[4], legend =F, add=T)
plot(s5, col=colors[5], legend =F, add=T)


# Add a custom legend
legend("topright", legend=c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"), fill=colors, bty="n")


#########################################
#FIGURE 4
#########################################
#boxplots
setwd("~/R/barley_collab/parental_selection/cross_validation/")
data2 <-read.csv('rrblup_GEBV_data_barley_plots_copy.csv')

#bio1
ggplot(data2, aes(x = HCPC_group, y = bio1, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio1 - Annual Mean Temperature",
       x = "HCPC Group",
       y = "GEAV bio4") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 3
ggplot(data2, aes(x = HCPC_group, y = bio3, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio3 - Isothermality (BIO2/BIO7) ",
       x = "HCPC Group",
       y = "GEAV bio3") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio4
ggplot(data2, aes(x = HCPC_group, y = bio4, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio4 - temperature seasonality ",
       x = "HCPC Group",
       y = "GEAV bio4") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 6
ggplot(data2, aes(x = HCPC_group, y = bio6, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio6 - Min Temperature of Coldest Month",
       x = "HCPC Group",
       y = "GEAV bio6") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 11
ggplot(data2, aes(x = HCPC_group, y = bio11, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio11 - Mean Temperature of Coldest Quarter",
       x = "HCPC Group",
       y = "GEAV bio11") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 14
ggplot(data2, aes(x = HCPC_group, y = bio14, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio14 - Precipitation of Driest Month",
       x = "HCPC Group",
       y = "GEAV bio14") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio17
ggplot(data2, aes(x = HCPC_group, y = bio17, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio17 - Precipitation of Driest Quarter",
       x = "HCPC Group",
       y = "GEAV bio17") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")



#########################################
#FIGURE 4B
#########################################
#Core collection identification from datasets
#########################
library(vcfR)
library(gdsfmt)
library(SNPRelate)
library(gwscaR)
library(gt)
library(snpsel)
library(ggplot2)
library(poppr)
library(reshape2)
library(ggplot2)
library(dartR)
library(adegenet)
library(poppr)
library(ape)
library(RColorBrewer)
library(adegenet)
library(rCNV)

#read in vcf

setwd("~/R/barley_collab/parental_selection/")
vcf <- read.vcfR("Supplemental_dataset_1.vcf")

####convert genlight object
x <- vcfR2genlight(vcf)

#create distance matrix
x.dist <- dist(x)
ploidy(x) <- 2
#tree <- aboot(gl.rubi, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

#change to genid
x<- vcfR2genind(vcf)
#set ploidy
ploidy(x) <- 2
xmlg<-mlg(x, quiet = FALSE)

#genetic distance
x.dist <- poppr::bitwise.dist(x)
hist(x.dist)
heatmap(as.matrix(x.dist)) 

#write.csv(x.dist,"lei_vcf1_distance_matrix.csv")

# Ward Hierarchical Clustering
fit <- hclust(x.dist, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=11) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=11, border="red")


#write.csv(as.matrix(x.dist), "distance_sug.csv")
thresholds<-mlg.filter(x, threshold = 0.05, distance = "bitwise.dist", missing="mean",threads = 1L)
cutoff_predictor(thresholds, fraction = 0.05)

pthresh  <- filter_stats(x, distance = "bitwise.dist",
                         plot = TRUE, stats = "THRESHOLD", threads = 1L)
cutoff_predictor(pthresh$farthest)

# prediction for all algorithms
sapply(pthresh, cutoff_predictor) #take note of the farthest number for later use
tab<-mlg.table(x)
diversity_boot(tab, 10L)
diversity_ci(tab, 10L, rarefy = F, raw = FALSE)
diversity_stats(tab)

#The threshold used to define MLGs was predicted using cutoff_predictor with the “farthest” algorithm. 
cc_test <- clonecorrect(x)
cc_test
nInd(x) # 1903
# How many after we clone corrected for country, year, and month?
nInd(cc_test) # 1190

#write.csv(x$tab, "test_1488.csv")

sbs<-as.genclone(x)
mll(sbs) 
mlg.filter(sbs, threshold = 0.05, distance = "bitwise.dist", missing="mean",threads = 1L)
raw_dist <- function(sbs){
  dist(genind2df(sbs, usepop = FALSE))
}
(xdis <- raw_dist(sbs))

mlg.filter(sbs)<- 0.01380191  #Barley

dups_sug = mlg.id(sbs)

for (i in dups_sug){ # for each element in the list object
  if (length(dups_sug[i]) > 1){ # if the length is greater than 1
    print(i) # print individuals that are duplicates
  }
}


library(corehunter)
# precomputed distance matrix
my.data <- distances(as.matrix(x.dist))
#core <- sampleCore(my.data, size = 25)
core <- sampleCore(my.data, size = 100)
#core <- sampleCore(my.data, size = 10)
core

write.csv(core, 'barley_core_lei.csv')


#NEW CORE 100
#boxplots
setwd("~/R/barley_collab/parental_selection/new_core_100/")
#data2 <-read.csv('rrblup_GEBV_data_barley_plots_copy.csv')
data2 <-read.csv('rrblup_GEBV_data_barley_100_abbreviatedbio.csv')

#bio1
ggplot(data2, aes(x = HCPC_group, y = bio1, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio1 - Annual Mean Temperature",
       x = "HCPC Group",
       y = "GEAV bio4") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 3
ggplot(data2, aes(x = HCPC_group, y = bio3, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio3 - Isothermality (BIO2/BIO7) ",
       x = "HCPC Group",
       y = "GEAV bio3") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio4
ggplot(data2, aes(x = HCPC_group, y = bio4, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio4 - temperature seasonality ",
       x = "HCPC Group",
       y = "GEAV bio4") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 6
ggplot(data2, aes(x = HCPC_group, y = bio6, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio6 - Min Temperature of Coldest Month",
       x = "HCPC Group",
       y = "GEAV bio6") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 11
ggplot(data2, aes(x = HCPC_group, y = bio11, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio11 - Mean Temperature of Coldest Quarter",
       x = "HCPC Group",
       y = "GEAV bio11") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 14
ggplot(data2, aes(x = HCPC_group, y = bio14, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio14 - Precipitation of Driest Month",
       x = "HCPC Group",
       y = "GEAV bio14") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio17
ggplot(data2, aes(x = HCPC_group, y = bio17, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio17 - Precipitation of Driest Quarter",
       x = "HCPC Group",
       y = "GEAV bio17") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")


#########################################
#FIGURE 5
#########################################
#Population level chromosome painting
#Example of run through for Group 2
#########
#Group 2
#########
# Load necessary libraries
library(dplyr)
library(circlize)

#Working directory (ON HPC)
setwd("/home/ahmccorm/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_SDM/chromosome_painting/sep_chroms/group2/")
dat <- read.csv("Barley_marker_effects_chromosome_map_expanded.csv", header = TRUE)
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

# Unique chromosomes
all_chroms <- unique(dat$chrom)
all_chroms <- all_chroms[all_chroms != "UN"]

# Initialize list to store outputs
all_chrom_data <- list()

# Loop through each chromosome
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[2]  #variable of interest
  i <- comps  # Direct assignment
  
  # Filter non-chromosome columns and convert them to numeric
  mark <- mat_chrom[, colnames(mat) != 'chrom', drop = FALSE]
  mark <- as.data.frame(lapply(mark, as.numeric))
  
  # Check dimensions
  cat("Dimensions of mark:", dim(mark), "\n")
  cat("Dimensions of dat_chrom component:", dim(dat_chrom[, i, drop = FALSE]), "\n")
  
  # Broadcast multiplication
  if (nrow(mark) == nrow(dat_chrom[, i, drop = FALSE])) {  # Ensure row counts match
    # Replicate the column to match the number of columns in mark
    multiplier <- dat_chrom[, i, drop = FALSE]
    mark <- mark * rep(multiplier, each = ncol(mark))
  } else {
    cat("Dimension mismatch error: cannot multiply\n")
    next  # Skip this iteration if dimensions do not match
  }
  
  # Store results
  mark$chrom <- mat_chrom$chrom
  mark$position <- dat_chrom$pos
  mark$rs <- dat_chrom$rs
  mark$var <- i
  
  # Prepare filename and save data frame
  filename <- paste(i, "Marker_Effects_chrom", chrom_num, sep = "_")
  assign(filename, mark)
  all_chrom_data[[filename]] <- mark  # Store in a list for easier access
}

barley_pca <- read.csv("barley_pca_copy.csv")
individuals_in_groups <- split(barley_pca$sample.id, barley_pca$HCPC_group)
for(i in seq_along(individuals_in_groups)) {
  assign(paste0("Group", i), individuals_in_groups[[i]])
}
bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  select(position, all_of(Group2)) ########CHANGE ACCORDINGLY

bio1.chrom2.plot <- bio1_Marker_Effects_chrom_2H %>%
  select(position, all_of(Group2)) ########CHANGE ACCORDINGLY

bio1.chrom3.plot <- bio1_Marker_Effects_chrom_3H %>%
  select(position, all_of(Group2)) ########CHANGE ACCORDINGLY

bio1.chrom4.plot <- bio1_Marker_Effects_chrom_4H %>%
  select(position, all_of(Group2)) ########CHANGE ACCORDINGLY

bio1.chrom5.plot <- bio1_Marker_Effects_chrom_5H %>%
  select(position, all_of(Group2)) ########CHANGE ACCORDINGLY

bio1.chrom6.plot <- bio1_Marker_Effects_chrom_6H %>%
  select(position, all_of(Group2)) ########CHANGE ACCORDINGLY

bio1.chrom7.plot <- bio1_Marker_Effects_chrom_7H %>%
  select(position, all_of(Group2)) ########CHANGE ACCORDINGLY

col_fun1 = colorRamp2(c(-0.12, -0.06, 0, 0.06, 0.12), c("#701130", "#e06e85", "white","#5d94cb", "#163670")) #chrom1
circos.clear()
circos.par(start.degree = 65, gap.after = 45)
circos.heatmap(bio1.chrom1.plot[,2:90], col = col_fun1, track.height = 0.8, bg.border = "black", bg.lwd = 1, bg.lty = 1)

circos.clear()
circos.par(start.degree = 65, gap.after = 45)
circos.heatmap(bio1.chrom2.plot[,2:90], col = col_fun1, track.height = 0.8, bg.border = "black", bg.lwd = 1, bg.lty = 1)

circos.clear()
circos.par(start.degree = 65, gap.after = 45)
circos.heatmap(bio1.chrom3.plot[,2:90], col = col_fun1, track.height = 0.8, bg.border = "black", bg.lwd = 1, bg.lty = 1)

circos.clear()
circos.par(start.degree = 65, gap.after = 45)
circos.heatmap(bio1.chrom4.plot[,2:90], col = col_fun1, track.height = 0.8, bg.border = "black", bg.lwd = 1, bg.lty = 1)

circos.clear()
circos.par(start.degree = 65, gap.after = 45)
circos.heatmap(bio1.chrom5.plot[,2:90], col = col_fun1, track.height = 0.8, bg.border = "black", bg.lwd = 1, bg.lty = 1)

circos.clear()
circos.par(start.degree = 65, gap.after = 45)
circos.heatmap(bio1.chrom6.plot[,2:90], col = col_fun1, track.height = 0.8, bg.border = "black", bg.lwd = 1, bg.lty = 1)

circos.clear()
circos.par(start.degree = 65, gap.after = 45)
circos.heatmap(bio1.chrom7.plot[,2:90], col = col_fun1, track.height = 0.8, bg.border = "black", bg.lwd = 1, bg.lty = 1)

#
#
#
#
#########################################
#SUPPLEMENTAL FIGURES
#########################################

#########################################
#SUPPLEMENTAL FIGURE 1
#########################################
#ON HPC
########################################
###Run all k-fold cross-validations       
##########################################
library(rrBLUP)
library(hibayes)
library(dplyr)

setwd("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_2/")
# load all functions from separate R file
source("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_2/xval_kfold_functions.R") 

gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = T)
envdat <- read.csv('Supplemental_dataset_5worldclim_all.csv', head = T) # full environmental dataset
trainingset <- read.csv("Supplemental_dataset_5_worldclim_training100.csv", head = T)

gd2 <- gd1[gd1$taxa %in% envdat$Taxa,]
row.names(gd2) <- gd2$taxa # set taxa as rownames
g.in <- gd2[,-1]
g.in <- as.matrix(g.in)

### Load Environmental Data ###
row.names(envdat) <- envdat$Taxa # set gen_id as rownames
row.names(trainingset) <- trainingset$Taxa

y.trainset <- trainingset[,c(1,5:(ncol(trainingset)))] # select unique identifier and environmental data only 
y.in <- envdat[,c(1,4:ncol(envdat))] 

################################################################################
## Run functions

############
### RR-BLUP
############
y.trainset.rr <- trainingset[,c(6:(ncol(trainingset)))]
y.in.rr <- envdat[,5:ncol(envdat)]

y.trainset.mat <- as.matrix(y.trainset.rr)
y.in.mat <- as.matrix(y.in.rr)

xval_k10_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 10, reps = 50)
saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10.RData")

############
# Gaussian Kernel
############
K <- A.mat(g.in)
k_dist <- dist(K) # Calculate Relationship Matrix

xval_k10_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_GAUSS, "xval_GAUSS_kfold_10.RData")

############
### Exponential Kernel
############
K.Exp=Kernel_computation(X=g.in, name="exponential", degree=NULL, nL=NULL)
row.names(K.Exp) <- rownames(g.in)
colnames(K.Exp) <- rownames(g.in)
exp_dist <- dist(K.Exp) # Calculate Relationship Matrix\

xval_k10_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_EXP, "xval_EXP_kfold_10.RData")

#write.table(exp_dist, "distance_matrix_exponential_kernel.txt")
#write.table(k_dist, "distance_matrix_gaussian_kernel.txt")

############
#BayesCPi
############
xval_k10_BayesCpi <- k.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 10, reps = 50, niter=3000,nburn=1200)
saveRDS(xval_k10_BayesCpi, "xval_BayesCpi_kfold_10.RData")

############
#BayesLASSO
############
#xval_k10_BayesLASSO <- k.xval.BayesLASSO(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 10, reps = 50)
#saveRDS(xval_k10_BayesLASSO, "xval_BayesLASSO_kfold_10.RData")


#############################################################
##### Cross Validation - Data Organization and Plotting #####
#############################################################

library(ggplot2)
setwd("~/R/barley_collab/parental_selection/cross_validation/")
#############################################################################
#### RR-BLUP - K-Fold Cross-Validation ####

rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)

#############################################################################
# Gaussian Kernel - K-fold
gauss_kfold_10 <- readRDS('xval_GAUSS_kfold_10.RData')
gauss_kfold_10 <- gauss_kfold_10$xval.result
gauss_kfold_10$r.mean <- as.numeric(gauss_kfold_10$r.mean)

################################################################################
# Exponential Kernel - K-fold
EXP_kfold_10 <- readRDS('xval_EXP_kfold_10.RData')
EXP_kfold_10 <- EXP_kfold_10$xval.result
EXP_kfold_10$r.mean <- as.numeric(EXP_kfold_10$r.mean)

# BayesCpi - fractional
#bayescpi_frac90 <- readRDS("xval_BayesCpi_frac90.RData")
#bayescpi_frac90 <- bayescpi_frac90$xval.result
#bayescpi_frac90$r.mean <- as.numeric(bayescpi_frac90$r.mean)

# BayesCpi - K-fold
bayescpi_kfold_10 <- readRDS('xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10 <- bayescpi_kfold_10$xval.result
bayescpi_kfold_10$r.mean <- as.numeric(bayescpi_kfold_10$r.mean)

################################################################################
# Organize all dataframes for merging

## Rename model names
#rrblup_frac90$model <- "rrBLUP"
rrblup_kfold10$model <- "rrBLUP"

#gauss_frac90$model <- "Gaussian Kernel"
gauss_kfold_10$model <- "Gaussian Kernel"

#EXP_frac90$model <- "Exponential Kernel"
EXP_kfold_10$model <- "Exponential Kernel"

#bayescpi_frac90$model <- "BayesCpi"
bayescpi_kfold_10$model <- "BayesCpi"

## Input xval type
#rrblup_frac90$xval <- "Fractional 90"
rrblup_kfold10$xval <- "Ten-Fold"

#gauss_frac90$xval <- "Fractional 90"
gauss_kfold_10$xval <- "Ten-Fold"

#EXP_frac90$xval <- "Fractional 90"
EXP_kfold_10$xval <- "Ten-Fold"

#bayescpi_frac90$xval <- "Fractional 90"
bayescpi_kfold_10$xval <- "Ten-Fold"

#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)
model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)

model_list1 <- lapply(model_list, na.omit)
all_models <- do.call("rbind", model_list1)
all_models$r.sd <- as.numeric(all_models$r.sd)
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))


all_bio <- all_models[all_models$trait %in% c('bio1', 'bio2', 'bio3', 'bio4', 'bio5','bio6',
                                              'bio7','bio8','bio9','bio10','bio11',
                                              'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19'),]

#################################################################################################################
#SUP FIGURE 1A
ggplot(all_bio, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1)

#################################################################################################################
###############################################
#new core 100 
###############################################
#ON HPC
########################################
###Run all k-fold cross-validations       
##########################################
library(rrBLUP)
library(hibayes)
library(dplyr)

setwd("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_2/")
source("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_2/xval_kfold_functions.R") # load all functions


gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = T)
envdat <- read.csv('Supplemental_dataset_5worldclim_all.csv', head = T) # full environmental dataset
trainingset <- read.csv("Supplemental_dataset_5_worldclim_training100.csv", head = T)


gd2 <- gd1[gd1$taxa %in% envdat$Taxa,]
row.names(gd2) <- gd2$taxa # set taxa as rownames
g.in <- gd2[,-1]
g.in <- as.matrix(g.in)

### Load Environmental Data ###
row.names(envdat) <- envdat$Taxa # set gen_id as rownames
row.names(trainingset) <- trainingset$Taxa

y.trainset <- trainingset[,c(1,5:(ncol(trainingset)))] # select unique identifier and environmental data only 
y.in <- envdat[,c(1,4:ncol(envdat))] 

################################################################################
## Run functions

############
### RR-BLUP
############
y.trainset.rr <- trainingset[,c(6:(ncol(trainingset)))]
y.in.rr <- envdat[,5:ncol(envdat)]

y.trainset.mat <- as.matrix(y.trainset.rr)
y.in.mat <- as.matrix(y.in.rr)

xval_k10_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 10, reps = 50)
saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10.RData")

############
# Gaussian Kernel
############
K <- A.mat(g.in)
k_dist <- dist(K) # Calculate Relationship Matrix

xval_k10_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_GAUSS, "xval_GAUSS_kfold_10.RData")

############
### Exponential Kernel
############
K.Exp=Kernel_computation(X=g.in, name="exponential", degree=NULL, nL=NULL)
row.names(K.Exp) <- rownames(g.in)
colnames(K.Exp) <- rownames(g.in)
exp_dist <- dist(K.Exp) # Calculate Relationship Matrix\

xval_k10_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_EXP, "xval_EXP_kfold_10.RData")

############
#BayesCPi
############
xval_k10_BayesCpi <- k.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 10, reps = 50, niter=3000,nburn=1200)
saveRDS(xval_k10_BayesCpi, "xval_BayesCpi_kfold_10.RData")


library(ggplot2)
setwd("~/R/barley_collab/parental_selection/new_core_100/")
#############################################################################
#### RR-BLUP - K-Fold Cross-Validation ####

rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)

#############################################################################
# Gaussian Kernel - K-fold
gauss_kfold_10 <- readRDS('xval_GAUSS_kfold_10.RData')
gauss_kfold_10 <- gauss_kfold_10$xval.result
gauss_kfold_10$r.mean <- as.numeric(gauss_kfold_10$r.mean)

################################################################################
# Exponential Kernel - K-fold
EXP_kfold_10 <- readRDS('xval_EXP_kfold_10.RData')
EXP_kfold_10 <- EXP_kfold_10$xval.result
EXP_kfold_10$r.mean <- as.numeric(EXP_kfold_10$r.mean)

# BayesCpi - fractional
#bayescpi_frac90 <- readRDS("xval_BayesCpi_frac90.RData")
#bayescpi_frac90 <- bayescpi_frac90$xval.result
#bayescpi_frac90$r.mean <- as.numeric(bayescpi_frac90$r.mean)

# BayesCpi - K-fold
bayescpi_kfold_10 <- readRDS('xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10 <- bayescpi_kfold_10$xval.result
bayescpi_kfold_10$r.mean <- as.numeric(bayescpi_kfold_10$r.mean)

################################################################################
# Organize all dataframes for merging

## Rename model names
#rrblup_frac90$model <- "rrBLUP"
rrblup_kfold10$model <- "rrBLUP"

#gauss_frac90$model <- "Gaussian Kernel"
gauss_kfold_10$model <- "Gaussian Kernel"

#EXP_frac90$model <- "Exponential Kernel"
EXP_kfold_10$model <- "Exponential Kernel"

#bayescpi_frac90$model <- "BayesCpi"
bayescpi_kfold_10$model <- "BayesCpi"

## Input xval type
#rrblup_frac90$xval <- "Fractional 90"
rrblup_kfold10$xval <- "Ten-Fold"

#gauss_frac90$xval <- "Fractional 90"
gauss_kfold_10$xval <- "Ten-Fold"

#EXP_frac90$xval <- "Fractional 90"
EXP_kfold_10$xval <- "Ten-Fold"

#bayescpi_frac90$xval <- "Fractional 90"
bayescpi_kfold_10$xval <- "Ten-Fold"

#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)
#model_list <- list(rrblup_kfold10)
model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)

model_list1 <- lapply(model_list, na.omit)
all_models <- do.call("rbind", model_list1)
all_models$r.sd <- as.numeric(all_models$r.sd)
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))


all_bio <- all_models[all_models$trait %in% c('bio1', 'bio2', 'bio3', 'bio4', 'bio5','bio6',
                                              'bio7','bio8','bio9','bio10','bio11',
                                              'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19'),]

#####################################################################################################################################
##  X axis by Model Type
#SUP FIGURE 1B
ggplot(all_bio, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 0.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ylim(0, 1)

#IS PLOT FOR SUP FIG 7
#########################################
#comparing the 31 and 100 core predictions
#########################################
library(dplyr)
library(ggplot2)

# Set working directory to where the results are saved
setwd("~/R/barley_collab/parental_selection/core_comparisons/")

# Load results for Core 31
rrblup_kfold10_31 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/31/xval_rrblup_kfold_10.RData')
rrblup_kfold10_31 <- rrblup_kfold10_31$xval.result
rrblup_kfold10_31$r.mean <- as.numeric(rrblup_kfold10_31$r.mean)
rrblup_kfold10_31$r.sd <- as.numeric(rrblup_kfold10_31$r.sd)
rrblup_kfold10_31$core_size <- "Core 31"
rrblup_kfold10_31$model <- "RR-BLUP"

gauss_kfold_10_31 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/31/xval_GAUSS_kfold_10.RData')
gauss_kfold_10_31 <- gauss_kfold_10_31$xval.result
gauss_kfold_10_31$r.mean <- as.numeric(gauss_kfold_10_31$r.mean)
gauss_kfold_10_31$r.sd <- as.numeric(gauss_kfold_10_31$r.sd)
gauss_kfold_10_31$core_size <- "Core 31"
gauss_kfold_10_31$model <- "Gaussian Kernel"

EXP_kfold_10_31 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/31/xval_EXP_kfold_10.RData')
EXP_kfold_10_31 <- EXP_kfold_10_31$xval.result
EXP_kfold_10_31$r.mean <- as.numeric(EXP_kfold_10_31$r.mean)
EXP_kfold_10_31$r.sd <- as.numeric(EXP_kfold_10_31$r.sd)
EXP_kfold_10_31$core_size <- "Core 31"
EXP_kfold_10_31$model <- "Exponential Kernel"

bayescpi_kfold_10_31 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/31/xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10_31 <- bayescpi_kfold_10_31$xval.result
bayescpi_kfold_10_31$r.mean <- as.numeric(bayescpi_kfold_10_31$r.mean)
bayescpi_kfold_10_31$r.sd <- as.numeric(bayescpi_kfold_10_31$r.sd)
bayescpi_kfold_10_31$core_size <- "Core 31"
bayescpi_kfold_10_31$model <- "BayesCpi"

# Load results for Core 100
rrblup_kfold10_100 <- readRDS("~/R/barley_collab/parental_selection/core_comparisons/100/xval_rrblup_kfold_10.RData")
rrblup_kfold10_100 <- rrblup_kfold10_100$xval.result
rrblup_kfold10_100$r.mean <- as.numeric(rrblup_kfold10_100$r.mean)
rrblup_kfold10_100$r.sd <- as.numeric(rrblup_kfold10_100$r.sd)
rrblup_kfold10_100$core_size <- "Core 100"
rrblup_kfold10_100$model <- "RR-BLUP"

gauss_kfold_10_100 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/100/xval_GAUSS_kfold_10.RData')
gauss_kfold_10_100 <- gauss_kfold_10_100$xval.result
gauss_kfold_10_100$r.mean <- as.numeric(gauss_kfold_10_100$r.mean)
gauss_kfold_10_100$r.sd <- as.numeric(gauss_kfold_10_100$r.sd)
gauss_kfold_10_100$core_size <- "Core 100"
gauss_kfold_10_100$model <- "Gaussian Kernel"

EXP_kfold_10_100 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/100/xval_EXP_kfold_10.RData')
EXP_kfold_10_100 <- EXP_kfold_10_100$xval.result
EXP_kfold_10_100$r.mean <- as.numeric(EXP_kfold_10_100$r.mean)
EXP_kfold_10_100$r.sd <- as.numeric(EXP_kfold_10_100$r.sd)
EXP_kfold_10_100$core_size <- "Core 100"
EXP_kfold_10_100$model <- "Exponential Kernel"

bayescpi_kfold_10_100 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/100/xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10_100 <- bayescpi_kfold_10_100$xval.result
bayescpi_kfold_10_100$r.mean <- as.numeric(bayescpi_kfold_10_100$r.mean)
bayescpi_kfold_10_100$r.sd <- as.numeric(bayescpi_kfold_10_100$r.sd)
bayescpi_kfold_10_100$core_size <- "Core 100"
bayescpi_kfold_10_100$model <- "BayesCpi"

# Combine results into a single dataframe
combined_results <- bind_rows(
  rrblup_kfold10_31,
  gauss_kfold_10_31,
  EXP_kfold_10_31,
  bayescpi_kfold_10_31,
  rrblup_kfold10_100,
  gauss_kfold_10_100,
  EXP_kfold_10_100,
  bayescpi_kfold_10_100
)

# Organize dataframe for plotting
combined_results$model <- factor(combined_results$model, levels = c("RR-BLUP", "Gaussian Kernel", "Exponential Kernel", "BayesCpi"))
combined_results$core_size <- factor(combined_results$core_size, levels = c("Core 31", "Core 100"))

# Filter the relevant traits
all_bio <- combined_results %>% filter(trait %in% c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19'))

# Ensure all r.sd are numeric
all_bio$r.sd <- as.numeric(all_bio$r.sd)

# Plot the results
ggplot(all_bio, aes(y = r.mean, x = model, color = core_size)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(color = "Core Size", y = "Prediction Accuracy", x = "Model")



#########################################
#SUPPLEMENTAL FIGURE 2
#########################################
library(dismo)
library(raster)
library(maptools)
library(rasterVis)
data(wrld_simpl)
library(readr)

setwd("/Users/annamccormick/R/barley_collab/parental_selection/SDM/")
barley_pop1 <-read.csv("pop1.csv",header=T, sep=",")
barley_pop2 <-read.csv("pop2.csv",header=T, sep=",")
barley_pop3 <-read.csv("pop3.csv",header=T, sep=",")
barley_pop4 <-read.csv("pop4.csv",header=T, sep=",")
barley_pop5 <-read.csv("pop5.csv",header=T, sep=",")

# Combine the datasets
combined_barley <- rbind(barley_pop1, barley_pop2, barley_pop3, barley_pop4, barley_pop5)

# Write the combined dataframe to a new CSV file
write.csv(combined_barley, "combined_barley.csv", row.names = FALSE)

species_distribution1 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop1.tif")
species_distribution2 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop2.tif")
species_distribution3 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop3.tif")
species_distribution4 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop4.tif")
species_distribution5 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop5.tif")

coordinates <- combined_barley[, c("lon", "lat")]

# Convert to SpatialPointsDataFrame if not already
coordinates_sp <- SpatialPoints(coordinates, proj4string = CRS(proj4string(species_distribution1)))

# Now extract suitability scores using the SpatialPointsDataFrame
scores_pop1 <- extract(species_distribution1, coordinates_sp)
scores_pop2 <- extract(species_distribution2, coordinates_sp)
scores_pop3 <- extract(species_distribution3, coordinates_sp)
scores_pop4 <- extract(species_distribution4, coordinates_sp)
scores_pop5 <- extract(species_distribution5, coordinates_sp)

# Combine the scores into a data frame (optional, depending on your needs)
suitability_scores <- data.frame(scores_pop1, scores_pop2, scores_pop3, scores_pop4, scores_pop5)

suitability_scores$id <- combined_barley$Taxa

# Add lat and lon columns back to the suitability scores data frame
suitability_scores$lat <- combined_barley$lat
suitability_scores$lon <- combined_barley$lon

# Write the suitability scores to a new CSV file
write.csv(suitability_scores, "suitability_scores.csv", row.names = FALSE)

########################################
###Run all k-fold cross-validations       
##########################################
library(rrBLUP)
library(hibayes)
library(dplyr)

setwd("/Users/annamccormick/R/barley_collab/parental_selection/cross_validation_SDMtrial/")
source("/Users/annamccormick/R/barley_collab/parental_selection/cross_validation_SDMtrial/xval_kfold_functions.R") # load all functions

gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = T)
envdat <- read.csv('suitability_scores.csv', head = T) # full environmental dataset
trainingset <- read.csv("suitability_scores_training31.csv", head = T)

gd2 <- gd1[gd1$taxa %in% envdat$Taxa,]
row.names(gd2) <- gd2$taxa # set taxa as rownames (careful with capitals!!)
g.in <- gd2[,-1]
g.in <- as.matrix(g.in)

### Load Environmental Data ###
row.names(envdat) <- envdat$Taxa # set gen_id as rownames
row.names(trainingset) <- trainingset$Taxa

y.trainset <- trainingset[,c(1,5:(ncol(trainingset)))] # select unique identifier and environmental data only 
y.in <- envdat[,c(1,4:ncol(envdat))] 

########################
#Run functions

############
# RR-BLUP
############
y.trainset.rr <- trainingset[,c(6:(ncol(trainingset)))]
y.in.rr <- envdat[,5:ncol(envdat)]

y.trainset.mat <- as.matrix(y.trainset.rr)
y.in.mat <- as.matrix(y.in.rr)

#dim(g.in)
#dim(y.in.mat)
#dim(y.trainset.mat)

xval_k10_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 10, reps = 50)
saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10.RData")

############
# Gaussian Kernel
############
K <- A.mat(g.in)
k_dist <- dist(K) # Calculate Relationship Matrix

xval_k10_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_GAUSS, "xval_GAUSS_kfold_10.RData")

############
# Exponential Kernel
############
K.Exp=Kernel_computation(X=g.in, name="exponential", degree=NULL, nL=NULL)
row.names(K.Exp) <- rownames(g.in)
colnames(K.Exp) <- rownames(g.in)
exp_dist <- dist(K.Exp) # Calculate Relationship Matrix\

xval_k10_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_EXP, "xval_EXP_kfold_10.RData")

#write.table(exp_dist, "distance_matrix_exponential_kernel.txt")
#write.table(k_dist, "distance_matrix_gaussian_kernel.txt")

############
# BayesCPi
############
xval_k10_BayesCpi <- k.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 10, reps = 50, niter=3000,nburn=1200)
saveRDS(xval_k10_BayesCpi, "xval_BayesCpi_kfold_10.RData")

###############################################
#Plotting resuts from cross validation 
###############################################

library(ggplot2)
setwd("/Users/annamccormick/R/barley_collab/parental_selection/cross_validation_SDMtrial/")
#############################################################################
# RR-BLUP - K-Fold Cross-Validation ####
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)
#############################################################################
# Gaussian Kernel - K-fold
gauss_kfold_10 <- readRDS('xval_GAUSS_kfold_10.RData')
gauss_kfold_10 <- gauss_kfold_10$xval.result
gauss_kfold_10$r.mean <- as.numeric(gauss_kfold_10$r.mean)
################################################################################
# Exponential Kernel - K-fold
EXP_kfold_10 <- readRDS('xval_EXP_kfold_10.RData')
EXP_kfold_10 <- EXP_kfold_10$xval.result
EXP_kfold_10$r.mean <- as.numeric(EXP_kfold_10$r.mean)
################################################################################
# BayesCpi - K-fold
bayescpi_kfold_10 <- readRDS('xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10 <- bayescpi_kfold_10$xval.result
bayescpi_kfold_10$r.mean <- as.numeric(bayescpi_kfold_10$r.mean)
################################################################################
# Organize all dataframes for merging

## Rename model names
rrblup_kfold10$model <- "rrBLUP"
gauss_kfold_10$model <- "Gaussian Kernel"
EXP_kfold_10$model <- "Exponential Kernel"
bayescpi_kfold_10$model <- "BayesCpi"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"
gauss_kfold_10$xval <- "Ten-Fold"
EXP_kfold_10$xval <- "Ten-Fold"
bayescpi_kfold_10$xval <- "Ten-Fold"

#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)
model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)

model_list1 <- lapply(model_list, na.omit)
all_models <- do.call("rbind", model_list1)
all_models$r.sd <- as.numeric(all_models$r.sd)
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

all_bio <- all_models[all_models$trait %in% c('scores_pop1', 'scores_pop2', 'scores_pop3', 'scores_pop4', 'scores_pop5'),]

#############################################################################################
##  X axis by Model Type
ggplot(all_bio, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 0.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ylim(0, 1)


###############################################
#Genomic selection
###############################################
library(rrBLUP)
library(dplyr)

setwd("/Users/annamccormick/R/barley_collab/parental_selection/cross_validation_SDMtrial/")
envdat <- read.csv('suitability_scores_training31.csv', head = T) # full environmental dataset
gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = T)

row.names(gd1) <- gd1$taxa
#gd2 <- gd1[gd1$taxa %in% envdat$gen_id,]
gd3 <- gd1[,-1]
#g.in <- as.matrix(gd1)
g.in <- as.matrix(gd3)

# RR-BLUP
row.names(envdat) <- envdat$Taxa
y.in.rr <- envdat[,5:ncol(envdat)]   #change where data starts? #was 5 to include core 
#if set to 6 use envdat for Core ==T

#minicore_entries <- which(y.in.rr$Core == T)
minicore_entries <- which(envdat$Core == T)
#minicore_entries <- which(y.in.rr$Core == T)
y.in.rr <- y.in.rr[minicore_entries,]
#y.in.rr <- y.in.rr[,-ncol(y.in.rr)]
y.in.rr <- y.in.rr[,-1]
y.in.mat <- as.matrix(y.in.rr)

train <- row.names(y.in.mat) # Names of training lines
g.train <- g.in[train,] # Set training genos

traits <- colnames(y.in.mat)

# Calculate marker effects
pred <- setdiff(row.names(g.in), train)
g.pred <- g.in[pred,] # Set prediction genos
y.pred <- y.in.rr[pred,]

marker.list <- list()
gebv_df <- data.frame(matrix(nrow = nrow(envdat), ncol = length(traits)))
#row.names(gebv_df) <- row.names(envdat)
colnames(gebv_df) <- traits

for(t in 1:length(traits)){
  trait <- traits[t]
  y.train <- as.matrix(y.in.mat[train,trait])
  
  ###### Run Model: RR-BLUP ######
  solve.out <- mixed.solve(y = y.train,  Z = g.train, SE = F, return.Hinv = F)
  u.hat <- solve.out$u
  
  # Calculate GEBVs and correlations
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  #outputs
  marker.list[[t]] <- u.hat
  gebv_df[,t] <- rbind(GEBV, GEBV_train)
}

row.names(gebv_df)[1:nrow(y.pred)] <- row.names(GEBV)
row.names(gebv_df)[(nrow(y.pred) + 1):(nrow(y.pred) + nrow(y.train))] <- row.names(GEBV_train)

write.csv(gebv_df, 'rrblup_GEBV_data_barley_31.csv')
saveRDS(marker.list, 'rrblup_markereffects_data_barley_31.csv')

#########################################
#SUPPLEMENTAL FIGURE 2B
#########################################
#Plotting results by population of GEBVs
#########################################

setwd("/Users/annamccormick/R/barley_collab/parental_selection/cross_validation_SDMtrial/")
gebv <-read.csv('rrblup_GEBV_data_barley_31_plots.csv') #with HCPC designations added in

# Extracting Line and GEBV columns
line_column <- gebv$line
gebv_columns <- gebv[, -1]  # Exclude the 'Line' column

# Number of traits
num_traits <- ncol(gebv_columns)

# Set the size of the plotting area
par(mfrow = c(ceiling(num_traits / 2), 2), mar = c(4, 4, 2, 1))

for (i in 1:num_traits) {
  trait_values <- gebv_columns[, i]
  bar_colors <- ifelse(trait_values < 0, "red", "blue")
  
  # Create a bar plot
  barplot(trait_values, names.arg = line_column, col = bar_colors,
          main = paste("Genomic Estimated Breeding Values for", colnames(gebv_columns)[i]),
          xlab = "Line", ylab = "GEBV")
}
#############
#SUP FIG 2B (Boxplots)
#############
data2 <-read.csv('rrblup_GEBV_data_barley_31_plots.csv')

library(ggplot2)
# Pop1
a <- ggplot(data2, aes(x = HCPC_group, y = scores_pop1, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "Population 1 raster",
       x = "HCPC Group",
       y = "GEBV") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
a
# Pop2
b <-ggplot(data2, aes(x = HCPC_group, y = scores_pop2, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "Population 2 raster",
       x = "HCPC Group",
       y = "GEBV") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
b
# Pop3
c<- ggplot(data2, aes(x = HCPC_group, y = scores_pop3, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "Population 3 raster",
       x = "HCPC Group",
       y = "GEBV") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
c
# Pop4
d<- ggplot(data2, aes(x = HCPC_group, y = scores_pop4, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "Population 4 raster",
       x = "HCPC Group",
       y = "GEBV") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
d
# Pop5
e<-ggplot(data2, aes(x = HCPC_group, y = scores_pop5, color = HCPC_group)) +
  geom_boxplot() +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "Population 5 raster",
       x = "HCPC Group",
       y = "GEBV") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
e

library(gridExtra)
grid.arrange(a, b, c, d, e, nrow=3)


#############################
#Generating chromosome map expanded
#takes marker effect_data from rrBLUP
#takes env data

#output = Barley_marker_effects_chromosome_map_expanded.csv

library(dplyr)
setwd("~/R/barley_collab/parental_selection/cross_validation/")
file <- 'Barley_hapmap_second.hmp.txt'

full <- read.delim(file, header = T)
full_df <- full[,c(1,3,4)]

#write and add in X for downstream merging
full_df

write.csv(full_df, file = "full_df_editing.csv", row.names = FALSE)
full_df<-read.csv('full_df_editing.csv')

markers <- readRDS('rrblup_markereffects_data_barley.csv')

# need to give names to list of marker effects
#envdat <- read.csv('./input_data/envdat_wildsoybean_master.csv')
envdat <- read.csv('Supplemental_dataset_5worldclim_all_chromosome_plot_master.csv')
head(envdat) 

vars <- colnames(envdat)[8:ncol(envdat)]
names(markers) <- vars
head(markers)

bio1<- markers[['bio1']]
bio3<- markers[['bio3']]
bio4 <- markers[["bio4"]]
bio6 <- markers[['bio6']]
bio11 <- markers[['bio11']]
bio14 <- markers[['bio14']]
bio17 <- markers[['bio17']]


#markers_sel <- cbind(names(markers[[1]]), ann.prc, prc.dry.m, prc.seas, ann.mean.tmp, mean.tmp.wrm.q, tmp.seas)
markers_sel <- cbind(names(markers[[1]]), bio1, bio3, bio4, bio6, bio11, bio14, bio17)

markers_sel <- as.data.frame(markers_sel)
markers_sel <- markers_sel %>% mutate_at(2:7, as.numeric)

write.csv(markers_sel, file = "markers_sel_editing.csv", row.names = FALSE)
markers_sel <- read.csv('markers_sel_editing.csv') #needed to fix the X label here 
colnames(markers_sel)[1] <- 'rs'
#there is an X in the markers in markers select have to fix
marker_chrom <- merge(markers_sel, full_df, by = 'rs')

# Check if 'rs' column exists in both data frames
'rs' %in% colnames(markers_sel)  # should return TRUE
'rs' %in% colnames(full_df)  # should return TRUE

write.csv(marker_chrom, 'Barley_marker_effects_chromosome_map_expanded.csv', row.names = F)
####################################


#########################################
# SUPPLEMENTAL FIGURE 3
#########################################
#Overlaps

#SUP FIGURE 3A
library(VennDiagram)
file_path <- "/Users/annamccormick/R/barley_collab/parental_selection/climate_line_overlaps/temp_overlap_lines.csv"
overlaps <- read.csv(file_path)

# Sample data
set1 <- overlaps$bio1
set2 <- overlaps$bio4
set3 <- overlaps$bio6
set4 <- overlaps$bio11

# Create Venn diagram
venn.plot <- venn.plot <- venn.diagram(
  x = list(Set1 = set1, Set2 = set2, Set3 = set3, Set4 = set4),
  category.names = c("Set 1- Bio1", "Set 2 - bio4 ", "Set 3 - bio6", "Set 4 - bio11"),
  filename = NULL,
  output = TRUE
)

# Plot the Venn diagram
grid.draw(venn.plot)

#SUP FIGURE 3B
library(VennDiagram)
file_path <- "/Users/annamccormick/R/barley_collab/parental_selection/climate_line_overlaps/ppt_overlap_lines.csv"
overlaps <- read.csv(file_path)
# Sample data
set1 <- overlaps$bio14
set2 <- overlaps$bio17


# Create Venn diagram
venn.plot <- venn.plot <- venn.diagram(
  x = list(Set1 = set1, Set2 = set2),
  category.names = c("Set 1 - bio14 ", "Set 2- bio17"),
  filename = NULL,
  output = TRUE
)

# Plot the Venn diagram
grid.draw(venn.plot)

#########################################
#SUPPLEMENTAL FIGURE 4
#########################################
library(dplyr)
library(circlize)

# Set working directory
#setwd("/home/ahmccorm/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_SDM/chromosome_painting/")
setwd("~/R/barley_collab/parental_selection/cross_validation/")
# Read data
dat <- read.csv("Barley_marker_effects_chromosome_map_expanded.csv", header = TRUE)
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

# Unique chromosomes
all_chroms <- unique(dat$chrom)
all_chroms <- all_chroms[all_chroms != "UN"]

# Initialize list to store outputs
all_chrom_data <- list()

# Loop through each chromosome
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[2]  # Variable of interest
  i <- comps  # Direct assignment
  
  # Filter non-chromosome columns and convert them to numeric
  mark <- mat_chrom[, colnames(mat) != 'chrom', drop = FALSE]
  mark <- as.data.frame(lapply(mark, as.numeric))
  
  # Check dimensions
  cat("Dimensions of mark:", dim(mark), "\n")
  cat("Dimensions of dat_chrom component:", dim(dat_chrom[, i, drop = FALSE]), "\n")
  
  # Broadcast multiplication
  if (nrow(mark) == nrow(dat_chrom[, i, drop = FALSE])) {  # Ensure row counts match
    # Replicate the column to match the number of columns in mark
    multiplier <- dat_chrom[, i, drop = FALSE]
    mark <- mark * rep(multiplier, each = ncol(mark))
  } else {
    cat("Dimension mismatch error: cannot multiply\n")
    next  # Skip this iteration if dimensions do not match
  }
  
  # Store results
  mark$chrom <- mat_chrom$chrom
  mark$position <- dat_chrom$pos
  mark$rs <- dat_chrom$rs
  mark$var <- i
  
  # Prepare filename and save data frame
  filename <- paste(i, "Marker_Effects_chrom", chrom_num, sep = "_")
  assign(filename, mark)
  all_chrom_data[[filename]] <- mark  # Store in a list for easier access
}


# Extract and sort the dataframe for Chromosome 1H
bio1_Marker_Effects_chrom_1H <- all_chrom_data[["bio1_Marker_Effects_chrom_1H"]]
#print(colnames(bio1_Marker_Effects_chrom_1H))

bio1_Marker_Effects_chrom_1H <- bio1_Marker_Effects_chrom_1H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_2H <- all_chrom_data[["bio1_Marker_Effects_chrom_2H"]]
#print(colnames(bio1_Marker_Effects_chrom_2H))

bio1_Marker_Effects_chrom_2H <- bio1_Marker_Effects_chrom_2H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_3H <- all_chrom_data[["bio1_Marker_Effects_chrom_3H"]]
#print(colnames(bio1_Marker_Effects_chrom_3H))

bio1_Marker_Effects_chrom_3H <- bio1_Marker_Effects_chrom_3H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_4H <- all_chrom_data[["bio1_Marker_Effects_chrom_4H"]]
#print(colnames(bio1_Marker_Effects_chrom_4H))

bio1_Marker_Effects_chrom_4H <- bio1_Marker_Effects_chrom_4H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_5H <- all_chrom_data[["bio1_Marker_Effects_chrom_5H"]]
#print(colnames(bio1_Marker_Effects_chrom_5H))

bio1_Marker_Effects_chrom_5H <- bio1_Marker_Effects_chrom_5H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_6H <- all_chrom_data[["bio1_Marker_Effects_chrom_6H"]]
#print(colnames(bio1_Marker_Effects_chrom_6H))

bio1_Marker_Effects_chrom_6H <- bio1_Marker_Effects_chrom_6H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_7H <- all_chrom_data[["bio1_Marker_Effects_chrom_7H"]]
#print(colnames(bio1_Marker_Effects_chrom_7H))

bio1_Marker_Effects_chrom_7H <- bio1_Marker_Effects_chrom_7H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

###### ALL CHROM
library(dplyr)
library(purrr)

# Assuming 'all_chrom_data' is a list with data frames for each chromosome
# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(all_chrom_data)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Calculate summary statistics for each SNP
summary_stats <- combined_data %>%
  group_by(rs) %>%
  summarise(min_effect = min(mean_effect, na.rm = TRUE),
            max_effect = max(mean_effect, na.rm = TRUE),
            median_effect = median(mean_effect, na.rm = TRUE))

# Plotting for mean effect of all SNPS
ggplot(summary_stats, aes(x = rs)) +
  geom_linerange(aes(ymin = min_effect, ymax = max_effect), color = "blue") +
  geom_point(aes(y = median_effect), color = "red", size = 2) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # Remove x-axis text for clarity
  labs(title = "Range and Median of Mean Effect per SNP",
       x = "SNP",
       y = "Mean Effect")


# Sort the data by mean_effect in descending order to get the largest effect sizes at the top
combined_sorted <- combined_data %>%
  arrange(desc(abs(mean_effect)))

# Select the top 10 SNPs with the highest mean effect
top_10_snps <- head(combined_sorted, 10)

# Select and display the relevant information for these top 10 SNPs
top_10_snps_selected <- top_10_snps %>%
  select(chrom, position, rs, var, mean_effect)

# Display the top 10 SNPs with selected information
print(top_10_snps_selected)


##############
# Allelic state of a SNP extraction

#First
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

# Assuming the rows are SNP identifiers and columns are line identifiers
if ("SCRI_RS_143091" %in% rownames(mat)) {
  allelic_states <- mat["SCRI_RS_143091", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP SCRI_RS_143091 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP SCRI_RS_143091 along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_SCRI_RS_143091.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP SCRI_RS_143091")

ggsave("SNP_Distribution_Map_1.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")

####################################
# Example data of significant SNPs
significant_snps <- data.frame(
  rs = c("SCRI_RS_143091", "SCRI_RS_140586", "SCRI_RS_176092", "11_20939",
         "SCRI_RS_127903", "SCRI_RS_230691", "11_20662", "12_31294", "SCRI_RS_123814", "SCRI_RS_236569")
)

# Adding a new column to easily merge
significant_snps$significant = TRUE

library(dplyr)

combined_data <- combined_data %>%
  left_join(significant_snps, by = "rs") %>%
  mutate(significant = ifelse(is.na(significant), FALSE, TRUE))


summary_stats <- combined_data %>%
  group_by(rs) %>%
  summarise(min_effect = min(mean_effect, na.rm = TRUE),
            max_effect = max(mean_effect, na.rm = TRUE),
            median_effect = median(mean_effect, na.rm = TRUE),
            significant = any(significant, na.rm = TRUE))  # Check if any of the grouped entries is significant

ggplot(summary_stats, aes(x = rs, ymin = min_effect, ymax = max_effect, color = significant)) +
  geom_linerange() +
  geom_point(aes(y = median_effect), size = 2) +
  scale_color_manual(values = c("blue", "red")) +  # Non-significant in blue, significant in red
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Range and Median of Mean Effect per SNP, Highlighting Significant SNPs",
       x = "SNP",
       y = "Mean Effect")

####################################

top_snps <- c("SCRI_RS_143091", "SCRI_RS_140586", "SCRI_RS_176092", "11_20939",
              "SCRI_RS_127903", "SCRI_RS_230691", "11_20662", "12_31294",
              "SCRI_RS_123814", "SCRI_RS_236569")

####################################
#Second
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")


# Assuming the rows are SNP identifiers and columns are line identifiers
if ("SCRI_RS_140586" %in% rownames(mat)) {
  allelic_states <- mat["SCRI_RS_140586", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP SCRI_RS_140586 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP SCRI_RS_143091 along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_SCRI_RS_140586.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP SCRI_RS_140586")
ggsave("SNP_Distribution_Map_2.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")

####################################
#Third - SCRI_RS_176092
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")


# Assuming the rows are SNP identifiers and columns are line identifiers
if ("SCRI_RS_176092" %in% rownames(mat)) {
  allelic_states <- mat["SCRI_RS_176092", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP SCRI_RS_176092 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP SCRI_RS_176092 along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_SCRI_RS_176092.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP SCRI_RS_176092")
ggsave("SNP_Distribution_Map_3.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")

####################################
#Fourth - 11_20939
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")


# Assuming the rows are SNP identifiers and columns are line identifiers
if ("11_20939" %in% rownames(mat)) {
  allelic_states <- mat["11_20939", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP 11_20939 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_11_20939.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP 11_20939")
ggsave("SNP_Distribution_Map_4.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")


####################################
#Fifth - SCRI_RS_127903
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")


# Assuming the rows are SNP identifiers and columns are line identifiers
if ("SCRI_RS_127903" %in% rownames(mat)) {
  allelic_states <- mat["SCRI_RS_127903", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP SCRI_RS_127903 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_SCRI_RS_127903.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP SCRI_RS_127903")
ggsave("SNP_Distribution_Map_5.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")



####################################
#Sixth
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")


# Assuming the rows are SNP identifiers and columns are line identifiers
if ("SCRI_RS_230691" %in% rownames(mat)) {
  allelic_states <- mat["SCRI_RS_230691", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP SCRI_RS_230691 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_SCRI_RS_230691.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP SCRI_RS_230691")
ggsave("SNP_Distribution_Map_6.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")


####################################
#Seventh
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")


# Assuming the rows are SNP identifiers and columns are line identifiers
if ("11_20662" %in% rownames(mat)) {
  allelic_states <- mat["11_20662", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP 11_20662 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_11_20662.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP 11_20662")
ggsave("SNP_Distribution_Map_7.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")

####################################
#Eight
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")


# Assuming the rows are SNP identifiers and columns are line identifiers
if ("12_31294" %in% rownames(mat)) {
  allelic_states <- mat["12_31294", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP 12_31294 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_12_31294.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP 12_31294")
ggsave("SNP_Distribution_Map_8.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")



####################################
#Nine
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")


# Assuming the rows are SNP identifiers and columns are line identifiers
if ("SCRI_RS_123814" %in% rownames(mat)) {
  allelic_states <- mat["SCRI_RS_123814", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP SCRI_RS_123814 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_SCRI_RS_123814.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP SCRI_RS_123814")
ggsave("SNP_Distribution_Map_9.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")

####################################
# TEN
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")


# Assuming the rows are SNP identifiers and columns are line identifiers
if ("SCRI_RS_236569" %in% rownames(mat)) {
  allelic_states <- mat["SCRI_RS_236569", , drop = FALSE]  # Extract the row for the SNP of interest
} else {
  cat("SNP SCRI_RS_236569 not found in the dataset.")
}

# Optionally convert the allelic states into a data frame for easier manipulation
allelic_states_df <- as.data.frame(t(allelic_states))

allelic_states_df$Line <- rownames(allelic_states_df)  # Set row names as a column
colnames(allelic_states_df)[1] <- "AllelicState"  # Rename the first column to 'AllelicState'

# Load the barley_pca data if not already loaded
barley_pca <- read.csv("barley_pca_copy.csv")

# Rename the 'Line' column in allelic_states_df to 'sample.id' to match the barley_pca
allelic_states_df <- allelic_states_df %>%
  dplyr::rename(sample.id = Line)

# Merge the allelic state data with the barley_pca data
barley_pca_enriched <- barley_pca %>%
  dplyr::left_join(allelic_states_df, by = "sample.id")

#arley_pca_enriched contains the allelic states for SNP along with the other data
print(barley_pca_enriched)

filtered_data <- barley_pca_enriched %>%
  dplyr::select(sample.id, Latitude, Longitude, AllelicState)

# Check the structure of the filtered data to ensure it has the correct columns
print(head(filtered_data))

#write.csv(filtered_data, "filtered_barley_data_SCRI_RS_123814.csv", row.names = FALSE)

###### MAP
library(ggplot2)
library(sf)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()


europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-10, 140), ylim = c(0, 60))


final_plot <- europe_asia_africa_plot +
  geom_point(data = filtered_data, 
             aes(x = Longitude, y = Latitude, fill = AllelicState), 
             pch = 21, size = 1, alpha = I(0.7)) + 
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "-1" = "green")) 

final_plot + labs(title = "Geographic Distribution of SNP SCRI_RS_236569")
ggsave("SNP_Distribution_Map_10.pdf", plot = final_plot, device = "pdf", width = 11, height = 8.5, units = "in")

#########################################
#SUPPLEMENTAL FIGURE 5
#########################################
# Load necessary libraries
library(dplyr)
library(circlize)

#Set working directory
setwd("~/R/barley_collab/parental_selection/cross_validation/")
dat <- read.csv("Barley_marker_effects_chromosome_map_expanded.csv", header = TRUE)
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

all_chroms <- unique(dat$chrom)
all_chroms <- all_chroms[all_chroms != "UN"]

# Initialize list to store outputs
all_chrom_data <- list()

# Loop through each group
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[2]  # Variable of interest
  i <- comps  # Direct assignment
  
  # Filter non-chromosome columns and convert them to numeric
  mark <- mat_chrom[, colnames(mat) != 'chrom', drop = FALSE]
  mark <- as.data.frame(lapply(mark, as.numeric))
  
  # Check dimensions
  cat("Dimensions of mark:", dim(mark), "\n")
  cat("Dimensions of dat_chrom component:", dim(dat_chrom[, i, drop = FALSE]), "\n")
  
  # Broadcast multiplication
  if (nrow(mark) == nrow(dat_chrom[, i, drop = FALSE])) {  # Ensure row counts match
    # Replicate the column to match the number of columns in mark
    multiplier <- dat_chrom[, i, drop = FALSE]
    mark <- mark * rep(multiplier, each = ncol(mark))
  } else {
    cat("Dimension mismatch error: cannot multiply\n")
    next  # Skip this iteration if dimensions do not match
  }
  
  # Store results
  mark$chrom <- mat_chrom$chrom
  mark$position <- dat_chrom$pos
  mark$rs <- dat_chrom$rs
  mark$var <- i
  
  # Prepare filename and save data frame
  filename <- paste(i, "Marker_Effects_chrom", chrom_num, sep = "_")
  assign(filename, mark)
  all_chrom_data[[filename]] <- mark  # Store in a list for easier access
}


# Extract and sort the dataframe for Chromosome 1H
bio1_Marker_Effects_chrom_1H <- all_chrom_data[["bio1_Marker_Effects_chrom_1H"]]
#print(colnames(bio1_Marker_Effects_chrom_1H))

bio1_Marker_Effects_chrom_1H <- bio1_Marker_Effects_chrom_1H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_2H <- all_chrom_data[["bio1_Marker_Effects_chrom_2H"]]
#print(colnames(bio1_Marker_Effects_chrom_2H))

bio1_Marker_Effects_chrom_2H <- bio1_Marker_Effects_chrom_2H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_3H <- all_chrom_data[["bio1_Marker_Effects_chrom_3H"]]
#print(colnames(bio1_Marker_Effects_chrom_3H))

bio1_Marker_Effects_chrom_3H <- bio1_Marker_Effects_chrom_3H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_4H <- all_chrom_data[["bio1_Marker_Effects_chrom_4H"]]
#print(colnames(bio1_Marker_Effects_chrom_4H))

bio1_Marker_Effects_chrom_4H <- bio1_Marker_Effects_chrom_4H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_5H <- all_chrom_data[["bio1_Marker_Effects_chrom_5H"]]
#print(colnames(bio1_Marker_Effects_chrom_5H))

bio1_Marker_Effects_chrom_5H <- bio1_Marker_Effects_chrom_5H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_6H <- all_chrom_data[["bio1_Marker_Effects_chrom_6H"]]
#print(colnames(bio1_Marker_Effects_chrom_6H))

bio1_Marker_Effects_chrom_6H <- bio1_Marker_Effects_chrom_6H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio1_Marker_Effects_chrom_7H <- all_chrom_data[["bio1_Marker_Effects_chrom_7H"]]
#print(colnames(bio1_Marker_Effects_chrom_7H))

bio1_Marker_Effects_chrom_7H <- bio1_Marker_Effects_chrom_7H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position


##############################
#Pop 1
##############################
barley_pca <- read.csv("barley_pca_copy.csv")
individuals_in_groups <- split(barley_pca$sample.id, barley_pca$HCPC_group)
for(i in seq_along(individuals_in_groups)) {
  assign(paste0("Group", i), individuals_in_groups[[i]])
}

names(individuals_in_groups)

# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group1_column_names <- c("chrom", "position","rs", individuals_in_groups[["Group 1"]])  # Assuming these are valid column names

bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  dplyr::select(all_of(group1_column_names))

bio1.chrom2.plot <- bio1_Marker_Effects_chrom_2H %>%
  dplyr::select(all_of(group1_column_names))

bio1.chrom3.plot <- bio1_Marker_Effects_chrom_3H %>%
  dplyr::select(all_of(group1_column_names))

bio1.chrom4.plot <- bio1_Marker_Effects_chrom_4H %>%
  dplyr::select(all_of(group1_column_names))

bio1.chrom5.plot <- bio1_Marker_Effects_chrom_5H %>%
  dplyr::select(all_of(group1_column_names))

bio1.chrom6.plot <- bio1_Marker_Effects_chrom_6H %>%
  dplyr::select(all_of(group1_column_names))

bio1.chrom7.plot <- bio1_Marker_Effects_chrom_7H %>%
  dplyr::select(all_of(group1_column_names))


library(dplyr)

# Merge all chromosome dataframes
bio1_all_chroms_plot <- bind_rows(
  bio1.chrom1.plot,
  bio1.chrom2.plot,
  bio1.chrom3.plot,
  bio1.chrom4.plot,
  bio1.chrom5.plot,
  bio1.chrom6.plot,
  bio1.chrom7.plot
)

# Check the combined dataframe
print(bio1_all_chroms_plot)

###### ALL CHROM
library(dplyr)
library(purrr)

# Assuming 'all_chrom_data' is a list with data frames for each chromosome
# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(bio1_all_chroms_plot)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Calculate summary statistics for each SNP
summary_stats <- combined_data %>%
  group_by(rs) %>%
  summarise(min_effect = min(mean_effect, na.rm = TRUE),
            max_effect = max(mean_effect, na.rm = TRUE),
            median_effect = median(mean_effect, na.rm = TRUE))

###################
library(dplyr)

# Assuming 'combined_data' contains the mean effects for each SNP
mean_effect <- mean(combined_data$mean_effect, na.rm = TRUE)
std_dev_effect <- sd(combined_data$mean_effect, na.rm = TRUE)

# Calculate Z-scores for each SNP
combined_data$z_score <- abs((combined_data$mean_effect - mean_effect) / std_dev_effect)

# Filter SNPs that are significant outliers, here defined as having a Z-score > 2
significant_snps_z <- combined_data %>%
  filter(z_score > 2) %>%
  arrange(desc(z_score))


# Calculate median and MAD (Median Absolute Deviation)
median_effect <- median(combined_data$mean_effect, na.rm = TRUE)
mad_effect <- mad(combined_data$mean_effect, na.rm = TRUE, constant = 1)

# Calculate robust Z-scores (modified Z-scores using MAD)
combined_data$robust_z_score <- abs((combined_data$mean_effect - median_effect) / mad_effect)

# Filter SNPs that are significant outliers, here defined as having a robust Z-score > 3
significant_snps_mad <- combined_data %>%
  filter(robust_z_score > 3) %>%
  arrange(desc(robust_z_score))

# Combine and compare significant SNPs identified by both methods
combined_significant_snps <- significant_snps_z %>%
  rename(z_score_significance = z_score) %>%
  inner_join(significant_snps_mad %>%
               rename(mad_score_significance = robust_z_score), by = c("chrom", "position", "rs", "mean_effect"))

# Print combined significant SNPs
print(combined_significant_snps)

# Plot mean effects with highlighting for outliers
library(ggplot2)
library(dplyr)

# First, add a new highlight column to the combined_data
combined_data$significant_highlight <- ifelse(combined_data$rs %in% combined_significant_snps$rs, "Significant", "Other")

# Now, create the plot
ggplot(combined_data, aes(x = position, y = mean_effect)) +
  geom_point(aes(color = significant_highlight, size = significant_highlight), alpha = 0.2) +  # Base layer for size and color
  geom_point(data = subset(combined_data, significant_highlight == "Significant"), 
             aes(color = significant_highlight, size = significant_highlight), 
             shape = 21,   # Shape 21 is a filled circle
             fill = "red", # Set fill color for significant SNPs
             stroke = 1.5, # Set border thickness for clarity
             alpha = 1) +  # Opaque for significant SNPs
  scale_color_manual(values = c("Significant" = "red", "Other" = "gray")) +
  scale_size_manual(values = c("Significant" = 0.2, "Other" = 1)) +
  facet_wrap(~ chrom, scales = "free_x") +  # Faceting by chromosome
  theme_minimal() +
  labs(title = "SNP Effects by Chromosome with Significant SNPs Highlighted",
       x = "Position on Chromosome",
       y = "Mean Effect Size") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Improve readability of x-axis labels
        legend.position = "bottom")

#print(combined_data)
##################
# Write the data frame to a CSV file
#write.csv(combined_significant_snps, "combined_significant_snps_pop1.csv", row.names = FALSE)
#write.csv(combined_significant_snps, "combined_significant_snps_pop2.csv", row.names = FALSE)
#write.csv(combined_significant_snps, "combined_significant_snps_pop3.csv", row.names = FALSE)
#write.csv(combined_significant_snps, "combined_significant_snps_pop4.csv", row.names = FALSE)
write.csv(combined_significant_snps, "combined_significant_snps_pop5.csv", row.names = FALSE)


##############################
#Run first segment then iterate through pops
##############################
#Pop 2
##############################
barley_pca <- read.csv("barley_pca_copy.csv")
individuals_in_groups <- split(barley_pca$sample.id, barley_pca$HCPC_group)
for(i in seq_along(individuals_in_groups)) {
  assign(paste0("Group", i), individuals_in_groups[[i]])
}

names(individuals_in_groups)

# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group2_column_names <- c("chrom", "position","rs", individuals_in_groups[["Group 2"]])  # Assuming these are valid column names

bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  dplyr::select(all_of(group2_column_names))

bio1.chrom2.plot <- bio1_Marker_Effects_chrom_2H %>%
  dplyr::select(all_of(group2_column_names))

bio1.chrom3.plot <- bio1_Marker_Effects_chrom_3H %>%
  dplyr::select(all_of(group2_column_names))

bio1.chrom4.plot <- bio1_Marker_Effects_chrom_4H %>%
  dplyr::select(all_of(group2_column_names))

bio1.chrom5.plot <- bio1_Marker_Effects_chrom_5H %>%
  dplyr::select(all_of(group2_column_names))

bio1.chrom6.plot <- bio1_Marker_Effects_chrom_6H %>%
  dplyr::select(all_of(group2_column_names))

bio1.chrom7.plot <- bio1_Marker_Effects_chrom_7H %>%
  dplyr::select(all_of(group2_column_names))


library(dplyr)

# Merge all chromosome dataframes
bio1_all_chroms_plot <- bind_rows(
  bio1.chrom1.plot,
  bio1.chrom2.plot,
  bio1.chrom3.plot,
  bio1.chrom4.plot,
  bio1.chrom5.plot,
  bio1.chrom6.plot,
  bio1.chrom7.plot
)

# Check the combined dataframe
print(bio1_all_chroms_plot)

###### ALL CHROM
library(dplyr)
library(purrr)

# Assuming 'all_chrom_data' is a list with data frames for each chromosome
# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(bio1_all_chroms_plot)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Calculate summary statistics for each SNP
summary_stats <- combined_data %>%
  group_by(rs) %>%
  summarise(min_effect = min(mean_effect, na.rm = TRUE),
            max_effect = max(mean_effect, na.rm = TRUE),
            median_effect = median(mean_effect, na.rm = TRUE))

library(dplyr)

# Assuming 'combined_data' contains the mean effects for each SNP
mean_effect <- mean(combined_data$mean_effect, na.rm = TRUE)
std_dev_effect <- sd(combined_data$mean_effect, na.rm = TRUE)

# Calculate Z-scores for each SNP
combined_data$z_score <- abs((combined_data$mean_effect - mean_effect) / std_dev_effect)

# Filter SNPs that are significant outliers, here defined as having a Z-score > 2
significant_snps_z <- combined_data %>%
  filter(z_score > 2) %>%
  arrange(desc(z_score))


# Calculate median and MAD (Median Absolute Deviation)
median_effect <- median(combined_data$mean_effect, na.rm = TRUE)
mad_effect <- mad(combined_data$mean_effect, na.rm = TRUE, constant = 1)

# Calculate robust Z-scores (modified Z-scores using MAD)
combined_data$robust_z_score <- abs((combined_data$mean_effect - median_effect) / mad_effect)

# Filter SNPs that are significant outliers, here defined as having a robust Z-score > 3
significant_snps_mad <- combined_data %>%
  filter(robust_z_score > 3) %>%
  arrange(desc(robust_z_score))

# Combine and compare significant SNPs identified by both methods
combined_significant_snps <- significant_snps_z %>%
  rename(z_score_significance = z_score) %>%
  inner_join(significant_snps_mad %>%
               rename(mad_score_significance = robust_z_score), by = c("chrom", "position", "rs", "mean_effect"))

# Print combined significant SNPs
print(combined_significant_snps)

###################
# First, add a new highlight column to the combined_data
combined_data$significant_highlight <- ifelse(combined_data$rs %in% combined_significant_snps$rs, "Significant", "Other")

# Now, create the plot
ggplot(combined_data, aes(x = position, y = mean_effect)) +
  geom_point(aes(color = significant_highlight, size = significant_highlight), alpha = 0.2) +  # Base layer for size and color
  geom_point(data = subset(combined_data, significant_highlight == "Significant"), 
             aes(color = significant_highlight, size = significant_highlight), 
             shape = 21,   # Shape 21 is a filled circle
             fill = "red", # Set fill color for significant SNPs
             stroke = 1.5, # Set border thickness for clarity
             alpha = 1) +  # Opaque for significant SNPs
  scale_color_manual(values = c("Significant" = "red", "Other" = "gray")) +
  scale_size_manual(values = c("Significant" = 0.2, "Other" = 1)) +
  facet_wrap(~ chrom, scales = "free_x") +  # Faceting by chromosome
  theme_minimal() +
  labs(title = "SNP Effects by Chromosome with Significant SNPs Highlighted",
       x = "Position on Chromosome",
       y = "Mean Effect Size") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Improve readability of x-axis labels
        legend.position = "bottom")


##############################
#Pop 3
##############################
barley_pca <- read.csv("barley_pca_copy.csv")
individuals_in_groups <- split(barley_pca$sample.id, barley_pca$HCPC_group)
for(i in seq_along(individuals_in_groups)) {
  assign(paste0("Group", i), individuals_in_groups[[i]])
}


names(individuals_in_groups)
# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group3_column_names <- c("chrom", "position","rs", individuals_in_groups[["Group 3"]])  # Assuming these are valid column names

bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  dplyr::select(all_of(group3_column_names))

bio1.chrom2.plot <- bio1_Marker_Effects_chrom_2H %>%
  dplyr::select(all_of(group3_column_names))

bio1.chrom3.plot <- bio1_Marker_Effects_chrom_3H %>%
  dplyr::select(all_of(group3_column_names))

bio1.chrom4.plot <- bio1_Marker_Effects_chrom_4H %>%
  dplyr::select(all_of(group3_column_names))

bio1.chrom5.plot <- bio1_Marker_Effects_chrom_5H %>%
  dplyr::select(all_of(group3_column_names))

bio1.chrom6.plot <- bio1_Marker_Effects_chrom_6H %>%
  dplyr::select(all_of(group3_column_names))

bio1.chrom7.plot <- bio1_Marker_Effects_chrom_7H %>%
  dplyr::select(all_of(group3_column_names))


library(dplyr)

# Merge all chromosome dataframes
bio1_all_chroms_plot <- bind_rows(
  bio1.chrom1.plot,
  bio1.chrom2.plot,
  bio1.chrom3.plot,
  bio1.chrom4.plot,
  bio1.chrom5.plot,
  bio1.chrom6.plot,
  bio1.chrom7.plot
)

# Check the combined dataframe
print(bio1_all_chroms_plot)

###### ALL CHROM
library(dplyr)
library(purrr)

# Assuming 'all_chrom_data' is a list with data frames for each chromosome
# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(bio1_all_chroms_plot)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Calculate summary statistics for each SNP
summary_stats <- combined_data %>%
  group_by(rs) %>%
  summarise(min_effect = min(mean_effect, na.rm = TRUE),
            max_effect = max(mean_effect, na.rm = TRUE),
            median_effect = median(mean_effect, na.rm = TRUE))

###################
library(dplyr)

# Assuming 'combined_data' contains the mean effects for each SNP
mean_effect <- mean(combined_data$mean_effect, na.rm = TRUE)
std_dev_effect <- sd(combined_data$mean_effect, na.rm = TRUE)

# Calculate Z-scores for each SNP
combined_data$z_score <- abs((combined_data$mean_effect - mean_effect) / std_dev_effect)

# Filter SNPs that are significant outliers, here defined as having a Z-score > 2
significant_snps_z <- combined_data %>%
  filter(z_score > 2) %>%
  arrange(desc(z_score))


# Calculate median and MAD (Median Absolute Deviation)
median_effect <- median(combined_data$mean_effect, na.rm = TRUE)
mad_effect <- mad(combined_data$mean_effect, na.rm = TRUE, constant = 1)

# Calculate robust Z-scores (modified Z-scores using MAD)
combined_data$robust_z_score <- abs((combined_data$mean_effect - median_effect) / mad_effect)

# Filter SNPs that are significant outliers, here defined as having a robust Z-score > 3
significant_snps_mad <- combined_data %>%
  filter(robust_z_score > 3) %>%
  arrange(desc(robust_z_score))

# Combine and compare significant SNPs identified by both methods
combined_significant_snps <- significant_snps_z %>%
  rename(z_score_significance = z_score) %>%
  inner_join(significant_snps_mad %>%
               rename(mad_score_significance = robust_z_score), by = c("chrom", "position", "rs", "mean_effect"))

# Print combined significant SNPs
print(combined_significant_snps)

###################
# First, add a new highlight column to the combined_data
combined_data$significant_highlight <- ifelse(combined_data$rs %in% combined_significant_snps$rs, "Significant", "Other")

# Now, create the plot
ggplot(combined_data, aes(x = position, y = mean_effect)) +
  geom_point(aes(color = significant_highlight, size = significant_highlight), alpha = 0.2) +  # Base layer for size and color
  geom_point(data = subset(combined_data, significant_highlight == "Significant"), 
             aes(color = significant_highlight, size = significant_highlight), 
             shape = 21,   # Shape 21 is a filled circle
             fill = "red", # Set fill color for significant SNPs
             stroke = 1.5, # Set border thickness for clarity
             alpha = 1) +  # Opaque for significant SNPs
  scale_color_manual(values = c("Significant" = "red", "Other" = "gray")) +
  scale_size_manual(values = c("Significant" = 0.2, "Other" = 1)) +
  facet_wrap(~ chrom, scales = "free_x") +  # Faceting by chromosome
  theme_minimal() +
  labs(title = "SNP Effects by Chromosome with Significant SNPs Highlighted",
       x = "Position on Chromosome",
       y = "Mean Effect Size") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Improve readability of x-axis labels
        legend.position = "bottom")


##############################
#Pop 4
##############################
barley_pca <- read.csv("barley_pca_copy.csv")
individuals_in_groups <- split(barley_pca$sample.id, barley_pca$HCPC_group)
for(i in seq_along(individuals_in_groups)) {
  assign(paste0("Group", i), individuals_in_groups[[i]])
}

names(individuals_in_groups)
# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group4_column_names <- c("chrom", "position","rs", individuals_in_groups[["Group 4"]])  # Assuming these are valid column names

bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  dplyr::select(all_of(group4_column_names))

bio1.chrom2.plot <- bio1_Marker_Effects_chrom_2H %>%
  dplyr::select(all_of(group4_column_names))

bio1.chrom3.plot <- bio1_Marker_Effects_chrom_3H %>%
  dplyr::select(all_of(group4_column_names))

bio1.chrom4.plot <- bio1_Marker_Effects_chrom_4H %>%
  dplyr::select(all_of(group4_column_names))

bio1.chrom5.plot <- bio1_Marker_Effects_chrom_5H %>%
  dplyr::select(all_of(group4_column_names))

bio1.chrom6.plot <- bio1_Marker_Effects_chrom_6H %>%
  dplyr::select(all_of(group4_column_names))

bio1.chrom7.plot <- bio1_Marker_Effects_chrom_7H %>%
  dplyr::select(all_of(group4_column_names))


library(dplyr)

# Merge all chromosome dataframes
bio1_all_chroms_plot <- bind_rows(
  bio1.chrom1.plot,
  bio1.chrom2.plot,
  bio1.chrom3.plot,
  bio1.chrom4.plot,
  bio1.chrom5.plot,
  bio1.chrom6.plot,
  bio1.chrom7.plot
)

# Check the combined dataframe
print(bio1_all_chroms_plot)

###### ALL CHROM
library(dplyr)
library(purrr)

# Assuming 'all_chrom_data' is a list with data frames for each chromosome
# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(bio1_all_chroms_plot)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Calculate summary statistics for each SNP
summary_stats <- combined_data %>%
  group_by(rs) %>%
  summarise(min_effect = min(mean_effect, na.rm = TRUE),
            max_effect = max(mean_effect, na.rm = TRUE),
            median_effect = median(mean_effect, na.rm = TRUE))

###################
library(dplyr)

# Assuming 'combined_data' contains the mean effects for each SNP
mean_effect <- mean(combined_data$mean_effect, na.rm = TRUE)
std_dev_effect <- sd(combined_data$mean_effect, na.rm = TRUE)

# Calculate Z-scores for each SNP
combined_data$z_score <- abs((combined_data$mean_effect - mean_effect) / std_dev_effect)

# Filter SNPs that are significant outliers, here defined as having a Z-score > 2
significant_snps_z <- combined_data %>%
  filter(z_score > 2) %>%
  arrange(desc(z_score))


# Calculate median and MAD (Median Absolute Deviation)
median_effect <- median(combined_data$mean_effect, na.rm = TRUE)
mad_effect <- mad(combined_data$mean_effect, na.rm = TRUE, constant = 1)

# Calculate robust Z-scores (modified Z-scores using MAD)
combined_data$robust_z_score <- abs((combined_data$mean_effect - median_effect) / mad_effect)

# Filter SNPs that are significant outliers, here defined as having a robust Z-score > 3
significant_snps_mad <- combined_data %>%
  filter(robust_z_score > 3) %>%
  arrange(desc(robust_z_score))

# Combine and compare significant SNPs identified by both methods
combined_significant_snps <- significant_snps_z %>%
  rename(z_score_significance = z_score) %>%
  inner_join(significant_snps_mad %>%
               rename(mad_score_significance = robust_z_score), by = c("chrom", "position", "rs", "mean_effect"))

# Print combined significant SNPs
print(combined_significant_snps)

###################
# First, add a new highlight column to the combined_data
combined_data$significant_highlight <- ifelse(combined_data$rs %in% combined_significant_snps$rs, "Significant", "Other")

# Now, create the plot
ggplot(combined_data, aes(x = position, y = mean_effect)) +
  geom_point(aes(color = significant_highlight, size = significant_highlight), alpha = 0.2) +  # Base layer for size and color
  geom_point(data = subset(combined_data, significant_highlight == "Significant"), 
             aes(color = significant_highlight, size = significant_highlight), 
             shape = 21,   # Shape 21 is a filled circle
             fill = "red", # Set fill color for significant SNPs
             stroke = 1.5, # Set border thickness for clarity
             alpha = 1) +  # Opaque for significant SNPs
  scale_color_manual(values = c("Significant" = "red", "Other" = "gray")) +
  scale_size_manual(values = c("Significant" = 0.2, "Other" = 1)) +
  facet_wrap(~ chrom, scales = "free_x") +  # Faceting by chromosome
  theme_minimal() +
  labs(title = "SNP Effects by Chromosome with Significant SNPs Highlighted",
       x = "Position on Chromosome",
       y = "Mean Effect Size") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Improve readability of x-axis labels
        legend.position = "bottom")

##############################
#Pop 5
##############################
barley_pca <- read.csv("barley_pca_copy.csv")
individuals_in_groups <- split(barley_pca$sample.id, barley_pca$HCPC_group)
for(i in seq_along(individuals_in_groups)) {
  assign(paste0("Group", i), individuals_in_groups[[i]])
}

names(individuals_in_groups)
# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group5_column_names <- c("chrom", "position","rs", individuals_in_groups[["Group 5"]])  # Assuming these are valid column names

bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  dplyr::select(all_of(group5_column_names))

bio1.chrom2.plot <- bio1_Marker_Effects_chrom_2H %>%
  dplyr::select(all_of(group5_column_names))

bio1.chrom3.plot <- bio1_Marker_Effects_chrom_3H %>%
  dplyr::select(all_of(group5_column_names))

bio1.chrom4.plot <- bio1_Marker_Effects_chrom_4H %>%
  dplyr::select(all_of(group5_column_names))

bio1.chrom5.plot <- bio1_Marker_Effects_chrom_5H %>%
  dplyr::select(all_of(group5_column_names))

bio1.chrom6.plot <- bio1_Marker_Effects_chrom_6H %>%
  dplyr::select(all_of(group5_column_names))

bio1.chrom7.plot <- bio1_Marker_Effects_chrom_7H %>%
  dplyr::select(all_of(group5_column_names))


library(dplyr)

# Merge all chromosome dataframes
bio1_all_chroms_plot <- bind_rows(
  bio1.chrom1.plot,
  bio1.chrom2.plot,
  bio1.chrom3.plot,
  bio1.chrom4.plot,
  bio1.chrom5.plot,
  bio1.chrom6.plot,
  bio1.chrom7.plot
)

# Check the combined dataframe
print(bio1_all_chroms_plot)

###### ALL CHROM
library(dplyr)
library(purrr)

# Assuming 'all_chrom_data' is a list with data frames for each chromosome
# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(bio1_all_chroms_plot)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Calculate summary statistics for each SNP
summary_stats <- combined_data %>%
  group_by(rs) %>%
  summarise(min_effect = min(mean_effect, na.rm = TRUE),
            max_effect = max(mean_effect, na.rm = TRUE),
            median_effect = median(mean_effect, na.rm = TRUE))

###################
library(dplyr)

# Assuming 'combined_data' contains the mean effects for each SNP
mean_effect <- mean(combined_data$mean_effect, na.rm = TRUE)
std_dev_effect <- sd(combined_data$mean_effect, na.rm = TRUE)

# Calculate Z-scores for each SNP
combined_data$z_score <- abs((combined_data$mean_effect - mean_effect) / std_dev_effect)

# Filter SNPs that are significant outliers, here defined as having a Z-score > 2
significant_snps_z <- combined_data %>%
  filter(z_score > 2) %>%
  arrange(desc(z_score))


# Calculate median and MAD (Median Absolute Deviation)
median_effect <- median(combined_data$mean_effect, na.rm = TRUE)
mad_effect <- mad(combined_data$mean_effect, na.rm = TRUE, constant = 1)

# Calculate robust Z-scores (modified Z-scores using MAD)
combined_data$robust_z_score <- abs((combined_data$mean_effect - median_effect) / mad_effect)

# Filter SNPs that are significant outliers, here defined as having a robust Z-score > 3
significant_snps_mad <- combined_data %>%
  filter(robust_z_score > 3) %>%
  arrange(desc(robust_z_score))

# Combine and compare significant SNPs identified by both methods
combined_significant_snps <- significant_snps_z %>%
  rename(z_score_significance = z_score) %>%
  inner_join(significant_snps_mad %>%
               rename(mad_score_significance = robust_z_score), by = c("chrom", "position", "rs", "mean_effect"))

# Print combined significant SNPs
print(combined_significant_snps)

###################
# First, add a new highlight column to the combined_data
combined_data$significant_highlight <- ifelse(combined_data$rs %in% combined_significant_snps$rs, "Significant", "Other")

# Now, create the plot
ggplot(combined_data, aes(x = position, y = mean_effect)) +
  geom_point(aes(color = significant_highlight, size = significant_highlight), alpha = 0.2) +  # Base layer for size and color
  geom_point(data = subset(combined_data, significant_highlight == "Significant"), 
             aes(color = significant_highlight, size = significant_highlight), 
             shape = 21,   # Shape 21 is a filled circle
             fill = "red", # Set fill color for significant SNPs
             stroke = 1.5, # Set border thickness for clarity
             alpha = 1) +  # Opaque for significant SNPs
  scale_color_manual(values = c("Significant" = "red", "Other" = "gray")) +
  scale_size_manual(values = c("Significant" = 0.2, "Other" = 1)) +
  facet_wrap(~ chrom, scales = "free_x") +  # Faceting by chromosome
  theme_minimal() +
  labs(title = "SNP Effects by Chromosome with Significant SNPs Highlighted",
       x = "Position on Chromosome",
       y = "Mean Effect Size") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Improve readability of x-axis labels
        legend.position = "bottom")

# Write the data frame to a CSV file
#write.csv(combined_significant_snps, "combined_significant_snps_pop1.csv", row.names = FALSE)
#write.csv(combined_significant_snps, "combined_significant_snps_pop2.csv", row.names = FALSE)
#write.csv(combined_significant_snps, "combined_significant_snps_pop3.csv", row.names = FALSE)
#write.csv(combined_significant_snps, "combined_significant_snps_pop4.csv", row.names = FALSE)
write.csv(combined_significant_snps, "combined_significant_snps_pop5.csv", row.names = FALSE)


#########################################
#SUPPLEMENTAL FIGURE 6
#########################################
#Overlaps in highest z-score SNPs across populations
setwd("~/R/barley_collab/parental_selection/cross_validation/")
#SUP FIG 3A
df <- read.csv("positive_GEBV_snp_overlaps.csv")
#df <- read.csv("negative_GEBV_snp_overlaps.csv")

# Load necessary library
library(VennDiagram)

pop1 <- unique(df$Pop1[!is.na(df$Pop1)])
pop2 <- unique(df$Pop2[!is.na(df$Pop2)])
pop3 <- unique(df$Pop3[!is.na(df$Pop3)])
pop4 <- unique(df$Pop4[!is.na(df$Pop4)])
pop5 <- unique(df$Pop5[!is.na(df$Pop5)])

# Use the list method in VennDiagram to generate a Venn diagram and draw on device
venn.plot <- venn.diagram(
  x = list(
    Pop1 = pop1,
    Pop2 = pop2,
    Pop3 = pop3,
    Pop4 = pop4,
    Pop5 = pop5
  ),
  category.names = c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"),
  filename = NULL
)

# Render the plot
grid.draw(venn.plot)

#SUP FIG 3B
df <- read.csv("negative_GEBV_snp_overlaps.csv")

# Load necessary library
library(VennDiagram)

pop1 <- unique(df$Pop1[!is.na(df$Pop1)])
pop2 <- unique(df$Pop2[!is.na(df$Pop2)])
pop3 <- unique(df$Pop3[!is.na(df$Pop3)])
pop4 <- unique(df$Pop4[!is.na(df$Pop4)])
pop5 <- unique(df$Pop5[!is.na(df$Pop5)])

# Use the list method in VennDiagram to generate a Venn diagram and draw on device
venn.plot <- venn.diagram(
  x = list(
    Pop1 = pop1,
    Pop2 = pop2,
    Pop3 = pop3,
    Pop4 = pop4,
    Pop5 = pop5
  ),
  category.names = c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"),
  filename = NULL
)

# Render the plot
grid.draw(venn.plot)


#########################################
#SUPPLEMENTAL FIGURE 7
#########################################

#########################################
#comparing the 31 and 100 core predictions
#########################################
library(dplyr)
library(ggplot2)

# Set working directory to where the results are saved
setwd("~/R/barley_collab/parental_selection/core_comparisons/")

# Load results for Core 31
rrblup_kfold10_31 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/31/xval_rrblup_kfold_10.RData')
rrblup_kfold10_31 <- rrblup_kfold10_31$xval.result
rrblup_kfold10_31$r.mean <- as.numeric(rrblup_kfold10_31$r.mean)
rrblup_kfold10_31$r.sd <- as.numeric(rrblup_kfold10_31$r.sd)
rrblup_kfold10_31$core_size <- "Core 31"
rrblup_kfold10_31$model <- "RR-BLUP"

gauss_kfold_10_31 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/31/xval_GAUSS_kfold_10.RData')
gauss_kfold_10_31 <- gauss_kfold_10_31$xval.result
gauss_kfold_10_31$r.mean <- as.numeric(gauss_kfold_10_31$r.mean)
gauss_kfold_10_31$r.sd <- as.numeric(gauss_kfold_10_31$r.sd)
gauss_kfold_10_31$core_size <- "Core 31"
gauss_kfold_10_31$model <- "Gaussian Kernel"

EXP_kfold_10_31 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/31/xval_EXP_kfold_10.RData')
EXP_kfold_10_31 <- EXP_kfold_10_31$xval.result
EXP_kfold_10_31$r.mean <- as.numeric(EXP_kfold_10_31$r.mean)
EXP_kfold_10_31$r.sd <- as.numeric(EXP_kfold_10_31$r.sd)
EXP_kfold_10_31$core_size <- "Core 31"
EXP_kfold_10_31$model <- "Exponential Kernel"

bayescpi_kfold_10_31 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/31/xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10_31 <- bayescpi_kfold_10_31$xval.result
bayescpi_kfold_10_31$r.mean <- as.numeric(bayescpi_kfold_10_31$r.mean)
bayescpi_kfold_10_31$r.sd <- as.numeric(bayescpi_kfold_10_31$r.sd)
bayescpi_kfold_10_31$core_size <- "Core 31"
bayescpi_kfold_10_31$model <- "BayesCpi"

# Load results for Core 100
rrblup_kfold10_100 <- readRDS("~/R/barley_collab/parental_selection/core_comparisons/100/xval_rrblup_kfold_10.RData")
rrblup_kfold10_100 <- rrblup_kfold10_100$xval.result
rrblup_kfold10_100$r.mean <- as.numeric(rrblup_kfold10_100$r.mean)
rrblup_kfold10_100$r.sd <- as.numeric(rrblup_kfold10_100$r.sd)
rrblup_kfold10_100$core_size <- "Core 100"
rrblup_kfold10_100$model <- "RR-BLUP"

gauss_kfold_10_100 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/100/xval_GAUSS_kfold_10.RData')
gauss_kfold_10_100 <- gauss_kfold_10_100$xval.result
gauss_kfold_10_100$r.mean <- as.numeric(gauss_kfold_10_100$r.mean)
gauss_kfold_10_100$r.sd <- as.numeric(gauss_kfold_10_100$r.sd)
gauss_kfold_10_100$core_size <- "Core 100"
gauss_kfold_10_100$model <- "Gaussian Kernel"

EXP_kfold_10_100 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/100/xval_EXP_kfold_10.RData')
EXP_kfold_10_100 <- EXP_kfold_10_100$xval.result
EXP_kfold_10_100$r.mean <- as.numeric(EXP_kfold_10_100$r.mean)
EXP_kfold_10_100$r.sd <- as.numeric(EXP_kfold_10_100$r.sd)
EXP_kfold_10_100$core_size <- "Core 100"
EXP_kfold_10_100$model <- "Exponential Kernel"

bayescpi_kfold_10_100 <- readRDS('~/R/barley_collab/parental_selection/core_comparisons/100/xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10_100 <- bayescpi_kfold_10_100$xval.result
bayescpi_kfold_10_100$r.mean <- as.numeric(bayescpi_kfold_10_100$r.mean)
bayescpi_kfold_10_100$r.sd <- as.numeric(bayescpi_kfold_10_100$r.sd)
bayescpi_kfold_10_100$core_size <- "Core 100"
bayescpi_kfold_10_100$model <- "BayesCpi"

# Combine results into a single dataframe
combined_results <- bind_rows(
  rrblup_kfold10_31,
  gauss_kfold_10_31,
  EXP_kfold_10_31,
  bayescpi_kfold_10_31,
  rrblup_kfold10_100,
  gauss_kfold_10_100,
  EXP_kfold_10_100,
  bayescpi_kfold_10_100
)

# Organize dataframe for plotting
combined_results$model <- factor(combined_results$model, levels = c("RR-BLUP", "Gaussian Kernel", "Exponential Kernel", "BayesCpi"))
combined_results$core_size <- factor(combined_results$core_size, levels = c("Core 31", "Core 100"))

# Filter the relevant traits
all_bio <- combined_results %>% filter(trait %in% c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19'))

# Ensure all r.sd are numeric
all_bio$r.sd <- as.numeric(all_bio$r.sd)

# Plot the results
ggplot(all_bio, aes(y = r.mean, x = model, color = core_size)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(color = "Core Size", y = "Prediction Accuracy", x = "Model")



#########################################
#SUPPLEMENTAL FIGURE 8
#########################################
#comparing the 31 and 100 core
#########################################
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory and load the dataset
setwd("~/R/barley_collab/parental_selection/core_comparisons/")
data <- read.csv('comparing_cores_geav.csv')

# Display the first few rows of the data to check
head(data)

# Reshape the data from wide to long format
long_data <- data %>%
  pivot_longer(cols = c(bio1_geav_n31, bio1_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio1_range), names_to = "Range", values_to = "RangeValue")

# Display the first few rows of the reshaped dataset to check
head(long_data)


ggplot(long_data, aes(x = RangeValue, y = GEBV, color = Core)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of GEAV between Core Collections for bio1",
       x = "bio1_range",
       y = "GEAV") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ HCPC_group, scales = "fixed") +
  scale_color_manual(values = c("bio1_geav_n31" = "blue", "bio1_geav_n100" = "red"))


ggplot(long_data, aes(x = RangeValue, y = GEBV, color = Core)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of GEAV between Core Collections for bio1",
       x = "bio1_range",
       y = "GEAV") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ HCPC_group, scales = "fixed")


###### BIO 3
# Set working directory and load the dataset
setwd("~/R/barley_collab/parental_selection/core_comparisons/")
data <- read.csv('comparing_cores_geav.csv')

# Display the first few rows of the data to check
head(data)

# Reshape the data from wide to long format
long_data <- data %>%
  pivot_longer(cols = c(bio3_geav_n31, bio3_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio3_range), names_to = "Range", values_to = "RangeValue")

# Display the first few rows of the reshaped dataset to check
head(long_data)


ggplot(long_data, aes(x = RangeValue, y = GEBV, color = Core)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of GEAV between Core Collections for bio3",
       x = "bio3_range",
       y = "GEAV") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ HCPC_group, scales = "fixed")


###### BIO 4
# Set working directory and load the dataset
setwd("~/R/barley_collab/parental_selection/core_comparisons/")
data <- read.csv('comparing_cores_geav.csv')

# Display the first few rows of the data to check
head(data)

# Reshape the data from wide to long format
long_data <- data %>%
  pivot_longer(cols = c(bio4_geav_n31, bio4_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio4_range), names_to = "Range", values_to = "RangeValue")

# Display the first few rows of the reshaped dataset to check
head(long_data)


ggplot(long_data, aes(x = RangeValue, y = GEBV, color = Core)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of GEAV between Core Collections for bio4",
       x = "bio4_range",
       y = "GEAV") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ HCPC_group, scales = "fixed")

###### BIO 6
# Set working directory and load the dataset
setwd("~/R/barley_collab/parental_selection/core_comparisons/")
data <- read.csv('comparing_cores_geav.csv')

# Display the first few rows of the data to check
head(data)

# Reshape the data from wide to long format
long_data <- data %>%
  pivot_longer(cols = c(bio6_geav_n31, bio6_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio6_range), names_to = "Range", values_to = "RangeValue")

# Display the first few rows of the reshaped dataset to check
head(long_data)


ggplot(long_data, aes(x = RangeValue, y = GEBV, color = Core)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of GEAV between Core Collections for bio6",
       x = "bio6_range",
       y = "GEAV") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ HCPC_group, scales = "fixed")


###### BIO 11
# Set working directory and load the dataset
setwd("~/R/barley_collab/parental_selection/core_comparisons/")
data <- read.csv('comparing_cores_geav.csv')

# Display the first few rows of the data to check
head(data)

# Reshape the data from wide to long format
long_data <- data %>%
  pivot_longer(cols = c(bio11_geav_n31, bio11_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio11_range), names_to = "Range", values_to = "RangeValue")

# Display the first few rows of the reshaped dataset to check
head(long_data)


ggplot(long_data, aes(x = RangeValue, y = GEBV, color = Core)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of GEAV between Core Collections for bio11",
       x = "bio11_range",
       y = "GEAV") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ HCPC_group, scales = "fixed")


###### BIO 14
# Set working directory and load the dataset
setwd("~/R/barley_collab/parental_selection/core_comparisons/")
data <- read.csv('comparing_cores_geav.csv')

# Display the first few rows of the data to check
head(data)

# Reshape the data from wide to long format
long_data <- data %>%
  pivot_longer(cols = c(bio14_geav_n31, bio14_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio14_range), names_to = "Range", values_to = "RangeValue")

# Display the first few rows of the reshaped dataset to check
head(long_data)


ggplot(long_data, aes(x = RangeValue, y = GEBV, color = Core)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of GEAV between Core Collections for bio14",
       x = "bio14_range",
       y = "GEAV") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ HCPC_group, scales = "fixed")


###### BIO 17
# Set working directory and load the dataset
setwd("~/R/barley_collab/parental_selection/core_comparisons/")
data <- read.csv('comparing_cores_geav.csv')

# Display the first few rows of the data to check
head(data)

# Reshape the data from wide to long format
long_data <- data %>%
  pivot_longer(cols = c(bio17_geav_n31, bio17_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio17_range), names_to = "Range", values_to = "RangeValue")

# Display the first few rows of the reshaped dataset to check
head(long_data)


ggplot(long_data, aes(x = RangeValue, y = GEBV, color = Core)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of GEAV between Core Collections for bio17",
       x = "bio17_range",
       y = "GEAV") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ HCPC_group, scales = "fixed")


#########################################


#########################################
#SPECIES DISTRIBUTION 
#########################################

###########
#Pop 1
###########
library(ggplot2)
library(viridis)
species_distribution <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop1.tif")
#species_distribution <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop2.tif")
#species_distribution <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop3.tif")
#species_distribution <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop4.tif")
#species_distribution <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop5.tif")

species_distribution_reduced <- aggregate(species_distribution, fact=5, fun=mean) # 'fact' is the aggregation factor

# Convert to data frame
df_reduced <- as.data.frame(rasterToPoints(species_distribution_reduced))
names(df_reduced) <- c("Longitude", "Latitude", "Value")

p4<- ggplot(data = df_reduced, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Pop 1") +
  theme_minimal()
p4
#ggsave("~annamccormick/R/barley_collab/parental_selection/SDM/pop1_heatmap.pdf", plot = p4, width = 12, height = 5.39, device = "pdf")

###########
#Pop 2
###########
species_distribution <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop2.tif")
species_distribution_reduced <- aggregate(species_distribution, fact=5, fun=mean) # 'fact' is the aggregation factor

# Convert to data frame
df_reduced <- as.data.frame(rasterToPoints(species_distribution_reduced))
names(df_reduced) <- c("Longitude", "Latitude", "Value")

p4<- ggplot(data = df_reduced, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Pop 2") +
  theme_minimal()
p4
#ggsave("~annamccormick/R/barley_collab/parental_selection/SDM/pop2_heatmap.pdf", plot = p4, width = 12, height = 5.39, device = "pdf")

###########
#Pop 3
###########
species_distribution <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop3.tif")
species_distribution_reduced <- aggregate(species_distribution, fact=5, fun=mean) # 'fact' is the aggregation factor

# Convert to data frame
df_reduced <- as.data.frame(rasterToPoints(species_distribution_reduced))
names(df_reduced) <- c("Longitude", "Latitude", "Value")

p4<- ggplot(data = df_reduced, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Pop 3") +
  theme_minimal()
p4
#ggsave("~annamccormick/R/barley_collab/parental_selection/SDM/pop3_heatmap.pdf", plot = p4, width = 12, height = 5.39, device = "pdf")

###########
#Pop 4
###########
species_distribution <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop4.tif")
species_distribution_reduced <- aggregate(species_distribution, fact=5, fun=mean) # 'fact' is the aggregation factor

# Convert to data frame
df_reduced <- as.data.frame(rasterToPoints(species_distribution_reduced))
names(df_reduced) <- c("Longitude", "Latitude", "Value")

p4<- ggplot(data = df_reduced, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Pop 4") +
  theme_minimal()
p4
#ggsave("~annamccormick/R/barley_collab/parental_selection/SDM/pop4_heatmap.pdf", plot = p4, width = 12, height = 5.39, device = "pdf")

###########
#Pop 5
###########
species_distribution <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop5.tif")
species_distribution_reduced <- aggregate(species_distribution, fact=5, fun=mean) # 'fact' is the aggregation factor

# Convert to data frame
df_reduced <- as.data.frame(rasterToPoints(species_distribution_reduced))
names(df_reduced) <- c("Longitude", "Latitude", "Value")

p4<- ggplot(data = df_reduced, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Pop 5") +
  theme_minimal()
p4
#ggsave("~annamccormick/R/barley_collab/parental_selection/SDM/pop5_heatmap.pdf", plot = p4, width = 12, height = 5.39, device = "pdf")

#Plotting the 5x population SDMs

library(raster)
library(rworldmap)
# Load your rasters
r1 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop1.tif")
r2 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop2.tif")
r3 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop3.tif")
r4 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop4.tif")
r5 <- raster("~annamccormick/R/barley_collab/parental_selection/SDM/SDM_pop5.tif")

# Create a mask for values greater than 0.2
s1 <- calc(r1, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )
s2 <- calc(r2, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )
s3 <- calc(r3, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )
s4 <- calc(r4, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )
s5 <- calc(r5, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)} )

# Prepare a map to plot on
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster Datasets")

# Define colors and add transparency
#colors <- c("#00FF7F44", "#1E90FF44", "#B2222244")  # RGBA values

# Define unique colors for each raster
#colors <- c('red', 'yellow', 'green', 'blue', 'purple')

colors <- c(adjustcolor('red', alpha=0.5), 
            adjustcolor('yellow', alpha=0.5),
            adjustcolor('green', alpha=0.5),
            adjustcolor('blue', alpha=0.2),
            adjustcolor('purple', alpha=0.2))

# Plot the masked rasters on the world map
plot(s1, col=colors[1], legend =F, add=T)
plot(s2, col=colors[2], legend =F, add=T)
plot(s3, col=colors[3], legend =F, add=T)
plot(s4, col=colors[4], legend =F, add=T)
plot(s5, col=colors[5], legend =F, add=T)


# Add a custom legend
legend("topright", legend=c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"), fill=colors, bty="n")


######EACH 
#Pop1
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster Datasets")

colors <- c(adjustcolor('red', alpha=0.5), 
            adjustcolor('yellow', alpha=0.7),
            adjustcolor('green', alpha=0.5),
            adjustcolor('blue', alpha=0.2),
            adjustcolor('purple', alpha=0.2))

#Pop2
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster Datasets")

plot(s2, col=colors[2], legend =F, add=T)
# Adding a legend specifically for Pop2 with yellow color
legend("topright", legend=c("Pop2"), fill=colors[2], bty="n")


#Pop3
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster Datasets")

plot(s3, col=colors[3], legend =F, add=T)
legend("topright", legend=c("Pop3"), fill=colors[3], bty="n")

#Pop4
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster Datasets")

plot(s4, col=colors[4], legend =F, add=T)
legend("topright", legend=c("Pop4"), fill=colors[4], bty="n")

#Pop5
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster Datasets")

plot(s5, col=colors[5], legend =F, add=T)
legend("topright", legend=c("Pop5"), fill=colors[5], bty="n")


#### FIN SDM
