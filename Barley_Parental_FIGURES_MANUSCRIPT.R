#######################################################################################################
#MAIN FIGURES
#######################################################################################################

#######################################################################################################
#FIGURE 1
######################################################################################################
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
setwd("~/R/barley_collab/parental_selection/")
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
f2<-f+ scale_color_manual(values = c( "purple",  "red",  "blue", "yellow1", "green", "yellow2", "grey"))

ggsave("~/R/barley_collab/parental_selection/hierarchical_cluster_dendrogram.pdf", 
       plot = f2, width = 10, height = 7)
#########################################
#FIGURE 1B
#########################################
#Plot 2D PCA 
library("ggrepel")
setwd("~/R/barley_collab/parental_selection/cross_validation/")
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
p + scale_color_manual(values = c( "red",  "yellow2",  "green", "blue", "purple"))

p <- ggplot(data=pca, aes(x=EV1, y=EV2, color=HCPC_group)) +
  geom_point() +
  scale_shape_manual(values=seq(0, 15)) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (9.84%)") + 
  ylab("PC 2 (6.96%)") +
  scale_color_manual(values = c("red", "yellow2", "green", "blue", "purple")) +
  theme(
    axis.title = element_text(size = 18),  # Increase axis titles size
    axis.text = element_text(size = 16),   # Increase axis labels size
    legend.text = element_text(size = 16), # Increase legend text size
    legend.title = element_text(size = 16) # Increase legend title size
  )
print(p)


#########################################
#FIGURE 1C
#########################################
#MAP
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

# Set working directory
setwd("~/R/barley_collab/parental_selection/cross_validation/")
pca <- read.csv("barley_pca_copy.csv")

# Load world map dataset using rnaturalearth
world <- ne_countries(scale = "medium", returnclass = "sf")

# Set up a base plot
base_plot <- ggplot() +
  geom_sf(data = world, fill = "white") +
  theme_minimal()

# Set the limits for longitude and latitude to focus on Europe, Asia, and Africa
europe_asia_africa_plot <- base_plot +
  coord_sf(xlim = c(-30, 160), ylim = c(-40, 80))

# Add the points with increased size and fill them
final_plot <- europe_asia_africa_plot +
  geom_point(data = pca, 
             aes(x = Longitude, y = Latitude, fill = HCPC_group, color = HCPC_group), 
             pch = 21, size = 2, alpha = 0.8) + # Increased size, and alpha for transparency
  scale_fill_manual(values = c("Population 1" = "red", "Population 2" = "yellow", 
                               "Population 3" = "green", "Population 4" = "blue", 
                               "Population 5" = "purple")) +
  scale_color_manual(values = c("Population 1" = "black", "Population 2" = "black", 
                                "Population 3" = "black", "Population 4" = "black", 
                                "Population 5" = "black")) + # Add black outline to the points
  theme(
    legend.title = element_text(size = 12), # Adjust legend title size
    legend.text = element_text(size = 10)   # Adjust legend text size
  )

# Display the final plot
print(final_plot)

# Save the plot as a PDF
ggsave("europe_asia_africa_pca_plot.pdf", plot = final_plot, width = 10, height = 7)

######################################################################################################
#FIGURE 2
######################################################################################################

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

######################################################################################################
#FIGURE 3
######################################################################################################
#core 31
######################
setwd("~/R/barley_collab/parental_selection/cross_validation/")
data2 <-read.csv('rrblup_GEBV_data_barley_plots_copy.csv')

#bio1
p1<-ggplot(data2, aes(x = HCPC_group, y = bio1, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio1 - Annual Mean Temperature",
       x = "HCPC Group",
       y = "GEAV bio1") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
p1


#bio 3
p2<-ggplot(data2, aes(x = HCPC_group, y = bio3, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio3 - Isothermality (BIO2/BIO7) ",
       x = "HCPC Group",
       y = "GEAV bio3") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio4
p3<-ggplot(data2, aes(x = HCPC_group, y = bio4, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio4 - temperature seasonality ",
       x = "HCPC Group",
       y = "GEAV bio4") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 6
p4<-ggplot(data2, aes(x = HCPC_group, y = bio6, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio6 - Min Temperature of Coldest Month",
       x = "HCPC Group",
       y = "GEAV bio6") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 11
p5<-ggplot(data2, aes(x = HCPC_group, y = bio11, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio11 - Mean Temperature of Coldest Quarter",
       x = "HCPC Group",
       y = "GEAV bio11") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 14
p6<-ggplot(data2, aes(x = HCPC_group, y = bio14, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio14 - Precipitation of Driest Month",
       x = "HCPC Group",
       y = "GEAV bio14") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio17
p7<-ggplot(data2, aes(x = HCPC_group, y = bio17, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio17 - Precipitation of Driest Quarter",
       x = "HCPC Group",
       y = "GEAV bio17") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")


library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 3)


#########################################
#Core collection identification from datasets
#########################################
library(vcfR)
library(poppr)
library(adegenet)

setwd("~/R/barley_collab/parental_selection/")
vcf <- read.vcfR("Supplemental_dataset_1.vcf")

####convert genlight object
x <- vcfR2genlight(vcf)

#create distance matrix
x.dist <- dist(x)
ploidy(x) <- 2

library(corehunter)
# precomputed distance matrix
my.data <- distances(as.matrix(x.dist))
core <- sampleCore(my.data, size = 100)
core

write.csv(core, 'barley_core_lei.csv')

######################
#CORE 100
######################
#boxplots
setwd("~/R/barley_collab/parental_selection/new_core_100/")
#data2 <-read.csv('rrblup_GEBV_data_barley_plots_copy.csv')
data2 <-read.csv('rrblup_GEBV_data_barley_100_abbreviatedbio.csv')

#bio1
p1a<-ggplot(data2, aes(x = HCPC_group, y = bio1, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio1 - Annual Mean Temperature",
       x = "HCPC Group",
       y = "GEAV bio1") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
p1a

#bio 3
p2a<-ggplot(data2, aes(x = HCPC_group, y = bio3, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio3 - Isothermality (BIO2/BIO7) ",
       x = "HCPC Group",
       y = "GEAV bio3") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio4
p3a<-ggplot(data2, aes(x = HCPC_group, y = bio4, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio4 - temperature seasonality ",
       x = "HCPC Group",
       y = "GEAV bio4") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 6
p4a<-ggplot(data2, aes(x = HCPC_group, y = bio6, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio6 - Min Temperature of Coldest Month",
       x = "HCPC Group",
       y = "GEAV bio6") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 11
p5a<-ggplot(data2, aes(x = HCPC_group, y = bio11, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio11 - Mean Temperature of Coldest Quarter",
       x = "HCPC Group",
       y = "GEAV bio11") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio 14
p6a<-ggplot(data2, aes(x = HCPC_group, y = bio14, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio14 - Precipitation of Driest Month",
       x = "HCPC Group",
       y = "GEAV bio14") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

#bio17
p7a<-ggplot(data2, aes(x = HCPC_group, y = bio17, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Group 1" = "red", "Group 2" = "yellow", "Group 3" = "green", "Group 4" = "blue", "Group 5" = "purple")) +
  labs(title = "GEAV bio17 - Precipitation of Driest Quarter",
       x = "HCPC Group",
       y = "GEAV bio17") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")

library(gridExtra)
grid.arrange(p1a, p2a, p3a, p4a, p5a, p6a, p7a, ncol = 3)


######################################################################################################
#FIGURE 4
######################################################################################################
##############
#FIGURE 4A 
##############
library(ggplot2)
library(factoextra)
###########
#All populations
###########
setwd("~/R/barley_collab/parental_selection/cross_validation/")
data <- read.csv("n784_climate_data_n7.csv")
env_variables <- data[, 9:15]

# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)

# Combine PCA scores with grouping variable
pca_data <- data.frame(pca_result$x, HCPC_group = data$HCPC_group)

# Create a biplot with color by 'HCPC_group'
p <- fviz_pca_biplot(pca_result, 
                     geom = "point", 
                     col.ind = "black", 
                     col.var = "blue", 
                     addEllipses = FALSE, 
                     habillage = pca_data$HCPC_group) +
  ggtitle("PCA Biplot All Populations") +  # Title
  theme_classic() +  
  
  # Adjust text sizes using theme()
  theme(
    plot.title = element_text(size = 18),     # Increase title size
    axis.title = element_text(size = 14),     # Increase axis title size
    axis.text = element_text(size = 12),      # Increase axis tick label size
    legend.text = element_text(size = 12),    # Increase legend text size
    legend.title = element_text(size = 14),   # Increase legend title size
    axis.line = element_blank()               # Remove x and y axis lines
  )

# Adjust x and y axis limits
p1 <- p + xlim(-8, 8) + ylim(-8, 8)

# Display the plot
print(p1)
ggsave("~/R/barley_collab/parental_selection/figure_pannels/pca_biplot_Figure2F.pdf", plot = p1, width = 8, height = 6)

##############
#FIGURE 4B
##############
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Set working directory
setwd("~/R/barley_collab/parental_selection/cross_validation/")

# Read in your data
data2 <- read.csv('geav_values_lat_long.csv')

# Convert the data frame to a spatial object using the Latitude and Longitude columns
data_sf <- st_as_sf(data2, coords = c("Longitude", "Latitude"), crs = 4326)

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the shape for each HCPC group
#shape_map <- c("Group 1" = 21, "Group 2" = 22, "Group 3" = 23, "Group 4" = 25, "Group 5" = 24)
shape_map <- c("Population 1" = 21, "Population 2" = 22, "Population 3" = 23, "Population 4" = 25, "Population 5" = 24)

# Plotting the map 
p1 <- ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "black") +  # Adding world boundaries with black border
  geom_sf(data = data_sf, aes(color = bio1, shape = HCPC_group), size = 1, show.legend = TRUE) +  # Smaller shapes with size = 2
  scale_color_gradient(low = "blue", high = "red", name = "Bio1 Value") +
  scale_shape_manual(values = shape_map) +
  coord_sf(xlim = c(-20, 150), ylim = c(-40, 80)) +  # Set the limits similar to the second map
  labs(title = "Bio1 (Annual mean temperature) GEAV Values") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 6),  # Adjust the title size to 10 (same as second plot)
        legend.text = element_text(size = 6),  # Adjust the legend text size to 10
        legend.title = element_text(size = 6), # Adjust the legend title size to 10
        legend.box.spacing = unit(0.2, "cm")    # Adjust the legend box spacing
  ) +
  guides(color = guide_colorbar(barwidth = 0.5, barheight = 5))  # Adjust color bar width and height

p1


#######
# Plotting the map with country boundaries and adjusted range
p2<-ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "black") +  # Adding world boundaries with black border
  geom_sf(data = data_sf, aes(color = bio3, shape = HCPC_group), size = 1, show.legend = TRUE) + # Smaller shapes with size = 2
  scale_color_gradient(low = "blue", high = "red", name = "Bio3 Value") +
  scale_shape_manual(values = shape_map) +
  coord_sf(xlim = c(-20, 150), ylim = c(-40, 80)) +  # Set the limits similar to the second map
  labs(title = "Bio3 (Isothermality) GEAV Values") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 6),  # Adjust the title size
        legend.text = element_text(size = 6),  # Adjust the legend text size
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(0.2, "cm")  # Adjust the legend title size
  )+
  guides(color = guide_colorbar(barwidth = 0.5, barheight = 5)) 
p2


# Plotting the map with country boundaries and adjusted range
p3<-ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "black") +  # Adding world boundaries with black border
  geom_sf(data = data_sf, aes(color = bio4, shape = HCPC_group), size = 1, show.legend = TRUE) + # Smaller shapes with size = 2
  scale_color_gradient(low = "blue", high = "red", name = "Bio4 Value") +
  scale_shape_manual(values = shape_map) +
  coord_sf(xlim = c(-20, 150), ylim = c(-40, 80)) +  # Set the limits similar to the second map
  labs(title = "Bio4 (Temperature seasonality) GEAV Values") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 6),  # Adjust the title size
        legend.text = element_text(size = 6),  # Adjust the legend text size
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(0.2, "cm")  # Adjust the legend title size
  )+
  guides(color = guide_colorbar(barwidth = 0.5, barheight = 5)) 
p3


# Plotting the map with country boundaries and adjusted range
p4<-ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "black") +  # Adding world boundaries with black border
  geom_sf(data = data_sf, aes(color = bio6, shape = HCPC_group), size = 1, show.legend = TRUE) + # Smaller shapes with size = 2
  scale_color_gradient(low = "blue", high = "red", name = "Bio6 Value") +
  scale_shape_manual(values = shape_map) +
  coord_sf(xlim = c(-20, 150), ylim = c(-40, 80)) +  # Set the limits similar to the second map
  labs(title = "Bio6 (Min temperature of the coldest month) GEAV Values") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 6),  # Adjust the title size
        legend.text = element_text(size = 6),  # Adjust the legend text size
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(0.2, "cm")  # Adjust the legend title size
  )+
  guides(color = guide_colorbar(barwidth = 0.5, barheight = 5)) 
p4


# Plotting the map with country boundaries and adjusted range
p5<-ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "black") +  # Adding world boundaries with black border
  geom_sf(data = data_sf, aes(color = bio11, shape = HCPC_group), size = 1, show.legend = TRUE) + # Smaller shapes with size = 2
  scale_color_gradient(low = "blue", high = "red", name = "Bio11 Value") +
  scale_shape_manual(values = shape_map) +
  coord_sf(xlim = c(-20, 150), ylim = c(-40, 80)) +  # Set the limits similar to the second map
  labs(title = "Bio11 (Mean temperature of the coldest quarter) GEAV Values") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 6),  # Adjust the title size
        legend.text = element_text(size = 6),  # Adjust the legend text size
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(0.2, "cm")  # Adjust the legend title size
  )+
  guides(color = guide_colorbar(barwidth = 0.5, barheight = 5)) 
p5


# Plotting the map with country boundaries and adjusted range
p6<-ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "black") +  # Adding world boundaries with black border
  geom_sf(data = data_sf, aes(color = bio14, shape = HCPC_group), size = 1, show.legend = TRUE) + # Smaller shapes with size = 2
  scale_color_gradient(low = "blue", high = "red", name = "Bio14 Value") +
  scale_shape_manual(values = shape_map) +
  coord_sf(xlim = c(-20, 150), ylim = c(-40, 80)) +  # Set the limits similar to the second map
  labs(title = "Bio14 (Precipitation of the driest month) GEAV Values") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 6),  # Adjust the title size
        legend.text = element_text(size = 6),  # Adjust the legend text size
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(0.2, "cm")  # Adjust the legend title size
  )+
  guides(color = guide_colorbar(barwidth = 0.5, barheight = 5)) 
p6


# Plotting the map with country boundaries and adjusted range
p7<-ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "black") +  # Adding world boundaries with black border
  geom_sf(data = data_sf, aes(color = bio17, shape = HCPC_group), size = 1, show.legend = TRUE) + # Smaller shapes with size = 2
  scale_color_gradient(low = "blue", high = "red", name = "Bio17 Value") +
  scale_shape_manual(values = shape_map) +
  coord_sf(xlim = c(-20, 150), ylim = c(-40, 80)) +  # Set the limits similar to the second map
  labs(title = "Bio17 (Precipitation of driest quarter) GEAV Values") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 6),  # Adjust the title size
        legend.text = element_text(size = 6),  # Adjust the legend text size
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(0.2, "cm")  # Adjust the legend title size
  )+
  guides(color = guide_colorbar(barwidth = 0.5, barheight = 5)) 
p7

# Load the gridExtra package
library(gridExtra)

# Arrange the plots in a grid layout
grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 3)



######################################################################################################
#FIGURE 5
######################################################################################################
##############
#FIGURE 5A 
##############
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


#############################################
#keep relevant points in SDM range and PRINT table 
#############################################
##############
#Figure 5B
##############
#Population 1
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r1 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop1.tif")

# Create a mask for values greater than 0.2
s1 <- calc(r1, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
data2 <- read.csv('~/R/barley_collab/parental_selection/cross_validation/geav_values_lat_long.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s1, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s1, col=adjustcolor('red', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 1"), fill=adjustcolor('red', alpha=0.5), bty="n")

##############
#Figure 5C
##############
#Population 2
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r2 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop2.tif")

# Create a mask for values greater than 0.2
s2 <- calc(r2, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
data2 <- read.csv('~/R/barley_collab/parental_selection/cross_validation/geav_values_lat_long.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s2, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points_pop2.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s2, col=adjustcolor('yellow', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 2"), fill=adjustcolor('yellow', alpha=0.5), bty="n")

##############
#Figure 5D
##############
#Population 3
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r3 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop3.tif") #change

# Create a mask for values greater than 0.2
s3 <- calc(r3, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
data2 <- read.csv('~/R/barley_collab/parental_selection/cross_validation/geav_values_lat_long.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s3, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points_pop3.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s3, col=adjustcolor('green', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 3"), fill=adjustcolor('green', alpha=0.5), bty="n")

##############
#Figure 5E
##############
#Population 4
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r4 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop4.tif") #change

# Create a mask for values greater than 0.2
s4 <- calc(r4, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
data2 <- read.csv('~/R/barley_collab/parental_selection/cross_validation/geav_values_lat_long.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s4, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points_pop4.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s4, col=adjustcolor('blue', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 4"), fill=adjustcolor('blue', alpha=0.5), bty="n")

##############
#Figure 5F
##############
#Population 5
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r5 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop5.tif") #change

# Create a mask for values greater than 0.2
s5 <- calc(r5, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
data2 <- read.csv('~/R/barley_collab/parental_selection/cross_validation/geav_values_lat_long.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s5, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points_pop5.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s5, col=adjustcolor('purple', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 5"), fill=adjustcolor('purple', alpha=0.5), bty="n")

#############################################
#keep relevant points in SDM range and PRINT table 
#Core 100 results
#############################################
#Population 1
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r1 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop1.tif")

# Create a mask for values greater than 0.2
s1 <- calc(r1, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
#data2 <- read.csv('~/R/barley_collab/parental_selection/cross_validation/geav_values_lat_long.csv')
data2 <- read.csv('~/R/barley_collab/parental_selection/new_core_100/rrblup_GEBV_data_barley_100_abbreviatedbio_tables.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s1, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points_core100_pop1.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s1, col=adjustcolor('red', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 1"), fill=adjustcolor('red', alpha=0.5), bty="n")

#Population 2
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r2 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop2.tif")

# Create a mask for values greater than 0.2
s2 <- calc(r2, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
#data2 <- read.csv('~/R/barley_collab/parental_selection/cross_validation/geav_values_lat_long.csv')
data2 <- read.csv('~/R/barley_collab/parental_selection/new_core_100/rrblup_GEBV_data_barley_100_abbreviatedbio_tables.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s2, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points_core100_pop2.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s2, col=adjustcolor('yellow', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 2"), fill=adjustcolor('yellow', alpha=0.5), bty="n")

#Population 3
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r3 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop3.tif") #change

# Create a mask for values greater than 0.2
s3 <- calc(r3, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
data2 <- read.csv('~/R/barley_collab/parental_selection/new_core_100/rrblup_GEBV_data_barley_100_abbreviatedbio_tables.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s3, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points_core100_pop3.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s3, col=adjustcolor('green', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 3"), fill=adjustcolor('green', alpha=0.5), bty="n")

#Population 4
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r4 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop4.tif") #change

# Create a mask for values greater than 0.2
s4 <- calc(r4, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
data2 <- read.csv('~/R/barley_collab/parental_selection/new_core_100/rrblup_GEBV_data_barley_100_abbreviatedbio_tables.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s4, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points_core100_pop4.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s4, col=adjustcolor('blue', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 4"), fill=adjustcolor('blue', alpha=0.5), bty="n")

#Population 5
library(raster)
library(rworldmap)
library(plotrix)

# Load your raster file (SDM suitability map)
r5 <- raster("~/R/barley_collab/parental_selection/SDM/SDM_pop5.tif") #change

# Create a mask for values greater than 0.2
s5 <- calc(r5, fun=function(x){ x[x > 0.2] <- 1; x[x <= 0.2] <- NA; return(x)})

# Load the GEAV data with lat/lon coordinates
data2 <- read.csv('~/R/barley_collab/parental_selection/new_core_100/rrblup_GEBV_data_barley_100_abbreviatedbio_tables.csv')

# Extract raster values at the GEAV data points (long/lat)
geav_raster_values <- extract(s5, data2[,c("Longitude", "Latitude")])

# Keep only the points where the raster value is not NA (i.e., inside the red zones)
valid_points <- !is.na(geav_raster_values)
data2_trimmed <- data2[valid_points, ]

# Save the filtered data (valid points) to a CSV file
write.csv(data2_trimmed, file = "~/R/barley_collab/parental_selection/cross_validation/valid_points_core100_pop5.csv", row.names = FALSE)

# Normalize bio1 (GEAV) values for color mapping
bio1_normalized <- (data2_trimmed$bio1 - min(data2_trimmed$bio1)) / (max(data2_trimmed$bio1) - min(data2_trimmed$bio1))

# Define a color gradient based on bio1 (GEAV) values
colors_gradient <- colorRampPalette(c("red", "white", "blue"))(100)

# Map the bio1 values to the color gradient
point_colors <- colors_gradient[as.numeric(cut(bio1_normalized, breaks = 100))]

# Prepare the world map to plot as the background
world_map <- getMap(resolution = "low")
plot(world_map, col="#f2f2f2", border="#a0a0a0", lwd=0.5, main="Overlay of Raster and GEAV values")

# Plot the masked raster data
plot(s5, col=adjustcolor('purple', alpha=0.5), legend = F, add = T)

# Plot the filtered points (only those inside the red zones)
points(data2_trimmed$Longitude, data2_trimmed$Latitude, pch=21, bg=point_colors, col='black', cex=0.5)

# Add a legend for GEAV gradient
legend("bottomleft", legend = "GEAV Values", pch = 21, pt.bg = point_colors, col='black', bty = "n")

# Add a color legend for the bio1 values
color.legend(-180, -60, -120, -40, legend = round(seq(min(data2_trimmed$bio1), max(data2_trimmed$bio1), length.out = 5), 1), 
             rect.col = colors_gradient, gradient = "x", align = "rb")

# Optionally add a custom legend for the raster
legend("topright", legend=c("Population 5"), fill=adjustcolor('purple', alpha=0.5), bty="n")


######################################################################################################
#FIGURE 6
######################################################################################################
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


#######################################################################################################
#SUPPLEMENTAL FIGURES
#######################################################################################################

#######################################################################################################
#SUPPLEMENTAL FIGURE 1
#######################################################################################################
library(MetBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggnewscale)

##############
#K2
##############
# Load the Hokusai palette
hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 6)]  # Colors for K1 and K2

# Load your data
test <- read.csv("/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/barley_k2.csv")

# Sort by K1 and K2 percentages
test <- test %>%
  arrange(desc(K1), desc(K2))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type
type_colors <- c("Group_1_East_African" = "red", 
                 "Group_2_Levant/Mediterranean" = "yellow", 
                 "Group_3_North_African/Mediterranean" = "green", 
                 "Group_4_Northern_Europe" = "blue", 
                 "Group_5_Asian" = "#9467BD")

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.05) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )

print(final_plot)

##############
#K3
##############
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(MetBrewer)

hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 4, 6)]

# Load your data
test <- read.csv("/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/barley_k3.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 ~ "K1",
    K2 >= K1 & K2 >= K3 ~ "K2",
    K3 >= K1 & K3 >= K2 ~ "K3"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c("Group_1_East_African" = "red", 
                 "Group_2_Levant/Mediterranean" = "yellow", 
                 "Group_3_North_African/Mediterranean" = "green", 
                 "Group_4_Northern_Europe" = "blue", 
                 "Group_5_Asian" = "#9467BD")

# Hiroshige palette for K categories
k_colors <- hiroshige_palette[c(1, 4, 6)]

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1, K2, and K3 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1, K2, and K3
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.05) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )

print(final_plot)

##############
#K4
##############
#load fastSTRUCTURE data
data <- read.table("/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/testoutput_simple4a.4.meanQ", header = FALSE)

#Assign column names
colnames(data) <- c("K1", "K2", "K3", "K4")  # Replace with appropriate names

#Export CSV 
write.csv(data, "/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/barley_k4.csv", row.names = FALSE)

hiroshige_palette <- met.brewer("Hiroshige", 8)
#selected_colors <- hiroshige_palette[c(1, 4,6, 8)]
k_colors <- hiroshige_palette[c(1,4,6,8)]  # Colors for K

test <- read.csv("/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/barley_k4.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 ~ "K4"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c("Group_1_East_African" = "red", 
                 "Group_2_Levant/Mediterranean" = "yellow", 
                 "Group_3_North_African/Mediterranean" = "green", 
                 "Group_4_Northern_Europe" = "blue", 
                 "Group_5_Asian" = "#9467BD")

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 4, 6, 8)]  # Choose appropriate colors for K1 to K4

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1, K2, K3, and K4 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 to K4
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.05) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )

print(final_plot)

##############
#K5
##############
#load fastSTRUCTURE data
data <- read.table("/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/testoutput_simple5a.5.meanQ", header = FALSE)

#Assign column names
colnames(data) <- c("K1", "K2", "K3", "K4", "K5")  # Replace with appropriate names

#Export CSV 
write.csv(data, "/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/barley_k5.csv", row.names = FALSE)

hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test <- read.csv("/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/barley_k5.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 & K1 >= K5 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 & K2 >= K5 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 & K3 >= K5 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 & K4 >= K5 ~ "K4",
    K5 >= K1 & K5 >= K2 & K5 >= K3 & K5 >= K4 ~ "K5"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4), desc(K5))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c("Group_1_East_African" = "red", 
                 "Group_2_Levant/Mediterranean" = "yellow", 
                 "Group_3_North_African/Mediterranean" = "green", 
                 "Group_4_Northern_Europe" = "blue", 
                 "Group_5_Asian" = "#9467BD")

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 3, 5, 6, 8)]  # Choose appropriate colors for K1 to K5

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1, K2, K3, K4, and K5 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 to K5
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.05) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )

print(final_plot)

##############
#K6
##############
#load fastSTRUCTURE data
data <- read.table("/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/testoutput_simple6a.6.meanQ", header = FALSE)

#Assign column names
colnames(data) <- c("K1", "K2", "K3", "K4", "K5", "K6")  # Replace with appropriate names

#Export CSV 
write.csv(data, "/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/barley_k6.csv", row.names = FALSE)


# Load your data
test <- read.csv("/Users/annamccormick/R/barley_collab/parental_selection/fastSTRUCTURE/barley_k6.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 & K1 >= K5 & K1 >= K6 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 & K2 >= K5 & K2 >= K6 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 & K3 >= K5 & K3 >= K6 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 & K4 >= K5 & K4 >= K6 ~ "K4",
    K5 >= K1 & K5 >= K2 & K5 >= K3 & K5 >= K4 & K5 >= K6 ~ "K5",
    K6 >= K1 & K6 >= K2 & K6 >= K3 & K6 >= K4 & K6 >= K5 ~ "K6"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4), desc(K5), desc(K6))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5", "K6"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c("Group_1_East_African" = "red", 
                 "Group_2_Levant/Mediterranean" = "yellow", 
                 "Group_3_North_African/Mediterranean" = "green", 
                 "Group_4_Northern_Europe" = "blue", 
                 "Group_5_Asian" = "#9467BD")

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 3, 5, 6, 7, 8)]  # Choose appropriate colors for K1 to K6

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1, K2, K3, K4, K5, and K6 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 to K6
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.05) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )

print(final_plot)



#######################################################################################################
#SUPPLEMENTAL FIGURE 2
#######################################################################################################
#ON HPC
########################################
###Run all k-fold cross-validations       
##########################################
library(rrBLUP) #rrBLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) #hibayes package for Bayesian models
library(dplyr)

setwd("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_2/")

#load all functions from separate R file
source("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_2/xval_kfold_functions.R") 

#Read genotype data in rrBLUP format
gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = T)

#Read environmental data from .csv
envdat <- read.csv('Supplemental_dataset_5worldclim_all.csv', head = T) # full environmental dataset

#Read training dataset from .csv
trainingset <- read.csv("Supplemental_dataset_5_worldclim_training100.csv", head = T)

#Filter genotype data to include only rows where taxa are in the environmental data
gd2 <- gd1[gd1$taxa %in% envdat$Taxa,]
row.names(gd2) <- gd2$taxa # set taxa as rownames
g.in <- gd2[,-1] # Remove the first column (taxa) to create a genotype matrix
g.in <- as.matrix(g.in)

### Load Environmental Data ###
row.names(envdat) <- envdat$Taxa # set gen_id as rownames
row.names(trainingset) <- trainingset$Taxa

#Select unique identifier and environmental data columns for training set
y.trainset <- trainingset[,c(1,5:(ncol(trainingset)))] # select unique identifier and environmental data only 

#Select unique identifier and environmental data columns for the full dataset
y.in <- envdat[,c(1,4:ncol(envdat))] 

############
### RR-BLUP
############
#Prepare data for rrBLUP model
y.trainset.rr <- trainingset[,c(6:(ncol(trainingset)))]
y.in.rr <- envdat[,5:ncol(envdat)]

#Convert data to matrix format
y.trainset.mat <- as.matrix(y.trainset.rr)
y.in.mat <- as.matrix(y.in.rr)

#Run k-fold cross-validation for rrBLUP
xval_k10_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 10, reps = 50)
saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10.RData")

############
# Gaussian Kernel
############
K <- A.mat(g.in) #Compute relationship matrix using Gaussian Kernel
k_dist <- dist(K) #Calculate distance matrix from relationship matrix

#Run k-fold cross-validation for Gaussian Kernel
xval_k10_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_GAUSS, "xval_GAUSS_kfold_10.RData")

############
### Exponential Kernel
############
#Compute relationship matrix using Exponential Kernel
K.Exp=Kernel_computation(X=g.in, name="exponential", degree=NULL, nL=NULL)

#Set row and column names for the Exponential Kernel matrix
row.names(K.Exp) <- rownames(g.in)
colnames(K.Exp) <- rownames(g.in)

#Calculate distance matrix from Exponential Kernel relationship matrix
exp_dist <- dist(K.Exp) # Calculate Relationship Matrix\

#Run k-fold cross-validation for Exponential Kernel
xval_k10_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_EXP, "xval_EXP_kfold_10.RData")

#write.table(exp_dist, "distance_matrix_exponential_kernel.txt")
#write.table(k_dist, "distance_matrix_gaussian_kernel.txt")

############
#BayesCPi
############
#Run k-fold cross-validation for BayesCPi model
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
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)

#############################################################################
# Gaussian Kernel - K-fold
#Load cross-validation results for Gaussian Kernel
gauss_kfold_10 <- readRDS('xval_GAUSS_kfold_10.RData')
gauss_kfold_10 <- gauss_kfold_10$xval.result
gauss_kfold_10$r.mean <- as.numeric(gauss_kfold_10$r.mean)

################################################################################
# Exponential Kernel - K-fold
#Load cross-validation results for Exponential Kernel
EXP_kfold_10 <- readRDS('xval_EXP_kfold_10.RData')
EXP_kfold_10 <- EXP_kfold_10$xval.result
EXP_kfold_10$r.mean <- as.numeric(EXP_kfold_10$r.mean)

# BayesCpi - fractional
#bayescpi_frac90 <- readRDS("xval_BayesCpi_frac90.RData")
#bayescpi_frac90 <- bayescpi_frac90$xval.result
#bayescpi_frac90$r.mean <- as.numeric(bayescpi_frac90$r.mean)

# BayesCpi - K-fold
#Load cross-validation results for BayesCpi model
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

#Combine all model results into a single list
model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

#Filter for specific traits of interest
all_bio <- all_models[all_models$trait %in% c('bio1', 'bio2', 'bio3', 'bio4', 'bio5','bio6',
                                              'bio7','bio8','bio9','bio10','bio11',
                                              'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19'),]

#################################################################################################################
#Plot
#SUP FIGURE 1A
all_bio$trait <- factor(all_bio$trait, levels = paste0("bio", 1:19))
ggplot(all_bio, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1)

##########################
# Genomic selection (n=31 core)
##########################
library(rrBLUP)
library(dplyr)

setwd("~/R/barley_collab/parental_selection/cross_validation/")
envdat <- read.csv('Supplemental_dataset_5_worldclim_training31_copy.csv', head = TRUE) # full environmental dataset
gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = TRUE) #this it the -1/0/1 format 

# Set row names for genotype data
row.names(gd1) <- gd1$taxa
gd3 <- gd1[,-1]  # Remove 'taxa' column for analysis
g.in <- as.matrix(gd3)

# RR-BLUP
row.names(envdat) <- envdat$Taxa
y.in.rr <- envdat[,5:ncol(envdat)]  # Data starts at column 5 for Core 31

# Subset to core lines
minicore_entries <- which(envdat$Core == TRUE)
y.in.rr <- y.in.rr[minicore_entries,]
y.in.rr <- y.in.rr[,-1]  # Remove 'Core' column
y.in.mat <- as.matrix(y.in.rr)

# Set up training and genotype data
train <- row.names(y.in.mat) # Names of training lines
g.train <- g.in[train,]  # Set training genotypes

# List of traits to analyze
traits <- colnames(y.in.mat)

# Prediction setup
pred <- setdiff(row.names(g.in), train)  # Names of lines not in training set
g.pred <- g.in[pred,]

# Initialize objects for storing results
# Initialize objects for storing results
marker.list <- list()
gebv_df <- data.frame(matrix(nrow = nrow(g.in), ncol = length(traits)))  # Initialize with all lines (train + pred)
colnames(gebv_df) <- traits

# RR-BLUP loop for each trait
for(t in 1:length(traits)) {
  trait <- traits[t]
  y.train <- as.matrix(y.in.mat[train, trait])  # Training set for the trait
  
  # Run RR-BLUP model
  solve.out <- mixed.solve(y = y.train, Z = g.train, SE = FALSE, return.Hinv = FALSE)
  u.hat <- solve.out$u
  
  # Calculate GEBVs for both prediction and training set
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  
  # Use match to ensure correct row alignment
  pred_rows <- match(row.names(g.pred), row.names(g.in))  # Indices for pred rows
  train_rows <- match(row.names(g.train), row.names(g.in))  # Indices for train rows
  
  # Store the results in the combined dataframe
  gebv_df[pred_rows, t] <- GEBV  # Predictions for test lines
  gebv_df[train_rows, t] <- GEBV_train  # Predictions for training lines
}

# Set row names for the result to match the line names from `g.in`
row.names(gebv_df) <- row.names(g.in)

write.csv(gebv_df, 'rrblup_GEBV_data_barley.csv')
#saveRDS(marker.list, 'rrblup_markereffects_data_barley.csv')

###############################################
#core 100 
###############################################
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

#Local
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
all_bio$trait <- factor(all_bio$trait, levels = paste0("bio", 1:19))
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

#Trait column to a factor to order it
all_bio$trait <- factor(all_bio$trait, levels = paste0("bio", 1:19))
# Plot the results with ordered traits
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

##########################
# Genomic selection (n=100 core)
##########################
library(rrBLUP)
library(dplyr)

# Set the working directory for the 100-core dataset
setwd("~/R/barley_collab/parental_selection/new_core_100/")

# Load the full environmental dataset and genotype data
envdat <- read.csv('Supplemental_dataset_5_worldclim_training100.csv', head = TRUE)  # Full environmental dataset for 100-core
gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = TRUE)  # Genotype data

# Set row names for genotype data
row.names(gd1) <- gd1$taxa
gd3 <- gd1[,-1]  # Remove 'taxa' column for analysis
g.in <- as.matrix(gd3)

# RR-BLUP setup
row.names(envdat) <- envdat$Taxa
y.in.rr <- envdat[,5:ncol(envdat)]  # Data starts at column 5 for Core 100

# Subset to core lines
minicore_entries <- which(envdat$Core == TRUE)
y.in.rr <- y.in.rr[minicore_entries,]
y.in.rr <- y.in.rr[,-1]  # Remove 'Core' column
y.in.mat <- as.matrix(y.in.rr)

# Set up training and genotype data
train <- row.names(y.in.mat)  # Names of training lines
g.train <- g.in[train,]  # Set training genotypes

# List of traits to analyze
traits <- colnames(y.in.mat)

# Prediction setup
pred <- setdiff(row.names(g.in), train)  # Names of lines not in training set
g.pred <- g.in[pred,]

# Initialize objects for storing results
marker.list <- list()
gebv_df <- data.frame(matrix(nrow = nrow(g.in), ncol = length(traits)))  # Initialize with all lines (train + pred)
colnames(gebv_df) <- traits

# RR-BLUP loop for each trait
for(t in 1:length(traits)) {
  trait <- traits[t]
  y.train <- as.matrix(y.in.mat[train, trait])  # Training set for the trait
  
  # Run RR-BLUP model
  solve.out <- mixed.solve(y = y.train, Z = g.train, SE = FALSE, return.Hinv = FALSE)
  u.hat <- solve.out$u
  
  # Calculate GEBVs for both prediction and training set
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  
  # Use match to ensure correct row alignment
  pred_rows <- match(row.names(g.pred), row.names(g.in))  # Indices for pred rows
  train_rows <- match(row.names(g.train), row.names(g.in))  # Indices for train rows
  
  # Store the results in the combined dataframe
  gebv_df[pred_rows, t] <- GEBV  # Predictions for test lines
  gebv_df[train_rows, t] <- GEBV_train  # Predictions for training lines
}

# Set row names for the result to match the line names from `g.in`
row.names(gebv_df) <- row.names(g.in)

write.csv(gebv_df, 'rrblup_GEBV_data_barley_100.csv')
#saveRDS(marker.list, 'rrblup_markereffects_data_barley_100.csv')


#######################################################################################################
#SUPPLEMENTAL FIGURE 3
#######################################################################################################
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
  pivot_longer(cols = c(bio1_range), names_to = "Range", values_to = "RangeValue") %>%
  
  # Recode HCPC_group to rename 'Group 1' to 'Population 1', etc.
  mutate(HCPC_group = recode(HCPC_group,
                             "Group 1" = "Population 1",
                             "Group 2" = "Population 2",
                             "Group 3" = "Population 3",
                             "Group 4" = "Population 4",
                             "Group 5" = "Population 5"))

# Display the first few rows of the reshaped and recoded dataset to check
head(long_data)

# Plot the data
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
setwd("~/R/barley_collab/parental_selection/core_comparisons/")
data <- read.csv('comparing_cores_geav.csv')

# Reshape the data from wide to long format
long_data <- data %>%
  pivot_longer(cols = c(bio3_geav_n31, bio3_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio3_range), names_to = "Range", values_to = "RangeValue") %>%
  
  # Recode HCPC_group to rename 'Group 1' to 'Population 1', etc.
  mutate(HCPC_group = recode(HCPC_group,
                             "Group 1" = "Population 1",
                             "Group 2" = "Population 2",
                             "Group 3" = "Population 3",
                             "Group 4" = "Population 4",
                             "Group 5" = "Population 5"))

# Plot the data
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
# Reshape the data for BIO 4
long_data <- data %>%
  pivot_longer(cols = c(bio4_geav_n31, bio4_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio4_range), names_to = "Range", values_to = "RangeValue") %>%
  mutate(HCPC_group = recode(HCPC_group,
                             "Group 1" = "Population 1",
                             "Group 2" = "Population 2",
                             "Group 3" = "Population 3",
                             "Group 4" = "Population 4",
                             "Group 5" = "Population 5"))

# Plot for BIO 4
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
long_data <- data %>%
  pivot_longer(cols = c(bio6_geav_n31, bio6_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio6_range), names_to = "Range", values_to = "RangeValue") %>%
  mutate(HCPC_group = recode(HCPC_group,
                             "Group 1" = "Population 1",
                             "Group 2" = "Population 2",
                             "Group 3" = "Population 3",
                             "Group 4" = "Population 4",
                             "Group 5" = "Population 5"))

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
# Reshape the data for BIO 11
long_data <- data %>%
  pivot_longer(cols = c(bio11_geav_n31, bio11_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio11_range), names_to = "Range", values_to = "RangeValue") %>%
  mutate(HCPC_group = recode(HCPC_group,
                             "Group 1" = "Population 1",
                             "Group 2" = "Population 2",
                             "Group 3" = "Population 3",
                             "Group 4" = "Population 4",
                             "Group 5" = "Population 5"))

# Plot for BIO 11
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
long_data <- data %>%
  pivot_longer(cols = c(bio14_geav_n31, bio14_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio14_range), names_to = "Range", values_to = "RangeValue") %>%
  mutate(HCPC_group = recode(HCPC_group,
                             "Group 1" = "Population 1",
                             "Group 2" = "Population 2",
                             "Group 3" = "Population 3",
                             "Group 4" = "Population 4",
                             "Group 5" = "Population 5"))

# Plot for BIO 14
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
long_data <- data %>%
  pivot_longer(cols = c(bio17_geav_n31, bio17_geav_n100), names_to = "Core", values_to = "GEBV") %>%
  pivot_longer(cols = c(bio17_range), names_to = "Range", values_to = "RangeValue") %>%
  mutate(HCPC_group = recode(HCPC_group,
                             "Group 1" = "Population 1",
                             "Group 2" = "Population 2",
                             "Group 3" = "Population 3",
                             "Group 4" = "Population 4",
                             "Group 5" = "Population 5"))

# Plot for BIO 17
ggplot(long_data, aes(x = RangeValue, y = GEBV, color = Core)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Comparison of GEAV between Core Collections for bio17",
       x = "bio17_range",
       y = "GEAV") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ HCPC_group, scales = "fixed")


#######################################################################################################
#SUPPLEMENTAL FIGURE 4
#######################################################################################################
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

#######################################################################################################
#SUPPLEMENTAL FIGURE 5
#######################################################################################################
# Load necessary libraries
library(ggplot2)
library(factoextra)

###########
#All populations
###########
setwd("~/R/barley_collab/parental_selection/cross_validation/")
data <- read.csv("n784_climate_data_n7.csv")
env_variables <- data[, 9:15]

# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)

# Combine PCA scores with grouping variable
pca_data <- data.frame(pca_result$x, HCPC_group = data$HCPC_group)

# Create a biplot with color by 'HCPC_group'
p <- fviz_pca_biplot(pca_result, 
                     geom = "point", 
                     col.ind = "black", 
                     col.var = "blue", 
                     addEllipses = FALSE, 
                     habillage = pca_data$HCPC_group) +
  ggtitle("PCA Biplot All Populations") +  # Title
  theme_classic() +  
  
  # Adjust text sizes using theme()
  theme(
    plot.title = element_text(size = 18),     # Increase title size
    axis.title = element_text(size = 14),     # Increase axis title size
    axis.text = element_text(size = 12),      # Increase axis tick label size
    legend.text = element_text(size = 12),    # Increase legend text size
    legend.title = element_text(size = 14),   # Increase legend title size
    axis.line = element_blank()               # Remove x and y axis lines
  )

# Adjust x and y axis limits
p1 <- p + xlim(-8, 8) + ylim(-8, 8)

# Display the plot
print(p1)
ggsave("~/R/barley_collab/parental_selection/figure_pannels/pca_biplot_Figure2F.pdf", plot = p1, width = 8, height = 6)

###########
#Pop 1
###########
data <- read.csv("biodata_HCPC_Group_1.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)

# Create a biplot and adjust text size
biplot <- fviz_pca_biplot(pca_result, 
                          geom = "point", 
                          col.ind = "black", 
                          col.var = "blue", 
                          addEllipses = TRUE, 
                          labelsize = 5) + # Increase the size of the variable and individual labels
  ggtitle("PCA Biplot - Population 1") +
  
  # Use theme_classic and remove axis lines
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),  # Increase plot title size
    axis.title = element_text(size = 14),  # Increase axis titles size
    axis.text = element_text(size = 12),   # Increase axis tick labels size
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    axis.line = element_blank(),           # Remove axis lines
    panel.grid = element_blank()            # Remove all grid lines
  )

# Optionally set limits on x and y axes
p2 <- biplot + xlim(-8, 8) + ylim(-8, 8)

# Display the plot
print(p2)

# Save the plot as a PDF
ggsave("~/R/barley_collab/parental_selection/figure_pannels/pca_biplot_Figure2A.pdf", 
       plot = p2, width = 8, height = 6)

###########
#Pop 2
###########
data <- read.csv("biodata_HCPC_Group_2.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)

biplot <- fviz_pca_biplot(pca_result, 
                          geom = "point", 
                          col.ind = "black", 
                          col.var = "blue", 
                          addEllipses = TRUE, 
                          labelsize = 5) + # Increase the size of the variable and individual labels
  ggtitle("PCA Biplot - Population 2") +
  
  # Use theme_classic and remove axis lines
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),  # Increase plot title size
    axis.title = element_text(size = 14),  # Increase axis titles size
    axis.text = element_text(size = 12),   # Increase axis tick labels size
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    axis.line = element_blank(),           # Remove axis lines
    panel.grid = element_blank()            # Remove all grid lines
  )

# Optionally set limits on x and y axes
p3 <- biplot + xlim(-8, 8) + ylim(-8, 8)

# Display the plot
print(p3)

# Save the plot as a PDF
ggsave("~/R/barley_collab/parental_selection/figure_pannels/pca_biplot_Figure2B.pdf", 
       plot = p3, width = 8, height = 6)

###########
#Pop 3
###########
data <- read.csv("biodata_HCPC_Group_3.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)

biplot <- fviz_pca_biplot(pca_result, 
                          geom = "point", 
                          col.ind = "black", 
                          col.var = "blue", 
                          addEllipses = TRUE, 
                          labelsize = 5) + # Increase the size of the variable and individual labels
  ggtitle("PCA Biplot - Population 3") +
  
  # Use theme_classic and remove axis lines
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),  # Increase plot title size
    axis.title = element_text(size = 14),  # Increase axis titles size
    axis.text = element_text(size = 12),   # Increase axis tick labels size
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    axis.line = element_blank(),           # Remove axis lines
    panel.grid = element_blank()            # Remove all grid lines
  )

# Optionally set limits on x and y axes
p4 <- biplot + xlim(-8, 8) + ylim(-8, 8)

# Display the plot
print(p4)

# Save the plot as a PDF
ggsave("~/R/barley_collab/parental_selection/figure_pannels/pca_biplot_Figure2C.pdf", 
       plot = p4, width = 8, height = 6)

###########
#Pop 4
###########
data <- read.csv("biodata_HCPC_Group_4.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)
# Create a biplot
biplot <- fviz_pca_biplot(pca_result, 
                          geom = "point", 
                          col.ind = "black", 
                          col.var = "blue", 
                          addEllipses = TRUE, 
                          labelsize = 5) + # Increase the size of the variable and individual labels
  ggtitle("PCA Biplot - Population 4") +
  
  # Use theme_classic and remove axis lines
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),  # Increase plot title size
    axis.title = element_text(size = 14),  # Increase axis titles size
    axis.text = element_text(size = 12),   # Increase axis tick labels size
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    axis.line = element_blank(),           # Remove axis lines
    panel.grid = element_blank()            # Remove all grid lines
  )

# Optionally set limits on x and y axes
p5 <- biplot + xlim(-8, 8) + ylim(-8, 8)

# Display the plot
print(p5)

# Save the plot as a PDF
ggsave("~/R/barley_collab/parental_selection/figure_pannels/pca_biplot_Figure2D.pdf", 
       plot = p5, width = 8, height = 6)

###########
#Pop 5
###########
data <- read.csv("biodata_HCPC_Group_5.csv")
env_variables <- data[, 9:15]
#env_variables <- env_data[, -1]
# Perform PCA
pca_result <- prcomp(env_variables, scale. = TRUE)
# Create a biplot
biplot <- fviz_pca_biplot(pca_result, 
                          geom = "point", 
                          col.ind = "black", 
                          col.var = "blue", 
                          addEllipses = TRUE, 
                          labelsize = 5) + # Increase the size of the variable and individual labels
  ggtitle("PCA Biplot - Population 5") +
  
  # Use theme_classic and remove axis lines
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),  # Increase plot title size
    axis.title = element_text(size = 14),  # Increase axis titles size
    axis.text = element_text(size = 12),   # Increase axis tick labels size
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    axis.line = element_blank(),           # Remove axis lines
    panel.grid = element_blank()            # Remove all grid lines
  )

# Optionally set limits on x and y axes
p6 <- biplot + xlim(-8, 8) + ylim(-8, 8)

# Display the plot
print(p6)

# Save the plot as a PDF
ggsave("~/R/barley_collab/parental_selection/figure_pannels/pca_biplot_Figure2E.pdf", 
       plot = p6, width = 8, height = 6)


#######################################################################################################
#SUPPLEMENTAL FIGURE 6
#######################################################################################################
#SPECIES DISTRIBUTION 
#########################################
#Example of SDM modelling with Population 1 - iterate through various Populations
#Run on HPC

# Set CRAN mirror for installing packages
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Loading necessary libraries for spatial analysis, raster data, and plotting
library(dismo)         # For species distribution modeling and MaxEnt
library(raster)        # For working with raster data
library(maptools)      # For working with spatial objects
library(rasterVis)     # For visualizing raster data
data(wrld_simpl)       # World map data (simplified version for plotting)
library(readr)         # For reading CSV files

# Load your coordinate data (longitude/latitude) for the barley population
can <- read.csv("/home/ahmccorm/kantar_koastore/anna/Barley_collab/barley_parental/SDM_2/pop1.csv", header=T, sep=",") 

# Filter out rows with missing longitude and latitude values
can <- can[complete.cases(can$lon), ]  # Remove rows with missing longitude
can <- can[complete.cases(can$lat), ]  # Remove rows with missing latitude
summary(can)

# Remove duplicate coordinates (same lat/lon)
dups <- duplicated(can[, c("lon", "lat")])
sum(dups)  # Check how many duplicates exist
can <- can[!dups, ]  # Keep only unique records
can <- can[, 2:3]  # Keep only longitude and latitude columns

# Plot world map and overlay barley population points
plot(wrld_simpl, axes=TRUE, col="light yellow")  # Plot world map
box()  # Add a border around the plot
points(can$lon, can$lat, col='orange', pch=20, cex=0.75)  # Add barley population points

# Load climate data (WorldClim layers) - .tif files for climate variables
files <- list.files(path="~/kantar_koastore/anna/Barley_collab/barley_parental/SDM/wc2.0_30s_bio/", pattern='tif', full.names=TRUE)
predictors <- stack(files)  # Stack raster layers for multiple climate variables

# Plot the first layer of the climate data raster stack
plot(predictors, 1)  # Plot the first climate variable
plot(wrld_simpl, add=TRUE)  # Overlay the world map
points(can$lon, can$lat, col='blue')  # Overlay population points in blue

# Extract climate data (predictor values) for the barley population locations
presvals <- extract(predictors, can)

# Generate random background points for species distribution modeling (absence data)
set.seed(0)
backgr <- randomPoints(predictors, 500)  # Generate 500 random background points
absvals <- extract(predictors, backgr)  # Extract climate data for background points

# Combine presence and absence data into a single data frame
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))  # Presence-absence labels
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))  # Combine data
head(sdmdata)

# Drop a specific layer (e.g., 'biome') from the raster stack if needed
pred_nf <- dropLayer(predictors, 'biome')

# Split the data into training (2/3) and testing (1/3) using k-fold cross-validation
group <- kfold(can, 3)  # Split data into 3 groups
pres_train <- can[group != 1, ]  # 2/3 for training
pres_test <- can[group == 1, ]  # 1/3 for testing

# Generate random background points for training and testing
backg <- randomPoints(pred_nf, n=1000, extf=1.25)  # Generate 1000 random points
colnames(backg) = c('lon', 'lat')  # Set column names

# Split background data into training (2/3) and testing (1/3)
group <- kfold(backg, 3)
backg_train <- backg[group != 1, ]  # 2/3 for training
backg_test <- backg[group == 1, ]  # 1/3 for testing

# Maxent Model Setup and Execution
library(rJava)  # Required for running Maxent
install.packages("snow")  # For parallel computing with Maxent (optional)

# Check if Maxent .jar file exists
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

# Run Maxent model if the .jar file is available
if (file.exists(jar)) {
  xm <- maxent(predictors, pres_train)  # Fit Maxent model
  plot(xm)  # Plot the Maxent model (variable contribution)
} else {
  cat('cannot run this example because maxent is not available')
  plot(1)  # Dummy plot if Maxent is unavailable
}

# World suitability map based on Maxent model
if (file.exists(jar)) {
  response(xm)  # Generate response curves for environmental variables
} else {
  cat('cannot run this example because maxent is not available')
  plot(1)  # Dummy plot if Maxent is unavailable
}

# Evaluate the Maxent model using the test data and generate predictions
if (file.exists(jar)) {
  xme <- evaluate(pres_test, backg_test, xm, predictors)  # Evaluate model
  px <- predict(predictors, xm, progress='')  # Predict suitability across the world
  trxm <- threshold(xme, 'spec_sens')  # Set presence/absence threshold
  
  par(mfrow=c(1,2))  # Set up for side-by-side plots
  plot(px, main='Maxent, raw values')  # Plot raw suitability values
  plot(wrld_simpl, add=TRUE, border='dark grey')  # Overlay world map
  
  plot(px > trxm, main='presence/absence')  # Plot presence/absence based on threshold
  plot(wrld_simpl, add=TRUE, border='dark grey')  # Overlay world map
  points(pres_train, pch='+')  # Add training points
} else {
  plot(1)  # Dummy plot if Maxent is unavailable
}

# Save the suitability map as a .tif file
writeRaster(px, '~/kantar_koastore/anna/Barley_collab/barley_parental/SDM_2/SDM_pop1.tif', options=c('TFW=YES'))

# Plot ROC curve (AUC) for model performance
par(mfrow=c(1, 1))  # Reset to single plot layout
plot(xme, 'ROC', cex.lab = 1, cex.axis =1, cex.main =2)  # Plot ROC curve

# Plot variable contribution (Maxent model output)
plot(xm, cex.lab = 2, cex.axis = 4, cex.main = 2, cex.sub = 3)

# Save the entire workspace for later use
save.image(file="~/kantar_koastore/anna/Barley_collab/barley_parental/SDM_2/my_workspace_pop1.RData")

#######################################################################################################
#SUPPLEMENTAL FIGURE 6
#######################################################################################################
#Pop 1
###########
library(ggplot2)
library(viridis)
library(raster)
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

#######################################################################################################
#SUPPLEMENTAL FIGURE 7
#######################################################################################################
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

#######################################################################################################
#SUPPLEMENTAL FIGURE 8
#######################################################################################################
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

#######################################################################################################
#SUPPLEMENTAL FIGURE 9
#######################################################################################################
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

#######################################################################################################
#SUPPLEMENTAL FIGURE 10
#######################################################################################################
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



#######################################################################################################
#SUPPLEMENTAL FIGURE 11
#######################################################################################################
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

all_bio <- all_bio %>%
  mutate(trait = factor(trait, levels = c('scores_pop1', 'scores_pop2', 'scores_pop3', 'scores_pop4', 'scores_pop5'),
                        labels = c('Population 1', 'Population 2', 'Population 3', 'Population 4', 'Population 5')))

ggplot(all_bio, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color = "red", size = 0.5, linetype = "longdash") +
  
  # Increase text sizes
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),  # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        axis.title.x = element_text(size = 16), # Increase x-axis title size
        axis.title.y = element_text(size = 16), # Increase y-axis title size
        strip.text = element_text(size = 14),   # Increase facet label size
        legend.text = element_text(size = 12),  # Increase legend text size
        legend.title = element_text(size = 14)) + # Increase legend title size
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
  
  # Calculate GEBVs for test and training
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
#SUPPLEMENTAL FIGURE 11B
#########################################
#Plotting results by population of GEBVs
#########################################
setwd("/Users/annamccormick/R/barley_collab/parental_selection/cross_validation_SDMtrial/")
data2 <-read.csv('rrblup_GEBV_data_barley_31_plots.csv')

library(ggplot2)
# Pop1
a <- ggplot(data2, aes(x = HCPC_group, y = scores_pop1, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Population 1" = "red", "Population 2" = "yellow", "Population 3" = "green", "Population 4" = "blue", "Population 5" = "purple")) +
  labs(title = "Population 1 raster",
       x = "HCPC Group",
       y = "GEAV") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
a


# Pop2
b <-ggplot(data2, aes(x = HCPC_group, y = scores_pop2, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Population 1" = "red", "Population 2" = "yellow", "Population 3" = "green", "Population 4" = "blue", "Population 5" = "purple")) +
  labs(title = "Population 2 raster",
       x = "HCPC Group",
       y = "GEAV") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
b
# Pop3
c<- ggplot(data2, aes(x = HCPC_group, y = scores_pop3, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Population 1" = "red", "Population 2" = "yellow", "Population 3" = "green", "Population 4" = "blue", "Population 5" = "purple")) +
  labs(title = "Population 3 raster",
       x = "HCPC Group",
       y = "GEAV") +
  theme_minimal() +
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
c
# Pop4
d<- ggplot(data2, aes(x = HCPC_group, y = scores_pop4, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Population 1" = "red", "Population 2" = "yellow", "Population 3" = "green", "Population 4" = "blue", "Population 5" = "purple")) +
  labs(title = "Population 4 raster",
       x = "HCPC Group",
       y = "GEAV") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
d
# Pop5
e<-ggplot(data2, aes(x = HCPC_group, y = scores_pop5, fill = HCPC_group)) +
  geom_boxplot(color = "black") +
  scale_color_manual(values = c("Population 1" = "red", "Population 2" = "yellow", "Population 3" = "green", "Population 4" = "blue", "Population 5" = "purple")) +
  labs(title = "Population 5 raster",
       x = "HCPC Group",
       y = "GEAV") +
  theme_minimal()+
  geom_hline(yintercept = 0, color="red", size = 0.5, linetype = "longdash")
e

library(gridExtra)
grid.arrange(a, b, c, d, e, ncol=3) #save pdf 15x20

#######################################################################################################
#SUPPLEMENTAL FIGURE 12
#######################################################################################################
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
setwd("~/R/barley_collab/parental_selection/cross_validation/")
barley_pca <- read.csv("barley_pca_copy.csv")
individuals_in_groups <- split(barley_pca$sample.id, barley_pca$HCPC_group)
for(i in seq_along(individuals_in_groups)) {
  assign(paste0("Population", i), individuals_in_groups[[i]])
}

names(individuals_in_groups)

# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group1_column_names <- c("chrom", "position","rs", individuals_in_groups[["Population 1"]])  # Assuming these are valid column names

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
  assign(paste0("Population", i), individuals_in_groups[[i]])
}

names(individuals_in_groups)

# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group2_column_names <- c("chrom", "position","rs", individuals_in_groups[["Population 2"]])  # Assuming these are valid column names

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
  assign(paste0("Population", i), individuals_in_groups[[i]])
}


names(individuals_in_groups)
# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group3_column_names <- c("chrom", "position","rs", individuals_in_groups[["Population 3"]])  # Assuming these are valid column names

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
  assign(paste0("Population", i), individuals_in_groups[[i]])
}

names(individuals_in_groups)
# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group4_column_names <- c("chrom", "position","rs", individuals_in_groups[["Population 4"]])  # Assuming these are valid column names

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
  assign(paste0("Population", i), individuals_in_groups[[i]])
}

names(individuals_in_groups)
# First, ensure group1_column_names only contains column names that exist in bio1_Marker_Effects_chrom_1H
group5_column_names <- c("chrom", "position","rs", individuals_in_groups[["Population 5"]])  # Assuming these are valid column names

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

#######################################################################################################
#SUPPLEMENTAL FIGURE 13
#######################################################################################################
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

#######################################################################################################
#SUPPLEMENTAL FIGURE 14
#######################################################################################################
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


#######################################################################################################
#SUPPLEMENTAL FIGURE 15
#######################################################################################################
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

#########################################
#Bio 1
#########################################

# Loop through each chromosome
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[2]  # Variable of interest - bio1
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



# Print all column names of the dat dataframe - to check for what variable to set
#print(colnames(dat))

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

#################
# Define the list of SNPs of interest
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter the combined data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the mean effect size for each of these SNPs
filtered_snps_summary <- filtered_snps %>%
  group_by(rs) %>%
  summarise(mean_effect = mean(mean_effect, na.rm = TRUE))

# Display the mean effect sizes for the SNPs of interest
print(filtered_snps_summary)


percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

print(percentile_95)
print(percentile_5)

summary(combined_data$mean_effect)

# Plot with both percentiles
iqr <- IQR(combined_data$mean_effect)
n <- length(combined_data$mean_effect)
bin_width <- 2 * iqr / (n^(1/3))

bin_width

ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.0018, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +  # Adjust `y` as needed
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +   # Adjust `y` as needed
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 1",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))


# Define the list of SNPs you care about
library(ggplot2)
library(dplyr)
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter combined_data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the 5th and 95th percentiles
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

# Add a column to indicate whether each SNP falls outside the 5th-95th percentile range
filtered_snps <- filtered_snps %>%
  mutate(outside_percentiles = ifelse(mean_effect < percentile_5 | mean_effect > percentile_95, TRUE, FALSE))

# Print the filtered SNPs and their status
print(filtered_snps)

# Plot the histogram with highlighted SNPs
# Plot the histogram with highlighted SNPs and their names
ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.0018, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +
  geom_point(data = filtered_snps %>% filter(outside_percentiles == TRUE), 
             aes(x = mean_effect, y = 50), color = "blue", size = 3) + # Highlight points outside percentiles
  geom_text(data = filtered_snps %>% filter(outside_percentiles == TRUE),
            aes(x = mean_effect, y = 60, label = rs), color = "blue", angle = 90, hjust = -0.2, size = 4) + # Add SNP names above dots
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 1",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))


########################################
# Print all column names of the dat dataframe - to check for what variable to set
#print(colnames(dat))
########################################
#########################################
#Bio 3
#########################################
library(dplyr)
library(circlize)
setwd("~/R/barley_collab/parental_selection/cross_validation/")
dat <- read.csv("Barley_marker_effects_chromosome_map_expanded.csv", header = TRUE)
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

# Unique chromosomes
all_chroms <- unique(dat$chrom)
all_chroms <- all_chroms[all_chroms != "UN"]

# Initialize list to store outputs
all_chrom_data <- list()

#########################################
#Bio 3
#########################################

# Loop through each chromosome
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[3]  # Variable of interest - bio 3
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


# Print all column names of the dat dataframe - to check for what variable to set
#print(colnames(dat))

# Extract and sort the dataframe for Chromosome 1H
bio3_Marker_Effects_chrom_1H <- all_chrom_data[["bio3_Marker_Effects_chrom_1H"]]


bio3_Marker_Effects_chrom_1H <- bio3_Marker_Effects_chrom_1H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio3_Marker_Effects_chrom_2H <- all_chrom_data[["bio3_Marker_Effects_chrom_2H"]]


bio3_Marker_Effects_chrom_2H <- bio3_Marker_Effects_chrom_2H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio3_Marker_Effects_chrom_3H <- all_chrom_data[["bio3_Marker_Effects_chrom_3H"]]


bio3_Marker_Effects_chrom_3H <- bio3_Marker_Effects_chrom_3H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio3_Marker_Effects_chrom_4H <- all_chrom_data[["bio3_Marker_Effects_chrom_4H"]]

bio3_Marker_Effects_chrom_4H <- bio3_Marker_Effects_chrom_4H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio3_Marker_Effects_chrom_5H <- all_chrom_data[["bio3_Marker_Effects_chrom_5H"]]

bio3_Marker_Effects_chrom_5H <- bio3_Marker_Effects_chrom_5H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio3_Marker_Effects_chrom_6H <- all_chrom_data[["bio3_Marker_Effects_chrom_6H"]]

bio3_Marker_Effects_chrom_6H <- bio3_Marker_Effects_chrom_6H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio3_Marker_Effects_chrom_7H <- all_chrom_data[["bio3_Marker_Effects_chrom_7H"]]

bio3_Marker_Effects_chrom_7H <- bio3_Marker_Effects_chrom_7H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

###### ALL CHROM
library(dplyr)
library(purrr)

# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(all_chrom_data)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Define the list of SNPs of interest
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter the combined data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the mean effect size for each of these SNPs
filtered_snps_summary <- filtered_snps %>%
  group_by(rs) %>%
  summarise(mean_effect = mean(mean_effect, na.rm = TRUE))

print(filtered_snps_summary)


percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

print(percentile_95)
print(percentile_5)

summary(combined_data$mean_effect)

# Plot with both percentiles
iqr <- IQR(combined_data$mean_effect)
n <- length(combined_data$mean_effect)
bin_width <- 2 * iqr / (n^(1/3))
bin_width

ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.00064, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +  # Adjust `y` as needed
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +   # Adjust `y` as needed
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 3",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))

library(ggplot2)
library(dplyr)
# Define the list of SNPs you care about
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter combined_data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the 5th and 95th percentiles
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

# Add a column to indicate whether each SNP falls outside the 5th-95th percentile range
filtered_snps <- filtered_snps %>%
  mutate(outside_percentiles = ifelse(mean_effect < percentile_5 | mean_effect > percentile_95, TRUE, FALSE))

# Print the filtered SNPs and their status
print(filtered_snps)

# Plot the histogram with highlighted SNPs
# Plot the histogram with highlighted SNPs and their names
ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.00064, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +
  geom_point(data = filtered_snps %>% filter(outside_percentiles == TRUE), 
             aes(x = mean_effect, y = 50), color = "blue", size = 3) + # Highlight points outside percentiles
  geom_text(data = filtered_snps %>% filter(outside_percentiles == TRUE),
            aes(x = mean_effect, y = 60, label = rs), color = "blue", angle = 90, hjust = -0.2, size = 4) + # Add SNP names above dots
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 3",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))

########################################
# Print all column names of the dat dataframe - to check for what variable to set
#print(colnames(dat))
########################################
#########################################
#Bio 4
#########################################
library(dplyr)
library(circlize)
setwd("~/R/barley_collab/parental_selection/cross_validation/")
dat <- read.csv("Barley_marker_effects_chromosome_map_expanded.csv", header = TRUE)
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

# Unique chromosomes
all_chroms <- unique(dat$chrom)
all_chroms <- all_chroms[all_chroms != "UN"]

# Initialize list to store outputs
all_chrom_data <- list()

#########################################
#Bio 4
#########################################

# Loop through each chromosome
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[4]  # Variable of interest - bio 4
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


# Print all column names of the dat dataframe - to check for what variable to set
#print(colnames(dat))

# Extract and sort the dataframe for Chromosome 1H
bio4_Marker_Effects_chrom_1H <- all_chrom_data[["bio4_Marker_Effects_chrom_1H"]]


bio4_Marker_Effects_chrom_1H <- bio4_Marker_Effects_chrom_1H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio4_Marker_Effects_chrom_2H <- all_chrom_data[["bio4_Marker_Effects_chrom_2H"]]


bio4_Marker_Effects_chrom_2H <- bio4_Marker_Effects_chrom_2H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio4_Marker_Effects_chrom_3H <- all_chrom_data[["bio4_Marker_Effects_chrom_3H"]]


bio4_Marker_Effects_chrom_3H <- bio4_Marker_Effects_chrom_3H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio4_Marker_Effects_chrom_4H <- all_chrom_data[["bio4_Marker_Effects_chrom_4H"]]

bio4_Marker_Effects_chrom_4H <- bio4_Marker_Effects_chrom_4H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio4_Marker_Effects_chrom_5H <- all_chrom_data[["bio4_Marker_Effects_chrom_5H"]]

bio4_Marker_Effects_chrom_5H <- bio4_Marker_Effects_chrom_5H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio4_Marker_Effects_chrom_6H <- all_chrom_data[["bio4_Marker_Effects_chrom_6H"]]

bio4_Marker_Effects_chrom_6H <- bio4_Marker_Effects_chrom_6H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio4_Marker_Effects_chrom_7H <- all_chrom_data[["bio4_Marker_Effects_chrom_7H"]]

bio4_Marker_Effects_chrom_7H <- bio4_Marker_Effects_chrom_7H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

###### ALL CHROM
library(dplyr)
library(purrr)

# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(all_chrom_data)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Define the list of SNPs of interest
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter the combined data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the mean effect size for each of these SNPs
filtered_snps_summary <- filtered_snps %>%
  group_by(rs) %>%
  summarise(mean_effect = mean(mean_effect, na.rm = TRUE))

# Display the mean effect sizes for the SNPs of interest
print(filtered_snps_summary)




# Sort the data by mean_effect in descending order to get the largest effect sizes at the top
combined_sorted <- combined_data %>%
  arrange(desc(abs(mean_effect)))


# Select the top 10 SNPs with the highest mean effect
top_10_snps <- head(combined_sorted, 20)

# Select and display the relevant information for these top 10 SNPs
top_10_snps_selected <- top_10_snps %>%
  select(chrom, position, rs, var, mean_effect)

# Display the top 10 SNPs with selected information
print(top_10_snps_selected)

#####
# Filter the combined data to include only positive mean effects
positive_effects <- combined_data %>%
  filter(mean_effect > 0)

# Sort the filtered data by mean_effect in descending order
positive_sorted <- positive_effects %>%
  arrange(desc(mean_effect))

# Select the top 10 SNPs with the highest positive mean effect
top_10_positive_snps <- head(positive_sorted, 20)

# Select and display the relevant information for these top 10 SNPs
top_10_positive_snps_selected <- top_10_positive_snps %>%
  select(chrom, position, rs, var, mean_effect)

# Display the top 10 SNPs with selected information
print(top_10_positive_snps_selected)

########histogram of marker effects
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

print(percentile_95)
print(percentile_5)

summary(combined_data$mean_effect)

# Plot with both percentiles
iqr <- IQR(combined_data$mean_effect)
n <- length(combined_data$mean_effect)
bin_width <- 2 * iqr / (n^(1/3))

bin_width

ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.1162, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +  # Adjust `y` as needed
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +   # Adjust `y` as needed
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 4",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))

library(ggplot2)
library(dplyr)
# Define the list of SNPs you care about
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter combined_data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the 5th and 95th percentiles
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

# Add a column to indicate whether each SNP falls outside the 5th-95th percentile range
filtered_snps <- filtered_snps %>%
  mutate(outside_percentiles = ifelse(mean_effect < percentile_5 | mean_effect > percentile_95, TRUE, FALSE))

# Print the filtered SNPs and their status
print(filtered_snps)

# Plot the histogram with highlighted SNPs
# Plot the histogram with highlighted SNPs and their names
ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.1162, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +
  geom_point(data = filtered_snps %>% filter(outside_percentiles == TRUE), 
             aes(x = mean_effect, y = 50), color = "blue", size = 3) + # Highlight points outside percentiles
  geom_text(data = filtered_snps %>% filter(outside_percentiles == TRUE),
            aes(x = mean_effect, y = 60, label = rs), color = "blue", angle = 90, hjust = -0.2, size = 4) + # Add SNP names above dots
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 4",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))



#########################################
#Bio 6
#########################################
library(dplyr)
library(circlize)
setwd("~/R/barley_collab/parental_selection/cross_validation/")
dat <- read.csv("Barley_marker_effects_chromosome_map_expanded.csv", header = TRUE)
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

# Unique chromosomes
all_chroms <- unique(dat$chrom)
all_chroms <- all_chroms[all_chroms != "UN"]

# Initialize list to store outputs
all_chrom_data <- list()

#########################################
#Bio 6
#########################################

# Loop through each chromosome
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[5]  # Variable of interest - bio 6
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


# Print all column names of the dat dataframe - to check for what variable to set
#print(colnames(dat))

# Extract and sort the dataframe for Chromosome 1H
bio6_Marker_Effects_chrom_1H <- all_chrom_data[["bio6_Marker_Effects_chrom_1H"]]


bio6_Marker_Effects_chrom_1H <- bio6_Marker_Effects_chrom_1H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio6_Marker_Effects_chrom_2H <- all_chrom_data[["bio6_Marker_Effects_chrom_2H"]]


bio6_Marker_Effects_chrom_2H <- bio6_Marker_Effects_chrom_2H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio6_Marker_Effects_chrom_3H <- all_chrom_data[["bio6_Marker_Effects_chrom_3H"]]


bio6_Marker_Effects_chrom_3H <- bio6_Marker_Effects_chrom_3H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio6_Marker_Effects_chrom_4H <- all_chrom_data[["bio6_Marker_Effects_chrom_4H"]]

bio6_Marker_Effects_chrom_4H <- bio6_Marker_Effects_chrom_4H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio6_Marker_Effects_chrom_5H <- all_chrom_data[["bio6_Marker_Effects_chrom_5H"]]

bio6_Marker_Effects_chrom_5H <- bio6_Marker_Effects_chrom_5H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio6_Marker_Effects_chrom_6H <- all_chrom_data[["bio6_Marker_Effects_chrom_6H"]]

bio6_Marker_Effects_chrom_6H <- bio6_Marker_Effects_chrom_6H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio6_Marker_Effects_chrom_7H <- all_chrom_data[["bio6_Marker_Effects_chrom_7H"]]

bio6_Marker_Effects_chrom_7H <- bio6_Marker_Effects_chrom_7H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

###### ALL CHROM
library(dplyr)
library(purrr)

# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(all_chrom_data)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Define the list of SNPs of interest
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter the combined data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the mean effect size for each of these SNPs
filtered_snps_summary <- filtered_snps %>%
  group_by(rs) %>%
  summarise(mean_effect = mean(mean_effect, na.rm = TRUE))

# Display the mean effect sizes for the SNPs of interest
print(filtered_snps_summary)

########histogram of marker effects
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

print(percentile_95)
print(percentile_5)

summary(combined_data$mean_effect)

# Plot with both percentiles
iqr <- IQR(combined_data$mean_effect)
n <- length(combined_data$mean_effect)
bin_width <- 2 * iqr / (n^(1/3))

bin_width

ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.00435, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +  # Adjust `y` as needed
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +   # Adjust `y` as needed
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 6",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))



# Define the list of SNPs you care about
library(ggplot2)
library(dplyr)
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter combined_data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the 5th and 95th percentiles
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

# Add a column to indicate whether each SNP falls outside the 5th-95th percentile range
filtered_snps <- filtered_snps %>%
  mutate(outside_percentiles = ifelse(mean_effect < percentile_5 | mean_effect > percentile_95, TRUE, FALSE))

# Print the filtered SNPs and their status
print(filtered_snps)

# Plot the histogram with highlighted SNPs
# Plot the histogram with highlighted SNPs and their names
ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.00435, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +
  geom_point(data = filtered_snps %>% filter(outside_percentiles == TRUE), 
             aes(x = mean_effect, y = 50), color = "blue", size = 3) + # Highlight points outside percentiles
  geom_text(data = filtered_snps %>% filter(outside_percentiles == TRUE),
            aes(x = mean_effect, y = 60, label = rs), color = "blue", angle = 90, hjust = -0.2, size = 4) + # Add SNP names above dots
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 6",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))

#########################################
#Bio 11
#########################################
library(dplyr)
library(circlize)
setwd("~/R/barley_collab/parental_selection/cross_validation/")
dat <- read.csv("Barley_marker_effects_chromosome_map_expanded.csv", header = TRUE)
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

# Unique chromosomes
all_chroms <- unique(dat$chrom)
all_chroms <- all_chroms[all_chroms != "UN"]

# Initialize list to store outputs
all_chrom_data <- list()

#########################################
#Bio 11
#########################################

# Loop through each chromosome
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[6]  # Variable of interest - bio 11
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


# Print all column names of the dat dataframe - to check for what variable to set
#print(colnames(dat))

# Extract and sort the dataframe for Chromosome 1H
bio11_Marker_Effects_chrom_1H <- all_chrom_data[["bio11_Marker_Effects_chrom_1H"]]


bio11_Marker_Effects_chrom_1H <- bio11_Marker_Effects_chrom_1H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio11_Marker_Effects_chrom_2H <- all_chrom_data[["bio11_Marker_Effects_chrom_2H"]]


bio11_Marker_Effects_chrom_2H <- bio11_Marker_Effects_chrom_2H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio11_Marker_Effects_chrom_3H <- all_chrom_data[["bio11_Marker_Effects_chrom_3H"]]


bio11_Marker_Effects_chrom_3H <- bio11_Marker_Effects_chrom_3H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio11_Marker_Effects_chrom_4H <- all_chrom_data[["bio11_Marker_Effects_chrom_4H"]]

bio11_Marker_Effects_chrom_4H <- bio11_Marker_Effects_chrom_4H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio11_Marker_Effects_chrom_5H <- all_chrom_data[["bio11_Marker_Effects_chrom_5H"]]

bio11_Marker_Effects_chrom_5H <- bio11_Marker_Effects_chrom_5H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio11_Marker_Effects_chrom_6H <- all_chrom_data[["bio11_Marker_Effects_chrom_6H"]]

bio11_Marker_Effects_chrom_6H <- bio11_Marker_Effects_chrom_6H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio11_Marker_Effects_chrom_7H <- all_chrom_data[["bio11_Marker_Effects_chrom_7H"]]

bio11_Marker_Effects_chrom_7H <- bio11_Marker_Effects_chrom_7H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

###### ALL CHROM
library(dplyr)
library(purrr)

# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(all_chrom_data)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Define the list of SNPs of interest
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter the combined data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the mean effect size for each of these SNPs
filtered_snps_summary <- filtered_snps %>%
  group_by(rs) %>%
  summarise(mean_effect = mean(mean_effect, na.rm = TRUE))

# Display the mean effect sizes for the SNPs of interest
print(filtered_snps_summary)


########histogram of marker effects
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

print(percentile_95)
print(percentile_5)

summary(combined_data$mean_effect)

# Plot with both percentiles
iqr <- IQR(combined_data$mean_effect)
n <- length(combined_data$mean_effect)
bin_width <- 2 * iqr / (n^(1/3))

bin_width

ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.0038, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +  # Adjust `y` as needed
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +   # Adjust `y` as needed
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 11",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))



# Define the list of SNPs you care about
library(ggplot2)
library(dplyr)
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter combined_data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the 5th and 95th percentiles
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

# Add a column to indicate whether each SNP falls outside the 5th-95th percentile range
filtered_snps <- filtered_snps %>%
  mutate(outside_percentiles = ifelse(mean_effect < percentile_5 | mean_effect > percentile_95, TRUE, FALSE))

# Print the filtered SNPs and their status
print(filtered_snps)

# Plot the histogram with highlighted SNPs
# Plot the histogram with highlighted SNPs and their names
ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.0038, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +
  geom_point(data = filtered_snps %>% filter(outside_percentiles == TRUE), 
             aes(x = mean_effect, y = 50), color = "blue", size = 3) + # Highlight points outside percentiles
  geom_text(data = filtered_snps %>% filter(outside_percentiles == TRUE),
            aes(x = mean_effect, y = 60, label = rs), color = "blue", angle = 90, hjust = -0.2, size = 4) + # Add SNP names above dots
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 11",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))


#########################################
#Bio 14
#########################################
library(dplyr)
library(circlize)
setwd("~/R/barley_collab/parental_selection/cross_validation/")
dat <- read.csv("Barley_marker_effects_chromosome_map_expanded.csv", header = TRUE)
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

# Unique chromosomes
all_chroms <- unique(dat$chrom)
all_chroms <- all_chroms[all_chroms != "UN"]

# Initialize list to store outputs
all_chrom_data <- list()

#########################################
#Bio 14
#########################################

# Loop through each chromosome
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[7]  # Variable of interest - bio 11
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


# Print all column names of the dat dataframe - to check for what variable to set
#print(colnames(dat))

# Extract and sort the dataframe for Chromosome 1H
bio14_Marker_Effects_chrom_1H <- all_chrom_data[["bio14_Marker_Effects_chrom_1H"]]


bio14_Marker_Effects_chrom_1H <- bio14_Marker_Effects_chrom_1H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio14_Marker_Effects_chrom_2H <- all_chrom_data[["bio14_Marker_Effects_chrom_2H"]]


bio14_Marker_Effects_chrom_2H <- bio14_Marker_Effects_chrom_2H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio14_Marker_Effects_chrom_3H <- all_chrom_data[["bio14_Marker_Effects_chrom_3H"]]


bio14_Marker_Effects_chrom_3H <- bio14_Marker_Effects_chrom_3H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio14_Marker_Effects_chrom_4H <- all_chrom_data[["bio14_Marker_Effects_chrom_4H"]]

bio14_Marker_Effects_chrom_4H <- bio14_Marker_Effects_chrom_4H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio14_Marker_Effects_chrom_5H <- all_chrom_data[["bio14_Marker_Effects_chrom_5H"]]

bio14_Marker_Effects_chrom_5H <- bio14_Marker_Effects_chrom_5H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio14_Marker_Effects_chrom_6H <- all_chrom_data[["bio14_Marker_Effects_chrom_6H"]]

bio14_Marker_Effects_chrom_6H <- bio14_Marker_Effects_chrom_6H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio14_Marker_Effects_chrom_7H <- all_chrom_data[["bio14_Marker_Effects_chrom_7H"]]

bio14_Marker_Effects_chrom_7H <- bio14_Marker_Effects_chrom_7H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

###### ALL CHROM
library(dplyr)
library(purrr)

# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(all_chrom_data)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Define the list of SNPs of interest
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter the combined data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the mean effect size for each of these SNPs
filtered_snps_summary <- filtered_snps %>%
  group_by(rs) %>%
  summarise(mean_effect = mean(mean_effect, na.rm = TRUE))

# Display the mean effect sizes for the SNPs of interest
print(filtered_snps_summary)


########histogram of marker effects
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

print(percentile_95)
print(percentile_5)

summary(combined_data$mean_effect)

# Plot with both percentiles
iqr <- IQR(combined_data$mean_effect)
n <- length(combined_data$mean_effect)
bin_width <- 2 * iqr / (n^(1/3))

bin_width

ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.0037, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +  # Adjust `y` as needed
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +   # Adjust `y` as needed
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 14",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))


# Define the list of SNPs you care about
library(ggplot2)
library(dplyr)
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter combined_data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the 5th and 95th percentiles
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

# Add a column to indicate whether each SNP falls outside the 5th-95th percentile range
filtered_snps <- filtered_snps %>%
  mutate(outside_percentiles = ifelse(mean_effect < percentile_5 | mean_effect > percentile_95, TRUE, FALSE))

# Print the filtered SNPs and their status
print(filtered_snps)

# Plot the histogram with highlighted SNPs
# Plot the histogram with highlighted SNPs and their names
ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.0037, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +
  geom_point(data = filtered_snps %>% filter(outside_percentiles == TRUE), 
             aes(x = mean_effect, y = 50), color = "blue", size = 3) + # Highlight points outside percentiles
  geom_text(data = filtered_snps %>% filter(outside_percentiles == TRUE),
            aes(x = mean_effect, y = 60, label = rs), color = "blue", angle = 90, hjust = -0.2, size = 4) + # Add SNP names above dots
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 14",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))


#########################################
#Bio 17
#########################################
library(dplyr)
library(circlize)
setwd("~/R/barley_collab/parental_selection/cross_validation/")
dat <- read.csv("Barley_marker_effects_chromosome_map_expanded.csv", header = TRUE)
mat <- readRDS("transposed_barley_genotype_matrix_for_chrommaps.rds")

# Unique chromosomes
all_chroms <- unique(dat$chrom)
all_chroms <- all_chroms[all_chroms != "UN"]

# Initialize list to store outputs
all_chrom_data <- list()

#########################################
#Bio 17
#########################################

# Loop through each chromosome
for (chrom_num in all_chroms) {
  mat_chrom <- mat[mat$chrom == chrom_num, ]
  dat_chrom <- dat[dat$chrom == chrom_num, ]
  
  comps <- colnames(dat)[8]  # Variable of interest - bio 11
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


# Print all column names of the dat dataframe - to check for what variable to set
#print(colnames(dat))

# Extract and sort the dataframe for Chromosome 1H
bio17_Marker_Effects_chrom_1H <- all_chrom_data[["bio17_Marker_Effects_chrom_1H"]]


bio17_Marker_Effects_chrom_1H <- bio17_Marker_Effects_chrom_1H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio17_Marker_Effects_chrom_2H <- all_chrom_data[["bio17_Marker_Effects_chrom_2H"]]


bio17_Marker_Effects_chrom_2H <- bio17_Marker_Effects_chrom_2H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio17_Marker_Effects_chrom_3H <- all_chrom_data[["bio17_Marker_Effects_chrom_3H"]]


bio17_Marker_Effects_chrom_3H <- bio17_Marker_Effects_chrom_3H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio17_Marker_Effects_chrom_4H <- all_chrom_data[["bio17_Marker_Effects_chrom_4H"]]

bio17_Marker_Effects_chrom_4H <- bio17_Marker_Effects_chrom_4H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio17_Marker_Effects_chrom_5H <- all_chrom_data[["bio17_Marker_Effects_chrom_5H"]]

bio17_Marker_Effects_chrom_5H <- bio17_Marker_Effects_chrom_5H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio17_Marker_Effects_chrom_6H <- all_chrom_data[["bio17_Marker_Effects_chrom_6H"]]

bio17_Marker_Effects_chrom_6H <- bio17_Marker_Effects_chrom_6H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

bio17_Marker_Effects_chrom_7H <- all_chrom_data[["bio17_Marker_Effects_chrom_7H"]]

bio17_Marker_Effects_chrom_7H <- bio17_Marker_Effects_chrom_7H %>%
  arrange(position)  # Use dplyr's arrange function to sort by position

###### ALL CHROM
library(dplyr)
library(purrr)

# Combine all chromosome data frames into one large data frame
combined_data <- bind_rows(all_chrom_data)

# Calculate the mean effect size for each SNP across all lines
combined_data$mean_effect <- rowMeans(combined_data[, -which(colnames(combined_data) %in% c("chrom", "position", "rs", "var"))], na.rm = TRUE)

# Define the list of SNPs of interest
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter the combined data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the mean effect size for each of these SNPs
filtered_snps_summary <- filtered_snps %>%
  group_by(rs) %>%
  summarise(mean_effect = mean(mean_effect, na.rm = TRUE))

# Display the mean effect sizes for the SNPs of interest
print(filtered_snps_summary)



# Check which SNPs from your list are missing in the dataset
missing_snps <- setdiff(snps_of_interest, filtered_snps$rs)
print(missing_snps)

# Remove trailing spaces and convert to lowercase for both
snps_of_interest <- tolower(trimws(snps_of_interest))
filtered_snps$rs <- tolower(trimws(filtered_snps$rs))


########histogram of marker effects
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

print(percentile_95)
print(percentile_5)

summary(combined_data$mean_effect)

# Plot with both percentiles
iqr <- IQR(combined_data$mean_effect)
n <- length(combined_data$mean_effect)
bin_width <- 2 * iqr / (n^(1/3))

bin_width

ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.0027, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +  # Adjust `y` as needed
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +   # Adjust `y` as needed
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 17",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))


# Define the list of SNPs you care about
library(ggplot2)
library(dplyr)
snps_of_interest <- c("SCRI_RS_154621", "12_30861", "12_30867", "12_30868", 
                      "SCRI_RS_172665", "12_30668", "12_30869", "BK_12", 
                      "BK_16", "BK_15", "12_30870", "12_30871", "BK_14", 
                      "12_30872", "BK_13", "SCRI_RS_142282", "11_10513", 
                      "SCRI_RS_137464", "12_30930", "BK_17", "11_20188", 
                      "12_30929", "12_30850", "12_30849", "BK_22", 
                      "11_11361", "11_10189", "11_20085", "11_11328", 
                      "12_30880")

# Filter combined_data for the SNPs of interest
filtered_snps <- combined_data %>%
  filter(rs %in% snps_of_interest)

# Calculate the 5th and 95th percentiles
percentile_95 <- quantile(combined_data$mean_effect, 0.95)
percentile_5 <- quantile(combined_data$mean_effect, 0.05)

# Add a column to indicate whether each SNP falls outside the 5th-95th percentile range
filtered_snps <- filtered_snps %>%
  mutate(outside_percentiles = ifelse(mean_effect < percentile_5 | mean_effect > percentile_95, TRUE, FALSE))

# Print the filtered SNPs and their status
print(filtered_snps)

# Plot the histogram with highlighted SNPs
# Plot the histogram with highlighted SNPs and their names
ggplot(combined_data, aes(x = mean_effect)) +
  geom_histogram(binwidth = 0.0027, color = "black", fill = "grey", alpha = 0.7) +
  geom_vline(aes(xintercept = percentile_95), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_vline(aes(xintercept = percentile_5), color = "red", linetype = "dashed", linewidth = 0.5) +  
  geom_text(aes(x = percentile_95, y = 200, label = round(percentile_95, 4)), color = "red", angle = 90, vjust = -1) +
  geom_text(aes(x = percentile_5, y = 200, label = round(percentile_5, 4)), color = "red", angle = 90, vjust = -1) +
  geom_point(data = filtered_snps %>% filter(outside_percentiles == TRUE), 
             aes(x = mean_effect, y = 50), color = "blue", size = 3) + # Highlight points outside percentiles
  geom_text(data = filtered_snps %>% filter(outside_percentiles == TRUE),
            aes(x = mean_effect, y = 60, label = rs), color = "blue", angle = 90, hjust = -0.2, size = 4) + # Add SNP names above dots
  theme_minimal() +
  labs(title = "Histogram of Marker Effects for Bio 17",
       x = "Mean Effect Size",
       y = "Frequency") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"))


#######################################################################################################
#SUPPLEMENTAL FIGURE 16
#######################################################################################################
################################
#GS with fixed effects for training n=31 showing model improves
################################
setwd("/Users/annamccormick/R/barley_collab/parental_selection/fixed_effects/")
library(rrBLUP)
# Load genotype data
genotype_data_original <- read.table("Barley_gd_numeric_rrblupformat_second.txt", header = TRUE)

# Set the first column as row names and remove it
row.names(genotype_data_original) <- genotype_data_original[,1]
genotype_data <- genotype_data_original[,-1]

# Convert the data frame to a matrix
genotype_matrix <- as.matrix(genotype_data)

# Load environmental data
Y <- read.csv('Supplemental_dataset_5_worldclim_training31_copy.csv', head = T) # full environmental dataset

# Define the SNPs of interest
snps_of_interest <- c("X11_10189", "X11_11328", "X11_11361", "X11_20085", 
                      "X12_30850", "X12_30880", "X11_10513", "X12_30668", 
                      "X12_30867", "X12_30869", "X12_30871", "X12_30872", 
                      "X12_30929", "X12_30930", "BK_12", "BK_13", "BK_14", 
                      "BK_15", "BK_16", "SCRI_RS_137464", "SCRI_RS_142282", 
                      "SCRI_RS_154621")

# Filter the genotype data to include only the SNPs of interest
filtered_genotype_data <- genotype_data_original[, c("taxa", snps_of_interest)]

# Remove the taxa column and convert to a matrix
fixed_genotype_matrix <- as.matrix(filtered_genotype_data[,-1])

# Define the training population based on the Core column
training_population <- Y[Y$Core == TRUE, ]
train_taxa <- training_population$Taxa

# Filter the genotype matrix to include only the training population
genotype_matrix_train <- genotype_matrix[train_taxa, ]

# Filter the phenotype data to include only the training population
Y_train <- training_population[, 6]  # trait bio1

# Filter the fixed effects SNPs data to include only the training population
fixed_genotype_matrix_train <- fixed_genotype_matrix[train_taxa, ]

# Ensure Y_train is numeric
Y_train <- as.numeric(Y_train)

# Baseline mixed model analysis
result_baseline <- mixed.solve(y = Y_train, 
                               Z = genotype_matrix_train, 
                               X = matrix(1, nrow = length(Y_train)), # Intercept-only model
                               K = NULL, 
                               method = "REML")

# Fixed Effects Analysis (With SNPs of Interest as Fixed Effects)
# Ensure there are no zero-variance SNPs
zero_variance_snps <- apply(fixed_genotype_matrix_train, 2, var) == 0
fixed_genotype_matrix_train_no_zero_var <- fixed_genotype_matrix_train[, !zero_variance_snps]

# Calculate the correlation matrix on the cleaned data
cor_matrix <- cor(fixed_genotype_matrix_train_no_zero_var)
library(caret)
# Identify highly correlated SNPs
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.9)

# Remove highly correlated SNPs
fixed_genotype_matrix_train_reduced <- fixed_genotype_matrix_train_no_zero_var[, -highly_correlated]

# Mixed Model with Fixed Effects
result_fixed <- mixed.solve(y = Y_train, 
                            Z = genotype_matrix_train, 
                            X = cbind(matrix(1, nrow = length(Y_train)), fixed_genotype_matrix_train_reduced), 
                            K = NULL, 
                            method = "REML")

# Compare Models
log_likelihood_baseline <- result_baseline$LL
log_likelihood_fixed <- result_fixed$LL
stat <- -2 * (log_likelihood_baseline - log_likelihood_fixed)
p_value <- pchisq(q = stat, df = ncol(fixed_genotype_matrix_train_reduced), lower.tail = FALSE)

# Print Results
cat("Log-likelihood of baseline model: ", log_likelihood_baseline, "\n")
cat("Log-likelihood of model with fixed effects: ", log_likelihood_fixed, "\n")
cat("Chi-square statistic: ", stat, "\n")
cat("p-value: ", p_value, "\n")

#######################################################################################################
#SUPPLEMENTAL FIGURE 16A
#######################################################################################################
#exploring cold tolerance snps only 
#########################################################
#Core n=31 rrBLUP with fixed effects
##########################################
library(rrBLUP)
library(hibayes)
library(dplyr)

#setwd("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_4/")
#source("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation/xval_kfold_functions_updated.R") # load all functions
setwd("~/R/barley_collab/parental_selection/cross_validation_5/")
source("~/R/barley_collab/parental_selection/cross_validation_5/xval_kfold_functions_updated_2.R") # load all functions


gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = T)
envdat <- read.csv('Supplemental_dataset_5worldclim_all.csv', head = T) # full environmental dataset
trainingset <- read.csv("Supplemental_dataset_5_worldclim_training31.csv", head = T)


gd2 <- gd1[gd1$taxa %in% envdat$Taxa,]
row.names(gd2) <- gd2$taxa # set taxa as rownames
g.in <- gd2[,-1]
g.in <- as.matrix(g.in)

### Load Environmental Data ###
row.names(envdat) <- envdat$Taxa # set gen_id as rownames
row.names(trainingset) <- trainingset$Taxa

y.trainset <- trainingset[,c(1,5:(ncol(trainingset)))] # select unique identifier and environmental data only
y.in <- envdat[,c(1,4:ncol(envdat))]

# Exclude "altitude" from y.in and y.trainset
y.in <- y.in[, !(colnames(y.in) %in% c("altitude", "Taxa"))]
y.trainset <- y.trainset[, !(colnames(y.trainset) %in% c("altitude", "Taxa"))]

print(colnames(y.in))  # Check the column names in y.in
print(colnames(y.trainset))  # Check the column names in y.trainset


snps_of_interest <- c("X11_10189", "X11_11328", "X11_11361", "X11_20085", "X12_30880")


# Subset the genotype matrix for the fixed SNPs
fixed.snps <- g.in[, snps_of_interest] # Ensure snps_of_interest are column names in g.in

# Ensure that rownames of fixed.snps are aligned with g.in
rownames(fixed.snps) <- rownames(g.in)

### RR-BLUP
y.trainset.rr <- trainingset[,c(5:(ncol(trainingset)))]
y.in.rr <- envdat[,5:ncol(envdat)]

y.trainset.mat <- as.matrix(y.trainset.rr)
y.in.mat <- as.matrix(y.in.rr)

# Run k-fold cross-validation using fixed.snps
xval_k10_rrblup <- k.xval(g.in = g.in, 
                          y.in = y.in.mat, 
                          y.trainset = y.trainset.mat, 
                          k.fold = 10, 
                          reps = 50, 
                          fixed.snps = fixed.snps) # Pass fixed.snps matrix here

saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10_fixed_cold31_1.RData")

#Core n=100 rrBLUP with fixed effects
##########################################
library(rrBLUP)
library(hibayes)
library(dplyr)

# Set working directory and load functions
setwd("~/R/barley_collab/parental_selection/cross_validation_5_100/")
source("~/R/barley_collab/parental_selection/cross_validation_5/xval_kfold_functions_updated_2.R") # load all functions

# Load genotype and environmental data
gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = T)
envdat <- read.csv('Supplemental_dataset_5worldclim_all.csv', head = T) # full environmental dataset
trainingset <- read.csv("Supplemental_dataset_5_worldclim_training100.csv", head = T)

# Filter the genotype data to match taxa in envdat
gd2 <- gd1[gd1$taxa %in% envdat$Taxa,]
row.names(gd2) <- gd2$taxa
g.in <- as.matrix(gd2[,-1])

# Load environmental data
row.names(envdat) <- envdat$Taxa
row.names(trainingset) <- trainingset$Taxa

# Subset y.trainset and y.in for environmental data only, excluding altitude and other unneeded columns
y.trainset <- trainingset[, !(colnames(trainingset) %in% c("altitude","Latitude","Longitude","Core" ,"Taxa", "bio2", "bio5", "bio7", "bio8", "bio9", "bio10", "bio12", "bio13", "bio15", "bio16", "bio18", "bio19"))]
y.in <- envdat[, !(colnames(envdat) %in% c("altitude", "Taxa","Latitude","Longitude"))]

print(colnames(y.in))      # Check the column names in y.in
print(colnames(y.trainset)) # Check the column names in y.trainset

# List of SNPs of interest
snps_of_interest <- c("X11_10189", "X11_11328", "X11_11361", "X11_20085", "X12_30880")

# Subset the genotype matrix for the fixed SNPs
fixed.snps <- g.in[, snps_of_interest] # Ensure snps_of_interest are column names in g.in
rownames(fixed.snps) <- rownames(g.in) # Ensure alignment of rownames

xval_k10_rrblup <- k.xval(g.in = g.in, 
                          y.in = as.matrix(y.in),      # Convert y.in to matrix
                          y.trainset = as.matrix(y.trainset),  # Convert y.trainset to matrix
                          k.fold = 10, 
                          reps = 50, 
                          fixed.snps = fixed.snps) # Pass fixed.snps matrix here

saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10_fixed_cold100_2.RData")

###############
#PLOTS
library(ggplot2)
#############################################################################
# Load cross-validation results for the first run (Core 31 - Original)

setwd("~/R/barley_collab/parental_selection/cross_validation/")

rrblup_kfold10_run1 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10_run1 <- rrblup_kfold10_run1$xval.result
rrblup_kfold10_run1$r.mean <- as.numeric(rrblup_kfold10_run1$r.mean)
rrblup_kfold10_run1$model <- "rrBLUP"
rrblup_kfold10_run1$xval <- "Ten-Fold"
rrblup_kfold10_run1$run <- "Core n=31 - Original"  # Add identifier for the first run

#############################################################################
# Load cross-validation results for the second run (Core 31 - Fixed SNPs)

setwd("~/R/barley_collab/parental_selection/cross_validation_5/")

rrblup_kfold10_run2 <- readRDS("xval_rrblup_kfold_10_fixed_cold31_1.RData")
rrblup_kfold10_run2 <- rrblup_kfold10_run2$xval.result
rrblup_kfold10_run2$r.mean <- as.numeric(rrblup_kfold10_run2$r.mean)
rrblup_kfold10_run2$model <- "rrBLUP"
rrblup_kfold10_run2$xval <- "Ten-Fold"
rrblup_kfold10_run2$run <- "Core n=31 - Fixed effect SNPs"  # Add identifier for the second run

#############################################################################
# Load cross-validation results for the third run (Core 100 - Original)

setwd("~/R/barley_collab/parental_selection/core_comparisons/100/")

rrblup_kfold10_run3 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10_run3 <- rrblup_kfold10_run3$xval.result
rrblup_kfold10_run3$r.mean <- as.numeric(rrblup_kfold10_run3$r.mean)
rrblup_kfold10_run3$model <- "rrBLUP"
rrblup_kfold10_run3$xval <- "Ten-Fold"
rrblup_kfold10_run3$run <- "Core n=100 - Original"  # Add identifier for the third run

#############################################################################
# Load cross-validation results for the fourth run (Core 100 - Fixed SNPs)

setwd("~/R/barley_collab/parental_selection/cross_validation_5_100/")

rrblup_kfold10_run4 <- readRDS("xval_rrblup_kfold_10_fixed_cold100_2.RData")
rrblup_kfold10_run4 <- rrblup_kfold10_run4$xval.result
rrblup_kfold10_run4$r.mean <- as.numeric(rrblup_kfold10_run4$r.mean)
rrblup_kfold10_run4$model <- "rrBLUP"
rrblup_kfold10_run4$xval <- "Ten-Fold"
rrblup_kfold10_run4$run <- "Core n=100 - Fixed effect SNPs"  # Add identifier for the fourth run

#############################################################################
# Combine the four datasets
all_runs <- rbind(rrblup_kfold10_run1, rrblup_kfold10_run2, rrblup_kfold10_run3, rrblup_kfold10_run4)

# Convert standard deviation values to numeric
all_runs$r.sd <- as.numeric(all_runs$r.sd)

# Ensure cross-validation type and run are factors with specified levels
all_runs$xval <- factor(all_runs$xval, levels = c("Ten-Fold"))
all_runs$run <- factor(all_runs$run, levels = c("Core n=31 - Original", "Core n=31 - Fixed effect SNPs", "Core n=100 - Original", "Core n=100 - Fixed effect SNPs"))

# Filter for specific traits of interest
all_bio <- all_runs[all_runs$trait %in% c('bio1', 'bio3', 'bio4', 'bio6', 'bio11', 'bio14', 'bio17'),]
all_bio$trait <- factor(all_bio$trait, levels = c('bio1', 'bio3', 'bio4', 'bio6', 'bio11', 'bio14', 'bio17'))

#################################################################################################################
# Plot - Comparing All Four Runs Side by Side
ggplot(all_bio, aes(y = r.mean, x = model, color = run)) +  # Color by 'run' to differentiate runs
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(color = "Run")  # Change the label for the legend

#######################################################################################################
#SUPPLEMENTAL FIGURE 16B
#######################################################################################################
#Flowering time related SNPs fixed ONLY
#########################################################
#Core n=31 rrBLUP with fixed effects
##########################################
library(rrBLUP)
library(hibayes)
library(dplyr)

#setwd("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_4/")
#source("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation/xval_kfold_functions_updated.R") # load all functions
setwd("~/R/barley_collab/parental_selection/cross_validation_5/")
source("~/R/barley_collab/parental_selection/cross_validation_5/xval_kfold_functions_updated_2.R") # load all functions


gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = T)
envdat <- read.csv('Supplemental_dataset_5worldclim_all.csv', head = T) # full environmental dataset
trainingset <- read.csv("Supplemental_dataset_5_worldclim_training31.csv", head = T)


gd2 <- gd1[gd1$taxa %in% envdat$Taxa,]
row.names(gd2) <- gd2$taxa # set taxa as rownames
g.in <- gd2[,-1]
g.in <- as.matrix(g.in)

### Load Environmental Data ###
row.names(envdat) <- envdat$Taxa # set gen_id as rownames
row.names(trainingset) <- trainingset$Taxa

y.trainset <- trainingset[,c(1,5:(ncol(trainingset)))] # select unique identifier and environmental data only
y.in <- envdat[,c(1,4:ncol(envdat))]

# Exclude "altitude" from y.in and y.trainset
y.in <- y.in[, !(colnames(y.in) %in% c("altitude", "Taxa"))]
y.trainset <- y.trainset[, !(colnames(y.trainset) %in% c("altitude", "Taxa"))]

print(colnames(y.in))  # Check the column names in y.in
print(colnames(y.trainset))  # Check the column names in y.trainset

# List of SNPs of interest (Fixed effects)
snps_of_interest <- c("X11_10513",
                      "X12_30867", "X12_30869", "X12_30929", "X12_30930", "BK_13",
                      "SCRI_RS_137464", "SCRI_RS_142282", "SCRI_RS_154621")


# Subset the genotype matrix for the fixed SNPs
fixed.snps <- g.in[, snps_of_interest] # Ensure snps_of_interest are column names in g.in

# Ensure that rownames of fixed.snps are aligned with g.in
rownames(fixed.snps) <- rownames(g.in)

### RR-BLUP
y.trainset.rr <- trainingset[,c(5:(ncol(trainingset)))]
y.in.rr <- envdat[,5:ncol(envdat)]

y.trainset.mat <- as.matrix(y.trainset.rr)
y.in.mat <- as.matrix(y.in.rr)

# Run k-fold cross-validation using fixed.snps
xval_k10_rrblup <- k.xval(g.in = g.in, 
                          y.in = y.in.mat, 
                          y.trainset = y.trainset.mat, 
                          k.fold = 10, 
                          reps = 50, 
                          fixed.snps = fixed.snps) # Pass fixed.snps matrix here

saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10_fixed_flowering31_1.RData")


#Core n=100 rrBLUP with fixed effects
##########################################
library(rrBLUP)
library(hibayes)
library(dplyr)

# Set working directory and load functions
setwd("~/R/barley_collab/parental_selection/cross_validation_5_100/")
source("~/R/barley_collab/parental_selection/cross_validation_5/xval_kfold_functions_updated_2.R") # load all functions

# Load genotype and environmental data
gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = T)
envdat <- read.csv('Supplemental_dataset_5worldclim_all.csv', head = T) # full environmental dataset
trainingset <- read.csv("Supplemental_dataset_5_worldclim_training100.csv", head = T)

# Filter the genotype data to match taxa in envdat
gd2 <- gd1[gd1$taxa %in% envdat$Taxa,]
row.names(gd2) <- gd2$taxa
g.in <- as.matrix(gd2[,-1])

# Load environmental data
row.names(envdat) <- envdat$Taxa
row.names(trainingset) <- trainingset$Taxa

# Subset y.trainset and y.in for environmental data only, excluding altitude and other unneeded columns
y.trainset <- trainingset[, !(colnames(trainingset) %in% c("altitude","Latitude","Longitude","Core" ,"Taxa", "bio2", "bio5", "bio7", "bio8", "bio9", "bio10", "bio12", "bio13", "bio15", "bio16", "bio18", "bio19"))]
y.in <- envdat[, !(colnames(envdat) %in% c("altitude", "Taxa","Latitude","Longitude"))]

print(colnames(y.in))      # Check the column names in y.in
print(colnames(y.trainset)) # Check the column names in y.trainset

# List of SNPs of interest
snps_of_interest <- c("X11_10513",
                      "X12_30867", "X12_30869", "X12_30929", "X12_30930", "BK_13",
                      "SCRI_RS_137464", "SCRI_RS_142282", "SCRI_RS_154621")

# Subset the genotype matrix for the fixed SNPs
fixed.snps <- g.in[, snps_of_interest] # Ensure snps_of_interest are column names in g.in
rownames(fixed.snps) <- rownames(g.in) # Ensure alignment of rownames

xval_k10_rrblup <- k.xval(g.in = g.in, 
                          y.in = as.matrix(y.in),      # Convert y.in to matrix
                          y.trainset = as.matrix(y.trainset),  # Convert y.trainset to matrix
                          k.fold = 10, 
                          reps = 50, 
                          fixed.snps = fixed.snps) # Pass fixed.snps matrix here

saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10_fixed_flowering100_2.RData")

###############
#PLOTS
library(ggplot2)

#############################################################################
# Load cross-validation results for the first run (Core 31 - Original)

setwd("~/R/barley_collab/parental_selection/cross_validation/")

rrblup_kfold10_run1 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10_run1 <- rrblup_kfold10_run1$xval.result
rrblup_kfold10_run1$r.mean <- as.numeric(rrblup_kfold10_run1$r.mean)
rrblup_kfold10_run1$model <- "rrBLUP"
rrblup_kfold10_run1$xval <- "Ten-Fold"
rrblup_kfold10_run1$run <- "Core n=31 - Original"  # Add identifier for the first run

#############################################################################
# Load cross-validation results for the second run (Core 31 - Fixed SNPs)

setwd("~/R/barley_collab/parental_selection/cross_validation_5/")

rrblup_kfold10_run2 <- readRDS("xval_rrblup_kfold_10_fixed_flowering31_1.RData")
rrblup_kfold10_run2 <- rrblup_kfold10_run2$xval.result
rrblup_kfold10_run2$r.mean <- as.numeric(rrblup_kfold10_run2$r.mean)
rrblup_kfold10_run2$model <- "rrBLUP"
rrblup_kfold10_run2$xval <- "Ten-Fold"
rrblup_kfold10_run2$run <- "Core n=31 - Fixed effect SNPs"  # Add identifier for the second run

#############################################################################
# Load cross-validation results for the third run (Core 100 - Original)

setwd("~/R/barley_collab/parental_selection/core_comparisons/100/")

rrblup_kfold10_run3 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10_run3 <- rrblup_kfold10_run3$xval.result
rrblup_kfold10_run3$r.mean <- as.numeric(rrblup_kfold10_run3$r.mean)
rrblup_kfold10_run3$model <- "rrBLUP"
rrblup_kfold10_run3$xval <- "Ten-Fold"
rrblup_kfold10_run3$run <- "Core n=100 - Original"  # Add identifier for the third run

#############################################################################
# Load cross-validation results for the fourth run (Core 100 - Fixed SNPs)

setwd("~/R/barley_collab/parental_selection/cross_validation_5_100/")

rrblup_kfold10_run4 <- readRDS("xval_rrblup_kfold_10_fixed_flowering100_2.RData")
rrblup_kfold10_run4 <- rrblup_kfold10_run4$xval.result
rrblup_kfold10_run4$r.mean <- as.numeric(rrblup_kfold10_run4$r.mean)
rrblup_kfold10_run4$model <- "rrBLUP"
rrblup_kfold10_run4$xval <- "Ten-Fold"
rrblup_kfold10_run4$run <- "Core n=100 - Fixed effect SNPs"  # Add identifier for the fourth run

#############################################################################
# Combine the four datasets
all_runs <- rbind(rrblup_kfold10_run1, rrblup_kfold10_run2, rrblup_kfold10_run3, rrblup_kfold10_run4)

# Convert standard deviation values to numeric
all_runs$r.sd <- as.numeric(all_runs$r.sd)

# Ensure cross-validation type and run are factors with specified levels
all_runs$xval <- factor(all_runs$xval, levels = c("Ten-Fold"))
all_runs$run <- factor(all_runs$run, levels = c("Core n=31 - Original", "Core n=31 - Fixed effect SNPs", "Core n=100 - Original", "Core n=100 - Fixed effect SNPs"))

# Filter for specific traits of interest
all_bio <- all_runs[all_runs$trait %in% c('bio1', 'bio3', 'bio4', 'bio6', 'bio11', 'bio14', 'bio17'),]
all_bio$trait <- factor(all_bio$trait, levels = c('bio1', 'bio3', 'bio4', 'bio6', 'bio11', 'bio14', 'bio17'))

#################################################################################################################
# Plot - Comparing All Four Runs Side by Side
ggplot(all_bio, aes(y = r.mean, x = model, color = run)) +  # Color by 'run' to differentiate runs
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(color = "Run")  # Change the label for the legend

#######################################################################################################
#SUPPLEMENTAL FIGURE 17
#######################################################################################################
#FIXED EFFECTS GENOMIC SELECTION
##########################
# Genomic selection (n=31 core) with SNPs as fixed effects
##########################
library(rrBLUP)
library(dplyr)

# Set working directory
setwd("~/R/barley_collab/parental_selection/cross_validation/")

# Load environmental and genotype data
envdat <- read.csv('Supplemental_dataset_5_worldclim_training31_copy.csv', head = TRUE)
gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = TRUE)

# Set row names for genotype data
row.names(gd1) <- gd1$taxa
gd3 <- gd1[,-1]  # Remove 'taxa' column for analysis
g.in <- as.matrix(gd3)

# RR-BLUP setup
row.names(envdat) <- envdat$Taxa
y.in.rr <- envdat[,5:ncol(envdat)]  # Data starts at column 5 for Core 31

# Subset to core lines
minicore_entries <- which(envdat$Core == TRUE)
y.in.rr <- y.in.rr[minicore_entries,]
y.in.rr <- y.in.rr[,-1]  # Remove 'Core' column
y.in.mat <- as.matrix(y.in.rr)

# Set up training and genotype data
train <- row.names(y.in.mat)  # Names of training lines
g.train <- g.in[train,]  # Set training genotypes

# List of SNPs of interest (Fixed effects)
snps_of_interest <- c("X11_10189", "X11_11328", "X11_11361", "X11_20085", "X12_30880", "X11_10513",
                      "X12_30867", "X12_30869", "X12_30929", "X12_30930", "BK_13",
                      "SCRI_RS_137464", "SCRI_RS_142282", "SCRI_RS_154621")

# Subset the genotype matrix for the fixed SNPs
fixed.snps <- g.in[, snps_of_interest]  # Ensure snps_of_interest are column names in g.in
rownames(fixed.snps) <- rownames(g.in)  # Ensure alignment of row names

# List of traits to analyze
traits <- colnames(y.in.mat)

# Prediction setup
pred <- setdiff(row.names(g.in), train)  # Names of lines not in training set
g.pred <- g.in[pred,]

# Initialize objects for storing results
marker.list <- list()
gebv_df <- data.frame(matrix(nrow = nrow(g.in), ncol = length(traits)))  # Initialize with all lines (train + pred)
colnames(gebv_df) <- traits

# RR-BLUP loop for each trait with fixed effects
for(t in 1:length(traits)) {
  trait <- traits[t]
  y.train <- as.matrix(y.in.mat[train, trait])  # Training set for the trait
  
  # Run RR-BLUP model with fixed SNPs
  solve.out <- mixed.solve(y = y.train, Z = g.train, X = fixed.snps[train,], SE = FALSE, return.Hinv = FALSE)
  u.hat <- solve.out$u
  
  # Calculate GEBVs for both prediction and training set
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  
  # Use match to ensure correct row alignment
  pred_rows <- match(row.names(g.pred), row.names(g.in))  # Indices for pred rows
  train_rows <- match(row.names(g.train), row.names(g.in))  # Indices for train rows
  
  # Store the results in the combined dataframe
  gebv_df[pred_rows, t] <- GEBV  # Predictions for test lines
  gebv_df[train_rows, t] <- GEBV_train  # Predictions for training lines
}

# Set row names for the result to match the line names from `g.in`
row.names(gebv_df) <- row.names(g.in)

# Save results
write.csv(gebv_df, 'rrblup_GEBV_data_barley_31_with_fixed_snps.csv')


##########################
# Genomic selection (n=100 core) with SNPs as fixed effects
##########################
library(rrBLUP)
library(dplyr)

# Set the working directory for the 100-core dataset
setwd("~/R/barley_collab/parental_selection/new_core_100/")

# Load the full environmental dataset and genotype data
envdat <- read.csv('Supplemental_dataset_5_worldclim_training100.csv', head = TRUE)  # Full environmental dataset for 100-core
gd1 <- read.table("Barley_gd_numeric_rrblupformat_second.txt", head = TRUE)  # Genotype data

# Set row names for genotype data
row.names(gd1) <- gd1$taxa
gd3 <- gd1[,-1]  # Remove 'taxa' column for analysis
g.in <- as.matrix(gd3)

# RR-BLUP setup
row.names(envdat) <- envdat$Taxa
y.in.rr <- envdat[,5:ncol(envdat)]  # Data starts at column 5 for Core 100

# Subset to core lines
minicore_entries <- which(envdat$Core == TRUE)
y.in.rr <- y.in.rr[minicore_entries,]
y.in.rr <- y.in.rr[,-1]  # Remove 'Core' column
y.in.mat <- as.matrix(y.in.rr)

# Set up training and genotype data
train <- row.names(y.in.mat)  # Names of training lines
g.train <- g.in[train,]  # Set training genotypes

# List of SNPs of interest (Fixed effects)
snps_of_interest <- c("X11_10189", "X11_11328", "X11_11361", "X11_20085", "X12_30880", "X11_10513",
                      "X12_30867", "X12_30869", "X12_30929", "X12_30930", "BK_13",
                      "SCRI_RS_137464", "SCRI_RS_142282", "SCRI_RS_154621")

# Subset the genotype matrix for the fixed SNPs
fixed.snps <- g.in[, snps_of_interest]  # Ensure snps_of_interest are column names in g.in
rownames(fixed.snps) <- rownames(g.in)  # Ensure alignment of row names

# List of traits to analyze
traits <- colnames(y.in.mat)

# Prediction setup
pred <- setdiff(row.names(g.in), train)  # Names of lines not in training set
g.pred <- g.in[pred,]

# Initialize objects for storing results
marker.list <- list()
gebv_df <- data.frame(matrix(nrow = nrow(g.in), ncol = length(traits)))  # Initialize with all lines (train + pred)
colnames(gebv_df) <- traits

# RR-BLUP loop for each trait with fixed effects
for(t in 1:length(traits)) {
  trait <- traits[t]
  y.train <- as.matrix(y.in.mat[train, trait])  # Training set for the trait
  
  # Run RR-BLUP model with fixed SNPs
  solve.out <- mixed.solve(y = y.train, Z = g.train, X = fixed.snps[train,], SE = FALSE, return.Hinv = FALSE)
  u.hat <- solve.out$u
  
  # Calculate GEBVs for both prediction and training set
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  
  # Use match to ensure correct row alignment
  pred_rows <- match(row.names(g.pred), row.names(g.in))  # Indices for pred rows
  train_rows <- match(row.names(g.train), row.names(g.in))  # Indices for train rows
  
  # Store the results in the combined dataframe
  gebv_df[pred_rows, t] <- GEBV  # Predictions for test lines
  gebv_df[train_rows, t] <- GEBV_train  # Predictions for training lines
}

# Set row names for the result to match the line names from `g.in`
row.names(gebv_df) <- row.names(g.in)

# Save results
write.csv(gebv_df, 'rrblup_GEBV_data_barley_100_with_fixed_snps.csv')

######
#Plots for rank changes
######
#BIO1
library(ggplot2)
library(ggalluvial)
library(dplyr)
setwd("~/R/barley_collab/parental_selection/ggsanky/bio1")
df <- read.csv("all_ranks_bio1.csv")

df_long <- df %>%
  pivot_longer(cols = starts_with("core"), names_to = "Stage", values_to = "Rank")

df_long$Stage <- factor(df_long$Stage, levels = c("core_31_fixed_snp", "core_31", "core_100", "core_100_fixed_snp"))

ggplot(df_long, aes(x = Stage, stratum = Rank, alluvium = lines, fill = as.numeric(Rank))) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(fill = as.numeric(Rank)), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_viridis_c(option = "viridis") +  # This will create the rainbow color gradient
  labs(title = "Rank Changes for Bio 1", x = "Stage", y = "Rank") +
  theme_minimal() +
  theme(legend.position = "none")

##############
#BIO3
setwd("~/R/barley_collab/parental_selection/ggsanky/bio3")
df <- read.csv("all_ranks_bio3.csv")

df_long <- df %>%
  pivot_longer(cols = starts_with("core"), names_to = "Stage", values_to = "Rank")

df_long$Stage <- factor(df_long$Stage, levels = c("core_31_fixed_snp", "core_31", "core_100", "core_100_fixed_snp"))

ggplot(df_long, aes(x = Stage, stratum = Rank, alluvium = lines, fill = as.numeric(Rank))) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(fill = as.numeric(Rank)), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_viridis_c(option = "viridis") +  # This will create the rainbow color gradient
  labs(title = "Rank Changes for Bio 3", x = "Stage", y = "Rank") +
  theme_minimal() +
  theme(legend.position = "none")

##############
#BIO4
setwd("~/R/barley_collab/parental_selection/ggsanky/bio4")
df <- read.csv("all_ranks_bio4.csv")

df_long <- df %>%
  pivot_longer(cols = starts_with("core"), names_to = "Stage", values_to = "Rank")

df_long$Stage <- factor(df_long$Stage, levels = c("core_31_fixed_snp", "core_31", "core_100", "core_100_fixed_snp"))

ggplot(df_long, aes(x = Stage, stratum = Rank, alluvium = lines, fill = as.numeric(Rank))) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(fill = as.numeric(Rank)), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_viridis_c(option = "viridis") +  # This will create the rainbow color gradient
  labs(title = "Rank Changes for Bio 4", x = "Stage", y = "Rank") +
  theme_minimal() +
  theme(legend.position = "none")


##############
#BIO6
setwd("~/R/barley_collab/parental_selection/ggsanky/bio6")
df <- read.csv("all_ranks_bio6.csv")

df_long <- df %>%
  pivot_longer(cols = starts_with("core"), names_to = "Stage", values_to = "Rank")

df_long$Stage <- factor(df_long$Stage, levels = c("core_31_fixed_snp", "core_31", "core_100", "core_100_fixed_snp"))

ggplot(df_long, aes(x = Stage, stratum = Rank, alluvium = lines, fill = as.numeric(Rank))) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(fill = as.numeric(Rank)), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_viridis_c(option = "viridis") +  # This will create the rainbow color gradient
  labs(title = "Rank Changes for Bio 6", x = "Stage", y = "Rank") +
  theme_minimal() +
  theme(legend.position = "none")

##############
#BIO11
setwd("~/R/barley_collab/parental_selection/ggsanky/bio11")
df <- read.csv("all_ranks_bio11.csv")

df_long <- df %>%
  pivot_longer(cols = starts_with("core"), names_to = "Stage", values_to = "Rank")

df_long$Stage <- factor(df_long$Stage, levels = c("core_31_fixed_snp", "core_31", "core_100", "core_100_fixed_snp"))

ggplot(df_long, aes(x = Stage, stratum = Rank, alluvium = lines, fill = as.numeric(Rank))) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(fill = as.numeric(Rank)), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_viridis_c(option = "viridis") +  # This will create the rainbow color gradient
  labs(title = "Rank Changes for Bio 11", x = "Stage", y = "Rank") +
  theme_minimal() +
  theme(legend.position = "none")


##############
#BIO14
setwd("~/R/barley_collab/parental_selection/ggsanky/bio14")
df <- read.csv("all_ranks_bio14.csv")

df_long <- df %>%
  pivot_longer(cols = starts_with("core"), names_to = "Stage", values_to = "Rank")

df_long$Stage <- factor(df_long$Stage, levels = c("core_31_fixed_snp", "core_31", "core_100", "core_100_fixed_snp"))

ggplot(df_long, aes(x = Stage, stratum = Rank, alluvium = lines, fill = as.numeric(Rank))) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(fill = as.numeric(Rank)), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_viridis_c(option = "viridis") +  # This will create the rainbow color gradient
  labs(title = "Rank Changes for Bio 14", x = "Stage", y = "Rank") +
  theme_minimal() +
  theme(legend.position = "none")

##############
#BIO17
setwd("~/R/barley_collab/parental_selection/ggsanky/bio17")
df <- read.csv("all_ranks_bio17.csv")

df_long <- df %>%
  pivot_longer(cols = starts_with("core"), names_to = "Stage", values_to = "Rank")

df_long$Stage <- factor(df_long$Stage, levels = c("core_31_fixed_snp", "core_31", "core_100", "core_100_fixed_snp"))

ggplot(df_long, aes(x = Stage, stratum = Rank, alluvium = lines, fill = as.numeric(Rank))) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(fill = as.numeric(Rank)), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_viridis_c(option = "viridis") +  # This will create the rainbow color gradient
  labs(title = "Rank Changes for Bio 17", x = "Stage", y = "Rank") +
  theme_minimal() +
  theme(legend.position = "none")

#################################################################################################################
# FIN
##################################################################################################################
