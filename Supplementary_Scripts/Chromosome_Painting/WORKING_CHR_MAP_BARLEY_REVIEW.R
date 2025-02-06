# Load necessary libraries
library(dplyr)
library(circlize)

# Set working directory
#setwd("/home/ahmccorm/kantar_koastore/anna/Barley_collab/barley_parental/chromosome_painting_redo/Group1/")
setwd("~annamccormick/R/barley_collab/parental_selection/Barley_chromosome_paint/")

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


# Define axis labels and positions for every 100 Mb
labs <- paste0(seq(from = 0, to = 500, by = 100), "Mb")  # Labels for 100 Mb intervals
ma <- c(0, 165, 171, 174, 265, 391)  # Approximate positions in the dataset closest to each interval


# Sort bio1.chrom1.plot by position
bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  select(position, all_of(Group2)) %>%
  arrange(position)  # Ensure SNPs are ordered by increasing position


dim(bio1.chrom1.plot)
head(bio1.chrom1.plot$position)
tail(bio1.chrom1.plot$position)


# Define the color function
col_fun1 = colorRamp2(c(-0.12, -0.06, 0, 0.06, 0.12), 
                      c("#701130", "#e06e85", "white", "#5d94cb", "#163670"))

# Save plot to PDF
pdf("circos_plot_chrom1_group2.pdf", width = 10, height = 10)  # Open PDF device

# Clear and set circos parameters
circos.clear()
circos.par(start.degree = 65, gap.after = 45)

# Draw the heatmap
circos.heatmap(bio1.chrom1.plot[, 2:90], 
               col = col_fun1, 
               track.height = 0.8, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1)

# Add position labels outside the plot
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(ma, rep(CELL_META$ylim[2], length(ma)), labs, 
              facing = "clockwise", 
              niceFacing = TRUE, 
              cex = 0.75, 
              adj = c(0.5, 0))
}, bg.border = NA)  # No border for the labels track

dev.off()  # Close the PDF device



############################
#Group 1 - Chrom 1
############################
library(dplyr)
library(circlize)

# Set working directory
#setwd("/home/ahmccorm/kantar_koastore/anna/Barley_collab/barley_parental/chromosome_painting_redo/Group1/")
setwd("~annamccormick/R/barley_collab/parental_selection/Barley_chromosome_paint/")

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


# Define axis labels and positions for every 100 Mb
labs <- paste0(seq(from = 0, to = 500, by = 100), "Mb")  # Labels for 100 Mb intervals
ma <- c(0, 165, 171, 174, 265, 391)  # Approximate positions in the dataset closest to each interval


# Sort bio1.chrom1.plot by position
bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  select(position, all_of(Group1)) %>%
  arrange(position)  # Ensure SNPs are ordered by increasing position


dim(bio1.chrom1.plot)
head(bio1.chrom1.plot$position)
tail(bio1.chrom1.plot$position)


# Define the color function
col_fun1 = colorRamp2(c(-0.12, -0.06, 0, 0.06, 0.12), 
                      c("#701130", "#e06e85", "white", "#5d94cb", "#163670"))

# Save plot to PDF
pdf("circos_plot_chrom1_group1a.pdf", width = 10, height = 10)  # Open PDF device

# Clear and set circos parameters
circos.clear()
circos.par(start.degree = 65, gap.after = 45)

# Draw the heatmap
circos.heatmap(bio1.chrom1.plot[, 2:90], 
               col = col_fun1, 
               track.height = 0.8, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1)


# Add position labels outside the plot with an increased offset
# Add position labels with staggered placement
# Dynamically adjust label positions to avoid overlap
# Add position labels with uniform height adjustment
circos.track(track.index = 1, panel.fun = function(x, y) {
  offset <- 5  # Base offset for labels
  stagger_increment <- 1  # Increment to stagger labels
  base_y <- CELL_META$ylim[2] + offset  # Base position
  max_stagger <- base_y + 0.5 * stagger_increment  # Maximum allowed height
  
  # Adjust y positions dynamically based on closeness
  staggered_y <- rep(base_y, length(ma))  # Start with base position
  for (i in seq_along(ma)) {
    if (i > 1 && abs(ma[i] - ma[i - 1]) < 80) {  # Adjust closeness threshold
      staggered_y[i] <- min(staggered_y[i - 1] + stagger_increment, max_stagger)
    }
  }
  
  # Add staggered labels with normalized height
  circos.text(ma, staggered_y, labs,
              facing = "clockwise",
              niceFacing = TRUE,
              cex = 0.6,
              adj = c(0.5, 0))
}, bg.border = NA)  # No border for the labels track

dev.off()  # Close the PDF device

############################
#Group 2 - Chrom 1
############################
library(dplyr)
library(circlize)

# Set working directory
#setwd("/home/ahmccorm/kantar_koastore/anna/Barley_collab/barley_parental/chromosome_painting_redo/Group1/")
setwd("~annamccormick/R/barley_collab/parental_selection/Barley_chromosome_paint/")

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


# Define axis labels and positions for every 100 Mb
labs <- paste0(seq(from = 0, to = 500, by = 100), "Mb")  # Labels for 100 Mb intervals
ma <- c(0, 165, 171, 174, 265, 391)  # Approximate positions in the dataset closest to each interval


# Sort bio1.chrom1.plot by position
bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  select(position, all_of(Group2)) %>%
  arrange(position)  # Ensure SNPs are ordered by increasing position


dim(bio1.chrom1.plot)
head(bio1.chrom1.plot$position)
tail(bio1.chrom1.plot$position)


# Define the color function
col_fun1 = colorRamp2(c(-0.12, -0.06, 0, 0.06, 0.12), 
                      c("#701130", "#e06e85", "white", "#5d94cb", "#163670"))

# Save plot to PDF
pdf("circos_plot_chrom1_group2a.pdf", width = 10, height = 10)  # Open PDF device

# Clear and set circos parameters
circos.clear()
circos.par(start.degree = 65, gap.after = 45)

# Draw the heatmap
circos.heatmap(bio1.chrom1.plot[, 2:90], 
               col = col_fun1, 
               track.height = 0.8, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1)


# Add position labels outside the plot with an increased offset
# Add position labels with staggered placement
# Dynamically adjust label positions to avoid overlap
# Add position labels with uniform height adjustment
circos.track(track.index = 1, panel.fun = function(x, y) {
  offset <- 5  # Base offset for labels
  stagger_increment <- 1  # Increment to stagger labels
  base_y <- CELL_META$ylim[2] + offset  # Base position
  max_stagger <- base_y + 0.5 * stagger_increment  # Maximum allowed height
  
  # Adjust y positions dynamically based on closeness
  staggered_y <- rep(base_y, length(ma))  # Start with base position
  for (i in seq_along(ma)) {
    if (i > 1 && abs(ma[i] - ma[i - 1]) < 80) {  # Adjust closeness threshold
      staggered_y[i] <- min(staggered_y[i - 1] + stagger_increment, max_stagger)
    }
  }
  
  # Add staggered labels with normalized height
  circos.text(ma, staggered_y, labs,
              facing = "clockwise",
              niceFacing = TRUE,
              cex = 0.6,
              adj = c(0.5, 0))
}, bg.border = NA)  # No border for the labels track

dev.off()  # Close the PDF device



############################
#Group 3 - Chrom 1
############################
library(dplyr)
library(circlize)

# Set working directory
#setwd("/home/ahmccorm/kantar_koastore/anna/Barley_collab/barley_parental/chromosome_painting_redo/Group1/")
setwd("~annamccormick/R/barley_collab/parental_selection/Barley_chromosome_paint/")

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


# Define axis labels and positions for every 100 Mb
labs <- paste0(seq(from = 0, to = 500, by = 100), "Mb")  # Labels for 100 Mb intervals
ma <- c(0, 165, 171, 174, 265, 391)  # Approximate positions in the dataset closest to each interval


# Sort bio1.chrom1.plot by position
bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  select(position, all_of(Group3)) %>%
  arrange(position)  # Ensure SNPs are ordered by increasing position


dim(bio1.chrom1.plot)
head(bio1.chrom1.plot$position)
tail(bio1.chrom1.plot$position)


# Define the color function
col_fun1 = colorRamp2(c(-0.12, -0.06, 0, 0.06, 0.12), 
                      c("#701130", "#e06e85", "white", "#5d94cb", "#163670"))

# Save plot to PDF
pdf("circos_plot_chrom1_group3a.pdf", width = 10, height = 10)  # Open PDF device

# Clear and set circos parameters
circos.clear()
circos.par(start.degree = 65, gap.after = 45)

# Draw the heatmap
circos.heatmap(bio1.chrom1.plot[, 2:90], 
               col = col_fun1, 
               track.height = 0.8, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1)


# Add position labels outside the plot with an increased offset
# Add position labels with staggered placement
# Dynamically adjust label positions to avoid overlap
# Add position labels with uniform height adjustment
circos.track(track.index = 1, panel.fun = function(x, y) {
  offset <- 5  # Base offset for labels
  stagger_increment <- 1  # Increment to stagger labels
  base_y <- CELL_META$ylim[2] + offset  # Base position
  max_stagger <- base_y + 0.5 * stagger_increment  # Maximum allowed height
  
  # Adjust y positions dynamically based on closeness
  staggered_y <- rep(base_y, length(ma))  # Start with base position
  for (i in seq_along(ma)) {
    if (i > 1 && abs(ma[i] - ma[i - 1]) < 80) {  # Adjust closeness threshold
      staggered_y[i] <- min(staggered_y[i - 1] + stagger_increment, max_stagger)
    }
  }
  
  # Add staggered labels with normalized height
  circos.text(ma, staggered_y, labs,
              facing = "clockwise",
              niceFacing = TRUE,
              cex = 0.6,
              adj = c(0.5, 0))
}, bg.border = NA)  # No border for the labels track

dev.off()  # Close the PDF device



############################
#Group 4 - Chrom 1
############################
library(dplyr)
library(circlize)

# Set working directory
#setwd("/home/ahmccorm/kantar_koastore/anna/Barley_collab/barley_parental/chromosome_painting_redo/Group1/")
setwd("~annamccormick/R/barley_collab/parental_selection/Barley_chromosome_paint/")

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


# Define axis labels and positions for every 100 Mb
labs <- paste0(seq(from = 0, to = 500, by = 100), "Mb")  # Labels for 100 Mb intervals
ma <- c(0, 165, 171, 174, 265, 391)  # Approximate positions in the dataset closest to each interval


# Sort bio1.chrom1.plot by position
bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  select(position, all_of(Group4)) %>%
  arrange(position)  # Ensure SNPs are ordered by increasing position


dim(bio1.chrom1.plot)
head(bio1.chrom1.plot$position)
tail(bio1.chrom1.plot$position)


# Define the color function
col_fun1 = colorRamp2(c(-0.12, -0.06, 0, 0.06, 0.12), 
                      c("#701130", "#e06e85", "white", "#5d94cb", "#163670"))

# Save plot to PDF
pdf("circos_plot_chrom1_group4a.pdf", width = 10, height = 10)  # Open PDF device

# Clear and set circos parameters
circos.clear()
circos.par(start.degree = 65, gap.after = 45)

# Draw the heatmap
circos.heatmap(bio1.chrom1.plot[, 2:90], 
               col = col_fun1, 
               track.height = 0.8, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1)


# Add position labels outside the plot with an increased offset
# Add position labels with staggered placement
# Dynamically adjust label positions to avoid overlap
# Add position labels with uniform height adjustment
circos.track(track.index = 1, panel.fun = function(x, y) {
  offset <- 5  # Base offset for labels
  stagger_increment <- 1  # Increment to stagger labels
  base_y <- CELL_META$ylim[2] + offset  # Base position
  max_stagger <- base_y + 0.5 * stagger_increment  # Maximum allowed height
  
  # Adjust y positions dynamically based on closeness
  staggered_y <- rep(base_y, length(ma))  # Start with base position
  for (i in seq_along(ma)) {
    if (i > 1 && abs(ma[i] - ma[i - 1]) < 80) {  # Adjust closeness threshold
      staggered_y[i] <- min(staggered_y[i - 1] + stagger_increment, max_stagger)
    }
  }
  
  # Add staggered labels with normalized height
  circos.text(ma, staggered_y, labs,
              facing = "clockwise",
              niceFacing = TRUE,
              cex = 0.6,
              adj = c(0.5, 0))
}, bg.border = NA)  # No border for the labels track

dev.off()  # Close the PDF device



############################
#Group 5 - Chrom 1
############################
library(dplyr)
library(circlize)

# Set working directory
#setwd("/home/ahmccorm/kantar_koastore/anna/Barley_collab/barley_parental/chromosome_painting_redo/Group1/")
setwd("~annamccormick/R/barley_collab/parental_selection/Barley_chromosome_paint/")

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


# Define axis labels and positions for every 100 Mb
labs <- paste0(seq(from = 0, to = 500, by = 100), "Mb")  # Labels for 100 Mb intervals
ma <- c(0, 165, 171, 174, 265, 391)  # Approximate positions in the dataset closest to each interval


# Sort bio1.chrom1.plot by position
bio1.chrom1.plot <- bio1_Marker_Effects_chrom_1H %>%
  select(position, all_of(Group5)) %>%
  arrange(position)  # Ensure SNPs are ordered by increasing position


dim(bio1.chrom1.plot)
head(bio1.chrom1.plot$position)
tail(bio1.chrom1.plot$position)


# Define the color function
col_fun1 = colorRamp2(c(-0.12, -0.06, 0, 0.06, 0.12), 
                      c("#701130", "#e06e85", "white", "#5d94cb", "#163670"))

# Save plot to PDF
pdf("circos_plot_chrom1_group5a.pdf", width = 10, height = 10)  # Open PDF device

# Clear and set circos parameters
circos.clear()
circos.par(start.degree = 65, gap.after = 45)

# Draw the heatmap
circos.heatmap(bio1.chrom1.plot[, 2:90], 
               col = col_fun1, 
               track.height = 0.8, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1)


# Add position labels outside the plot with an increased offset
# Add position labels with staggered placement
# Dynamically adjust label positions to avoid overlap
# Add position labels with uniform height adjustment
circos.track(track.index = 1, panel.fun = function(x, y) {
  offset <- 5  # Base offset for labels
  stagger_increment <- 1  # Increment to stagger labels
  base_y <- CELL_META$ylim[2] + offset  # Base position
  max_stagger <- base_y + 0.5 * stagger_increment  # Maximum allowed height
  
  # Adjust y positions dynamically based on closeness
  staggered_y <- rep(base_y, length(ma))  # Start with base position
  for (i in seq_along(ma)) {
    if (i > 1 && abs(ma[i] - ma[i - 1]) < 80) {  # Adjust closeness threshold
      staggered_y[i] <- min(staggered_y[i - 1] + stagger_increment, max_stagger)
    }
  }
  
  # Add staggered labels with normalized height
  circos.text(ma, staggered_y, labs,
              facing = "clockwise",
              niceFacing = TRUE,
              cex = 0.6,
              adj = c(0.5, 0))
}, bg.border = NA)  # No border for the labels track

dev.off()  # Close the PDF device



