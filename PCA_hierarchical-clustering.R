# Load necessary libraries
library(DESeq2)   # For differential expression analysis
library(dplyr)    # For data manipulation
library(ggplot2)  # For visualization

# Load the raw count matrix
load("counts.RData")
data.counts <- as.data.frame(data.counts)  # Convert to a data frame if necessary

# Filter out genes that have zero counts across all samples
data.counts <- data.counts[which(rowSums(data.counts) > 0),]

# Define experimental conditions for each sample
condition <- factor(c("NC", "NC", "NC",         # Negative control (NC)
                      "PRE_SAT", "PRE_SAT", "PRE_SAT", # Pre-evagination -TA
                      "PRE_AT", "PRE_AT", "PRE_AT",   # Pre-evagination +TA
                      "EV_SAT", "EV_SAT", "EV_SAT",   # Evagination -TA
                      "EV_AT", "EV_AT", "EV_AT",     # Evagination +TA
                      "POST_SAT", "POST_SAT", "POST_SAT", # Post-evagination -TA
                      "POST_AT", "POST_AT", "POST_AT"))  # Post-evagination +TA

# Create a metadata table linking samples to their conditions
coldata <- data.frame(row.names = colnames(data.counts), condition)

# Create a DESeqDataSet object for differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = data.counts,
                              colData = coldata,
                              design = ~condition)  # Model design based on condition

# Perform differential expression analysis
dds <- DESeq(dds)

# Save the DESeq2 object for future use
save(dds, file='dds.RData')

# Apply variance stabilizing transformation (VST) for normalized expression values
vsdata <- vst(dds, blind = FALSE)

# Generate a Principal Component Analysis (PCA) plot to visualize sample clustering
svg("plotPCA.svg", width = 8, height = 6) # open SVG device for plotPCA
plotPCA(vsdata, intgroup = "condition") # create PCA plot
dev.off() # close SVG device

# Extract transformed expression data
transformed_data <- assay(vsdata)

# Compute a sample-to-sample Euclidean distance matrix
dist_matrix <- dist(t(transformed_data), 
                    method = "euclidean", 
                    diag = TRUE,   # Include diagonal values
                    upper = TRUE)  # Display the upper triangle

# Convert the distance matrix to a square matrix format
dist_matrix_sq <- as.matrix(dist_matrix, header = TRUE, rownames = 1)

# Load pheatmap library for heatmap visualization
library(pheatmap)

# Save the heatmap as an SVG file
svg("SampleDist_heatmap.svg", width = 10, height = 8)

# Generate a hierarchical clustering heatmap of sample distances
pheatmap(dist_matrix_sq, 
         border_color = NA,  # Remove grid borders
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         color = colorRampPalette(rev(c("#deebf7", "#9ecae1", "#3182bd")))(100), # Blue color gradient
         angle_col = 45,  # Rotate column labels for readability
         fontsize_col = 7, 
         fontsize_row = 7, 
         legend_title_fontsize = 8, 
         legend_labes_fontsize = 6)  # Adjust font sizes for clarity

# Close the SVG device to finalize the plot
dev.off()
