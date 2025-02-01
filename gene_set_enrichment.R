# Load necessary libraries
library(dplyr)  # Data manipulation
library(ggplot2)  # Data visualization
library(DESeq2)  # Differential expression analysis
library(clusterProfiler)  # GO enrichment analysis
library(org.Tsolium.eg.db)  # Annotation database for Taenia solium

# Create directories for output files
dir.create("ORA/blast2go/")
dir.create("ORA/blast2go/EV_PRE/")
dir.create("ORA/blast2go/EV_AT_SAT/")

# Load the count data
load("counts.RData")
data.counts <- as.data.frame(data.counts)

# Filter out rows with no counts
data.counts <- data.counts[which(rowSums(data.counts) > 0),]

# Define experimental conditions
condition <- factor(c("NC", "NC", "NC",
                      "PRE_SAT", "PRE_SAT", "PRE_SAT",
                      "PRE_AT", "PRE_AT", "PRE_AT",
                      "EV_SAT", "EV_SAT", "EV_SAT",
                      "EV_AT", "EV_AT", "EV_AT",
                      "POST_SAT", "POST_SAT", "POST_SAT",
                      "POST_AT", "POST_AT", "POST_AT"))

# Create a dataframe for conditions
coldata <- data.frame(row.names = colnames(data.counts), condition)

# Perform differential expression analysis (DESeq2)
res <- results(dds, contrast = c("condition", "EV_AT", "PRE_AT"))
res <- as.data.frame(res)

# Sort the results by the 'stat' value
res <- res %>% 
  arrange(-stat)

# Extract the gene list for GSEA
gene_list <- res$stat
names(gene_list) <- rownames(res)

# Perform Gene Set Enrichment Analysis (GSEA) for Biological Process (BP)
gse <- gseGO(gene_list,
             ont = "BP",
             keyType = "GID",
             OrgDb = "org.Tsolium.eg.db",
             eps = 1e-300)

# Convert the GSEA result into a data frame
gse_bar <- as.data.frame(gse)

# Add additional columns for visualization (log-transformed p-values and group)
gse_bar <- gse_bar %>% 
  mutate(log_p = -log10(p.adjust)) %>% 
  mutate(group = ifelse(enrichmentScore > 0, "EV", "PRE"))

# Ensure the log_p and enrichmentScore are numeric
gse_bar <- gse_bar %>%
  mutate(log_p = as.numeric(log_p),
         enrichmentScore = as.numeric(enrichmentScore))

# Define cutoff for significance (adjusted p-value < 0.05)
cutoff <- -log10(0.05)

# Create the bar plot for GSEA results
bar_plot <- ggplot(gse_bar, aes(x = reorder(Description, log_p), 
                                y = ifelse(enrichmentScore > 0, log_p, -log_p), 
                                fill = group)) +
  geom_col(width = 0.8) +
  coord_flip() +
  # Add dashed lines for p.adjust = 0.05
  geom_hline(yintercept = c(cutoff, -cutoff), 
             linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("EV" = "red", "PRE" = "gray")) +
  scale_y_continuous(
    limits = c(-6, 6),  # Set axis limits to -6 to 6
    breaks = seq(-6, 6, 2),  # Tick marks every 2 units
    labels = abs(seq(-6, 6, 2)),  # Use absolute values for labels
    expand = c(0, 0)
  ) +
  theme_minimal() +
  labs(
    x = NULL,
    y = expression(-log[10](p-value)),
    fill = NULL,
    title = "Enrichment Analysis"
  ) +
  theme(
    axis.text.y = element_text(size = 10),  # Customize y-axis text size
    axis.text.x = element_text(size = 10, color = "black"),  # X-axis labels in black
    axis.title.x = element_text(size = 12, color = "black"),  # X-axis title in black
    axis.ticks.x = element_line(color = "black"),  # X-axis ticks in black
    axis.line.x = element_line(color = "black"),  # X-axis line in black
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )

# Save the bar plot as SVG
ggsave("ORA/blast2go/EV_PRE/enrichment_plot.svg", 
       plot = bar_plot, width = 9, height = 12, units = "in", 
       dpi = 600, device = "svg")

# Generate and save specific GSEA plots for gene sets of interest
svg("ORA/blast2go/EV_PRE/mRNA_splicing.svg", height = 8, width = 12)
gseaplot(gse, geneSetID = 2) # Plot for mRNA splicing
dev.off()

svg("ORA/blast2go/EV_PRE/translation.svg", height = 8, width = 12)
gseaplot(gse, geneSetID = 5) # Plot for translation
dev.off()

svg("ORA/blast2go/EV_PRE/protein_folding.svg", height = 8, width = 12)
gseaplot(gse, geneSetID = 6) # Plot for protein folding
dev.off()

svg("ORA/blast2go/EV_PRE/monoatomic_transport.svg", height = 8, width = 12)
gseaplot(gse, geneSetID = 15) # Plot for monoatomic ion transport
dev.off()

svg("ORA/blast2go/EV_PRE/intracellular_signal_transduction.svg", height = 8, width = 12)
gseaplot(gse, geneSetID = 12) # Plot for intracellular signal transduction
dev.off()
