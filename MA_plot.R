library(DESeq2)  # Load DESeq2 for differential expression analysis
library(dplyr)   # Load dplyr for data manipulation

# Load DESeq2 object containing RNA-seq data
load("dds.RData")

# Perform differential expression analysis for two conditions
resLFC1.true <- results(dds, contrast=c('condition', "EV_SAT", "PRE_SAT"))
resLFC2.true <- results(dds, contrast=c('condition', "EV_AT", "PRE_AT"))

# Generate MA plots and return data
LFC1.true <- plotMA(resLFC1.true, ylim=c(-3,3), alpha = 0.05, returnData=TRUE)
LFC2.true <- plotMA(resLFC2.true, ylim=c(-3,3), alpha = 0.05, returnData=TRUE)

# Add log-transformed mean expression values and significance indicators
LFC1.true <- LFC1.true %>%
  mutate(
    t.mean = log10(mean),  # Log10 transformation of mean expression
    SIG = abs(lfc) >= 1 & isDE  # Mark differentially expressed genes (DEG)
  )

LFC2.true <- LFC2.true %>%
  mutate(
    t.mean = log10(mean),
    SIG = abs(lfc) >= 1 & isDE
  )

# Apply shrinkage to log2 fold change values
resLFC1 <- lfcShrink(dds, contrast=c('condition', "EV_SAT", "PRE_SAT"), type="ashr")
resLFC2 <- lfcShrink(dds, contrast=c('condition', "EV_AT", "PRE_AT"), type="ashr")

# Generate MA plots with shrunken log2 fold changes
LFC1 <- plotMA(resLFC1, ylim=c(-3,3), alpha = 0.05, returnData=TRUE)
LFC2 <- plotMA(resLFC2, ylim=c(-3,3), alpha = 0.05, returnData=TRUE)

# Transfer significance labels from the unshrunken data
LFC1$SIG <- LFC1.true$SIG
LFC2$SIG <- LFC2.true$SIG

# Add log-transformed mean expression values
LFC1 <- LFC1 %>% mutate(t.mean = log10(mean))
LFC2 <- LFC2 %>% mutate(t.mean = log10(mean))

# Filtrar datos dentro y fuera del rango
LFC1_in_range <- LFC1[LFC1$lfc >= -3 & LFC1$lfc <= 3, ]
LFC1_out_range <- LFC1[LFC1$lfc < -3 | LFC1$lfc > 3, ]

LFC2_in_range <- LFC1[LFC2$lfc >= -3 & LFC2$lfc <= 3, ]
LFC2_out_range <- LFC2[LFC2$lfc < -3 | LFC2$lfc > 3, ]

# Save the MA plot data to CSV files
write.csv(LFC1, "PlotMA1.csv")
write.csv(LFC2, "PlotMA2.csv")

library(ggplot2)  # Load ggplot2 for visualization

# Create a SVG file for the MA plots
svg("plotMA_AT.svg", width=8, height=8)

# Generate MA plot for EV vs PRE (TA+)
ggplot() +
  # points in range
  geom_point(data = LFC2_in_range, aes(x = t.mean, y = lfc, color = SIG)) +
  
  # triangles if over/under ±3
  geom_point(data = LFC2_out_range, 
             aes(x = t.mean, y = ifelse(lfc > 3, 3, -3), color = SIG, 
                 shape = ifelse(lfc > 3, "up", "down")), 
             size = 3) +
  
  # up if over - down if under
  scale_shape_manual(values = c("up" = 2, "down" = 6)) +
  
  # colors
  scale_color_manual(values = c("TRUE" = "#5c5ffa", "FALSE" = "grey")) +
  
  # reference lines
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Threshold lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Baseline at 0
  
  # Y axis limits
  ylim(-3, 3) +
  
  # theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 24, face = "bold"),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  
  # labels
  labs(title = "EV vs PRE (TA+)", 
       x = "mean normalized counts (log10)", 
       y = "shrunk log2 fold change", 
       color = "DEG")

dev.off() # Close the SVG device

# Generate MA plot for EV vs PRE (TA+)
svg("plotMA_SAT.svg", width=8, height=8)

ggplot() +
  # points in range
  geom_point(data = LFC1_in_range, aes(x = t.mean, y = lfc, color = SIG)) +
  
  # triangles if over/under ±3
  geom_point(data = LFC1_out_range, 
             aes(x = t.mean, y = ifelse(lfc > 3, 3, -3), color = SIG, 
                 shape = ifelse(lfc > 3, "up", "down")), 
             size = 3) +
  
  # up if over - down if under
  scale_shape_manual(values = c("up" = 2, "down" = 6)) +
  
  # colors
  scale_color_manual(values = c("TRUE" = "#5c5ffa", "FALSE" = "grey")) +
  
  # reference lines
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Threshold lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Baseline at 0
  
  # Y axis limits
  ylim(-3, 3) +
  
  # theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 24, face = "bold"),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  
  # labels
  labs(title = "EV vs PRE (TA-)", 
       x = "mean normalized counts (log10)", 
       y = "shrunk log2 fold change", 
       color = "DEG")

dev.off()  # Close the SVG device

# Count differentially expressed genes (DEGs) in different conditions
DEG_count <- data.frame(
  comp=c("EV vs PRE (TA-)_neg", "EV vs PRE (TA-)_pos", "EV vs PRE (TA+)_neg", "EV vs PRE (TA+)_pos"),
  SIG_n =c(
    sum(LFC1$SIG == TRUE & LFC1$lfc <= 1), 
    sum(LFC1$SIG == TRUE & LFC1$lfc >= 1), 
    sum(LFC2$SIG == TRUE & LFC2$lfc <= 1), 
    sum(LFC2$SIG == TRUE & LFC2$lfc >= 1)
  )
)

# Print DEG counts
sum(LFC1$SIG == TRUE & LFC1$lfc <= 1)
sum(LFC1$SIG == TRUE & LFC1$lfc >= 1)
sum(LFC2$SIG == TRUE & LFC2$lfc <= 1)
sum(LFC2$SIG == TRUE & LFC2$lfc >= 1)

# Save DEG counts to CSV
write.csv(DEG_count, "DEG_count.csv")
