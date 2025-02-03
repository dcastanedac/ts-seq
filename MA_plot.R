library("DESeq2")

# Load DESeq2 results object
load("dds.RData")

# Perform log-fold change (LFC) shrinkage using the "ashr" method for two contrasts
resLFC1 <- lfcShrink(dds, contrast=c('Factor', ev_sat, pre_ev_sat), type="ashr")
resLFC2 <- lfcShrink(dds, contrast=c('Factor', ev_at, pre_ev_at), type="ashr")

# Generate MA plots and store the plot data for further processing
LFC1 <- plotMA(resLFC1, ylim=c(-3,3), alpha = 0.05, returnData=TRUE)
LFC2 <- plotMA(resLFC2, ylim=c(-3,3), alpha = 0.05, returnData=TRUE)

# Save the MA plot data as CSV files
write.csv(LFC1, "PlotMA1.csv")
write.csv(LFC2, "PlotMA2.csv")

# Load dplyr for data manipulation
library(dplyr)

# Add transformed mean expression and significance flag (lfc >= 1 and significant DE genes)
LFC1 <- LFC1 %>%
  mutate(
    t.mean = log10(mean),  # Log10 transformation of mean normalized counts
    SIG = abs(lfc) >= 1 & isDE  # Define significance based on log-fold change and differential expression
  )

LFC2 <- LFC2 %>%
  mutate(
    t.mean = log10(mean),  # Log10 transformation of mean normalized counts
    SIG = abs(lfc) >= 1 & isDE  # Define significance based on log-fold change and differential expression
  )

# Load ggplot2 for visualization
library(ggplot2)

# Open a PDF file to save plots
pdf("plotMA.pdf", width=11, height=8, onefile = TRUE)

# Generate first MA plot
ggplot(LFC1, aes(x = t.mean, y = lfc, color = SIG)) +
  geom_point() +
  scale_color_manual(values = c("TRUE" = "#DDA0DD", "FALSE" = "grey")) +  # Custom colors for significant/non-significant points
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Threshold lines at log2FC = -1 and 1
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Reference line at log2FC = 0
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border to plot
    axis.text = element_text(size = 20),  # Increase axis text size
    axis.title = element_text(size = 24, face = "bold"),  # Increase and bold axis title size
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # Centered, bold plot title
    legend.text = element_text(size = 20),  # Increase legend text size
    legend.title = element_text(size = 20, face = "bold")  # Bold legend title
  ) +
  labs(title = "EV vs PRE (TA-)", x = "log10 of mean normalized counts", y = "log2 of fold change (lfc)", color = "padj <= 0.05 and lfc >= 1")

# Generate second MA plot
ggplot(LFC2, aes(x = t.mean, y = lfc, color = SIG)) +
  geom_point() +
  scale_color_manual(values = c("TRUE" = "#DDA0DD", "FALSE" = "grey")) +  # Custom colors for significant/non-significant points
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Threshold lines at log2FC = -1 and 1
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Reference line at log2FC = 0
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border to plot
    axis.text = element_text(size = 20),  # Increase axis text size
    axis.title = element_text(size = 24, face = "bold"),  # Increase and bold axis title size
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # Centered, bold plot title
    legend.text = element_text(size = 20),  # Increase legend text size
    legend.title = element_text(size = 20, face = "bold")  # Bold legend title
  ) +
  labs(title = "EV vs PRE (TA-)", x = "log10 of mean normalized counts", y = "log2 of fold change (lfc)", color = "padj <= 0.05 and lfc >= 1")

# Close the PDF file
dev.off()
