library(dplyr) # to manipulate data
library(DESeq2) # to retrieve DEGs from a dds object
library(tidyplots) # to visualize DEGs as volcano plots

# Load the count data and remove genes with zero counts in any sample
load("counts.RData")
data.counts <- data.counts[which(rowSums(data.counts) > 0),]

# Create a condition matrix for the comparisons
condition <- factor(c("NC", "NC", "NC",
                      "PRE_SAT", "PRE_SAT", "PRE_SAT",
                      "PRE_AT", "PRE_AT", "PRE_AT",
                      "EV_SAT", "EV_SAT", "EV_SAT",
                      "EV_AT", "EV_AT", "EV_AT",
                      "POST_SAT", "POST_SAT", "POST_SAT",
                      "POST_AT", "POST_AT", "POST_AT"))

coldata <- data.frame(row.names = colnames(data.counts), condition)

# Load the previously obtained DESeq2 dds object
load("dds.RData")

# Data Usage: Visualize differentially expressed genes with volcano plots
# EV_AT vs PRE_AT comparison
res_at <- results(dds, contrast = c("condition", "EV_AT", "PRE_AT"))
res_at <- as.data.frame(res_at)
res_at <- res_at %>% 
  mutate(
    neg_log10_padj = -log10(padj),  # Calculate the negative log10 of adjusted p-values
    direction = if_else(log2FoldChange > 0, "up", "down", NA),  # Classify log2FoldChange as up or down
    candidate = abs(log2FoldChange) >= 1 & padj < 0.05,  # Identify significant candidates
    gene_id = rownames(res_at)  # Add gene IDs for labeling
  )

# EV_SAT vs PRE_SAT comparison
res_sat <- results(dds, contrast = c("condition", "EV_SAT", "PRE_SAT"))
res_sat <- as.data.frame(res_sat)
res_sat <- res_sat %>% 
  mutate(
    neg_log10_padj = -log10(padj),  # Calculate the negative log10 of adjusted p-values
    direction = if_else(log2FoldChange > 0, "up", "down", NA),  # Classify log2FoldChange as up or down
    candidate = abs(log2FoldChange) >= 1 & padj < 0.05,  # Identify significant candidates
    gene_id = rownames(res_sat)  # Add gene IDs for labeling
  )

# Generate volcano plots

# EV vs PRE (+TA) volcano plot
tidy_volcano <- res_at %>% 
  tidyplot(x = log2FoldChange, y = neg_log10_padj) %>% 
  add_data_points(data = filter_rows(!candidate), color = "lightgrey", rasterize = TRUE) %>% 
  add_data_points(data = filter_rows(candidate, direction == "up"), color = "#FF7777", alpha = 0.5) %>% 
  add_data_points(data = filter_rows(candidate, direction == "down"), color = "#7DA8E6", alpha = 0.5) %>% 
  add_reference_lines(x = c(-1, 1), y = -log10(0.05)) %>% 
  add_data_labels_repel(data = min_rows(padj, 5, by = direction), label = gene_id, color = "#000000", 
                        min.segment.length = 0, background = TRUE) %>% 
  adjust_x_axis_title("$Log[2]~fold~change$") %>% 
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") %>% 
  adjust_size(width = NA, height = NA) %>% 
  save_plot("ev_pre_at_volcano.svg", width = 16, height = 12, units = "cm")

# EV vs PRE (-TA) volcano plot
tidy_volcano <- res_sat %>% 
  tidyplot(x = log2FoldChange, y = neg_log10_padj) %>% 
  add_data_points(data = filter_rows(!candidate), color = "lightgrey", rasterize = TRUE) %>% 
  add_data_points(data = filter_rows(candidate, direction == "up"), color = "#FF7777", alpha = 0.5) %>% 
  add_data_points(data = filter_rows(candidate, direction == "down"), color = "#7DA8E6", alpha = 0.5) %>% 
  add_reference_lines(x = c(-1, 1), y = -log10(0.05)) %>% 
  add_data_labels_repel(data = min_rows(padj, 5, by = direction), label = gene_id, color = "#000000", 
                        min.segment.length = 0, background = TRUE) %>% 
  adjust_x_axis_title("$Log[2]~fold~change$") %>% 
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") %>% 
  adjust_size(width = NA, height = NA) %>% 
  save_plot("tidy_ev_pre_sat_v.svg", width = 16, height = 12, units = "cm")
