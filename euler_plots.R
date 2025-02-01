library(dplyr) # to manipulate dataframe content
library(DESeq2) # to get DEGs from dds object
library(eulerr) # for euler plots generation

# Load the counts file and filter out rows with zero counts
load("counts.RData")
data.counts <- data.counts[which(rowSums(data.counts) > 0),]

# Define conditions for comparisons and create a condition matrix
condition <- factor(c("NC", "NC", "NC",
                      "PRE_SAT", "PRE_SAT", "PRE_SAT",
                      "PRE_AT", "PRE_AT", "PRE_AT",
                      "EV_SAT", "EV_SAT", "EV_SAT",
                      "EV_AT", "EV_AT", "EV_AT",
                      "POST_SAT", "POST_SAT", "POST_SAT",
                      "POST_AT", "POST_AT", "POST_AT"))

coldata <- data.frame(row.names = colnames(data.counts), condition)

# Load the DESeq2 dds object
load("dds.RData")

# PRE vs NC comparisons
# Analyze PRE (+TA) vs NC
res_at <- results(dds, contrast = c("condition", "PRE_AT", "NC"))
res_at <- as.data.frame(res_at)

# Analyze PRE (-TA) vs NC
res_sat <- results(dds, contrast = c("condition", "PRE_SAT", "NC"))
res_sat <- as.data.frame(res_sat)

# Filter significant DEGs (padj < 0.05, |log2FC| > 1)
sig_at <- na.omit(res_at) %>% filter(padj < 0.05)
de_sig_at <- sig_at %>% filter(abs(log2FoldChange) > 1)

sig_sat <- na.omit(res_sat) %>% filter(padj < 0.05)
de_sig_sat <- sig_sat %>% filter(abs(log2FoldChange) > 1)

# Create gene sets for the +TA and -TA conditions
gene_sets <- list(
  "AT" = rownames(de_sig_at), # +TA gene set
  "SAT" = rownames(de_sig_sat) # -TA gene set
)

# Generate Euler plot for PRE vs NC comparison
fit <- euler(gene_sets)
svg("euler_pre_nc.svg", width = 11, height = 7)
plot(fit, fills = c("#3498db", "#e74c3c"), alpha = 0.5, labels = TRUE, quantities = TRUE, edges = TRUE)
dev.off()

# EV vs PRE comparisons
# Analyze EV vs PRE (+TA)
res_at <- results(dds, contrast = c("condition", "EV_AT", "PRE_AT"))
res_at <- as.data.frame(res_at)

# Analyze EV vs PRE (-TA)
res_sat <- results(dds, contrast = c("condition", "EV_SAT", "PRE_SAT"))
res_sat <- as.data.frame(res_sat)

# Filter significant DEGs (padj < 0.05, |log2FC| > 1)
sig_at <- na.omit(res_at) %>% filter(padj < 0.05)
de_sig_at <- sig_at %>% filter(abs(log2FoldChange) > 1)

sig_sat <- na.omit(res_sat) %>% filter(padj < 0.05)
de_sig_sat <- sig_sat %>% filter(abs(log2FoldChange) > 1)

# Create gene sets for the +TA and -TA conditions
gene_sets <- list(
  "AT" = rownames(de_sig_at), # +TA gene set
  "SAT" = rownames(de_sig_sat) # -TA gene set
)

# Generate Euler plot for EV vs PRE comparison
fit <- euler(gene_sets)
svg("euler_ev_pre.svg", width = 11, height = 7)
plot(fit, fills = c("#3498db", "#e74c3c"), alpha = 0.5, labels = TRUE, quantities = TRUE, edges = TRUE)
dev.off()

# POST vs EV comparisons
# Analyze POST vs EV (+TA)
res_at <- results(dds, contrast = c("condition", "POST_AT", "EV_AT"))
res_at <- as.data.frame(res_at)

# Analyze POST vs EV (-TA)
res_sat <- results(dds, contrast = c("condition", "POST_SAT", "EV_SAT"))
res_sat <- as.data.frame(res_sat)

# Filter significant DEGs (padj < 0.05, |log2FC| > 1)
sig_at <- na.omit(res_at) %>% filter(padj < 0.05)
de_sig_at <- sig_at %>% filter(abs(log2FoldChange) > 1)

sig_sat <- na.omit(res_sat) %>% filter(padj < 0.05)
de_sig_sat <- sig_sat %>% filter(abs(log2FoldChange) > 1)

# Create gene sets for the +TA and -TA conditions
gene_sets <- list(
  "AT" = rownames(de_sig_at), # +TA gene set
  "SAT" = rownames(de_sig_sat) # -TA gene set
)

# Generate Euler plot for POST vs EV comparison
fit <- euler(gene_sets)
svg("euler_post_ev.svg", width = 11, height = 7)
plot(fit, fills = c("#3498db", "#e74c3c"), alpha = 0.5, labels = TRUE, quantities = TRUE, edges = TRUE)
dev.off()

# +TA vs -TA comparisons across different stages
# Analyze POST +TA vs POST -TA
res_post <- results(dds, contrast = c("condition", "POST_AT", "POST_SAT"))
res_post <- as.data.frame(res_post)

# Analyze EV +TA vs EV -TA
res_ev <- results(dds, contrast = c("condition", "EV_AT", "EV_SAT"))
res_ev <- as.data.frame(res_ev)

# Analyze PRE +TA vs PRE -TA
res_pre <- results(dds, contrast = c("condition", "PRE_AT", "PRE_SAT"))
res_pre <- as.data.frame(res_pre)

# Filter significant DEGs (padj < 0.05, |log2FC| > 1)
sig_post <- na.omit(res_post) %>% filter(padj < 0.05)
de_sig_post <- sig_post %>% filter(abs(log2FoldChange) > 1)

sig_ev <- na.omit(res_ev) %>% filter(padj < 0.05)
de_sig_ev <- sig_ev %>% filter(abs(log2FoldChange) > 1)

sig_pre <- na.omit(res_pre) %>% filter(padj < 0.05)
de_sig_pre <- sig_pre %>% filter(abs(log2FoldChange) > 1)

# Create gene sets for POST, EV, and PRE conditions
gene_sets <- list(
  "POST" = rownames(de_sig_post),
  "EV" = rownames(de_sig_ev),
  "PRE" = rownames(de_sig_pre)
)

# Generate Euler plot for +TA vs -TA comparison across stages
fit <- euler(gene_sets)
svg("euler_at_sat.svg", width = 11, height = 7)
plot(fit, fills = c("#3498db", "#e74c3c", "#a3be5a"), alpha = 0.5, labels = TRUE, quantities = TRUE, edges = TRUE)
dev.off()

# Data usage example: Identify repressed and overexpressed genes in EV vs PRE comparisons
# Analyze EV vs PRE (+TA)
res_at <- results(dds, contrast = c("condition", "EV_AT", "PRE_AT"))
res_at <- as.data.frame(res_at)

# Analyze EV vs PRE (-TA)
res_sat <- results(dds, contrast = c("condition", "EV_SAT", "PRE_SAT"))
res_sat <- as.data.frame(res_sat)

# Filter significant DEGs (padj < 0.05, |log2FC| > 1)
sig_at <- na.omit(res_at) %>% filter(padj < 0.05)
de_sig_at <- sig_at %>% filter(abs(log2FoldChange) > 1)

sig_sat <- na.omit(res_sat) %>% filter(padj < 0.05)
de_sig_sat <- sig_sat %>% filter(abs(log2FoldChange) > 1)

# Identify repressed genes with log2FC < -1
de_sig_at_under <- sig_at %>% filter(log2FoldChange < -1)
de_sig_sat_under <- sig_sat %>% filter(log2FoldChange < -1)

# Create gene sets for repressed genes in both conditions
gene_sets <- list(
  "AT_under" = rownames(de_sig_at_under),
  "SAT_under" = rownames(de_sig_sat_under)
)

# Generate Euler plot for repressed genes
fit <- euler(gene_sets)
svg("euler_ev_pre_under.svg", width = 11, height = 7)
plot(fit, fills = c("#3498db", "#e74c3c"), alpha = 0.5, labels = TRUE, quantities = TRUE, edges = TRUE)
dev.off()

# Identify overexpressed genes with log2FC > 1
de_sig_at_over <- sig_at %>% filter(log2FoldChange > 1)
de_sig_sat_over <- sig_sat %>% filter(log2FoldChange > 1)

# Create gene sets for overexpressed genes in both conditions
gene_sets <- list(
  "AT_over" = rownames(de_sig_at_over),
  "SAT_over" = rownames(de_sig_sat_over)
)

# Generate Euler plot for overexpressed genes
fit <- euler(gene_sets)
svg("euler_ev_pre_over.svg", width = 11, height = 7)
plot(fit, fills = c("#3498db", "#e74c3c"), alpha = 0.5, labels = TRUE, quantities = TRUE, edges = TRUE)
dev.off()
