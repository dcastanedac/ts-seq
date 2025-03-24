# Open R and make sure you have the following packages installed:
# Rsubread, GenomicFeatures, DESeq2, ggplot2, tidyverse, dplyr, purrr, openxlsx, pheatmap, ggrepel, fgsea, biomaRt, tibble
# You will need to replace the directory '/home/user/Documents/wd' with your actual working directory
main_dir <- '/home/user/Documents/wd'

##### Build index for the Taenia solium genome
# Create a directory to store the index files, if it doesn't already exist
dir.create(file.path(main_dir, 'taenia_solium_index'), showWarnings = FALSE)

# Build the index for the Taenia solium genome
# .fa file was downloaded from WormBase ParaSite
buildindex(basename="taenia_solium_index/tsolium", reference="taenia_solium.PRJNA170813.WBPS18/taenia_solium.PRJNA170813.WBPS18.genomic.fa")

###### Identify trimmed files
# Create a directory to store trimmed data if it doesn't already exist
dir.create(file.path(main_dir, 'trimmed_data'), showWarnings = FALSE)
trimmed_dir <- "./trimmed_data"

# List paired-end read files (R1 and R2) from the trimmed data directory
reads1 <- list.files(path = file.path(trimmed_dir), pattern = "*1.fq.gz$", full.names = TRUE)
reads2 <- list.files(path = file.path(trimmed_dir), pattern = "*2.fq.gz$", full.names = TRUE)

##### Map reads to the Taenia solium genome using Rsubread package
library(Rsubread)

# Perform the read alignment using the Taenia solium index
Rsubread::align(index = "taenia_solium_index/tsolium",
      readfile1 = reads1,
      readfile2 = reads2,
      type = "rna", # RNA-seq data
      input_format = "gzFASTQ", # Input format is gzipped FASTQ
      output_format = "BAM", # Output format is BAM
      PE_orientation = "rf", # Paired-end orientation is RF
      nthreads = 24 ) # Use 24 threads for processing

# List BAM files generated from the alignment
bam.files <- list.files(trimmed_dir, pattern = "BAM$", full.names=TRUE)
bam.files

###### Get the feature data to obtain the counts for each feature
# First, make a TxDb object from the GFF file, which contains gene annotations
# GFF file was downloaded from WormBase ParaSite
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("taenia_solium.PRJNA170813.WBPS18/taenia_solium.PRJNA170813.WBPS18.annotations.gff3.gz")

# Extract exon data from the TxDb object by gene
ex <- exonsBy(txdb, "gene")

# Generate a SAF (Seq Annotation File) data.frame from the exon data
len <- lengths(ex)  # Get the lengths of each exon
ex2 <- unlist(ex)    # Flatten the exon data
saf <- data.frame(GeneID = rep(names(len), len), Chr = seqnames(ex2), Start = start(ex2), End = end(ex2), Strand = as.character(strand(ex2)))

# Perform feature counting using the featureCounts function
fc.str <- featureCounts(files = bam.files, annot.ext = saf, isPairedEnd = TRUE, requireBothEndsMapped = TRUE, nthreads = 24)

##### Get counts matrix and reorder the columns
# Extract the raw counts from the featureCounts output
data.counts <- fc.str$counts

# Reorder the columns to match the experimental conditions
data.counts <- data.counts[,c(10,11,12,19,20,21,7,8,9,13,14,15,1,2,3,16,17,18,4,5,6)]

# Assign meaningful column names based on the experimental conditions
colnames(data.counts) <-  c('CTRL_1', 'CTRL_2','CTRL_3',
				            'PRE_EV_SAT_1', 'PRE_EV_SAT_2', 'PRE_EV_SAT_3',
				            'PRE_EV_AT_1', 'PRE_EV_AT_2', 'PRE_EV_AT_3',
				            'EV_SAT_1', 'EV_SAT_2','EV_SAT_3',
				            'EV_AT_1', 'EV_AT_2', 'EV_AT_3',
				            'POST_EV_SAT_1', 'POST_EV_SAT_2', 'POST_EV_SAT_3',
				            'POST_EV_AT_1', 'POST_EV_AT_2', 'POST_EV_AT_3')

# Save the counts matrix as an RData file for future use
save(data.counts, file = 'counts.RData')

