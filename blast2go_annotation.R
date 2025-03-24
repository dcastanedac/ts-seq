library(AnnotationForge) # for annotation
library(AnnotationDbi) # for annotation
library(dplyr) # for data manipulation

# Define the path to the GO annotations file
path <- "GO annotations/GO_annotations.txt"

# Read the GO annotations file into a dataframe
background <- read.csv(path, sep = "\t")

# Separate rows and columns to parse GO information into distinct fields
background <- background %>% 
  separate_rows(GO.in.Extended.Format..GO.Category.GO.ID.GO.Term., 
                sep = ",") %>% 
  separate(GO.in.Extended.Format..GO.Category.GO.ID.GO.Term.,
           into = c("GO_Category", "GO_ID", "GO_Term"), 
           sep = " ", extra = "merge", fill = "right")

# Remove the 'Enzyme.Codes' column as it's not needed for further analysis
background <- background %>% 
  dplyr::select(-Enzyme.Codes)

# Filter out rows with missing or invalid GO categories
background <- background %>% 
  filter(!is.na(GO_Category) & GO_Category != "") %>%
  filter(GO_Category %in% c("P","F","C","p","f","c")) %>% 
  # Recode GO categories for consistency
  mutate(GO_Category = recode(GO_Category, 
                              P = "BP", 
                              F = "MF", 
                              C = "CC",
                              p = "BP", 
                              f = "MF", 
                              c = "CC"))

# Rename the columns for clarity
colnames(background) <- c("Gene.stable.ID","Description","GO.domain", 
                          "GO.term.accession", "GO.term.name")

# Create a unique gene information table with distinct gene IDs and descriptions
gene_info <- background %>%
  dplyr::select(GeneID = Gene.stable.ID, Description) %>%
  dplyr::distinct()

# Rename the columns for consistency
colnames(gene_info) <- c("GID", "gene_name")

# Create unique gene names by appending a number to duplicates
gene_info_unique <- gene_info %>%
  group_by(gene_name) %>%
  mutate(unique_name = paste(gene_name, row_number(), sep = "_")) %>%
  ungroup() %>%
  dplyr::select(GID, unique_name) %>%
  dplyr::rename(gene_name = unique_name)

# Update the gene_info dataframe with the unique names
gene_info <- gene_info_unique

# Extract the GO terms and associated information
go_terms <- background %>%
  dplyr::select(GID = Gene.stable.ID, Term = GO.term.name, 
                Ontology = GO.domain) %>%
  dplyr::distinct()

# Rename the columns for clarity
colnames(go_terms) <- c("GID", "TERM", "ONTOLOGY")

# Create a table mapping gene IDs to GO terms
gene2go <- background %>%
  dplyr::select(GeneID = Gene.stable.ID, GOID = GO.term.accession, 
                GOALL = GO.term.accession) %>%
  dplyr::distinct()

# Rename the columns for consistency
colnames(gene2go) <- c("GID", "GO", "GOALL")

# Create directories for storing the annotation output
dir.create("annotation_blast2go/")
dir.create("annotation_blast2go/tsolium/")

# Generate the OrgDb package for Taenia solium annotations
makeOrgPackage(
  gene_info = gene_info,
  go = go_terms,
  gene2go = gene2go,
  version = "0.1",
  maintainer = "David Castaneda-Carío <david.castaneda.c@upch.pe>",
  author = "David Castaneda-Carpio",
  outputDir = "annotation_blast2go/tsolium",
  genus = "Taenia",
  species = "solium",
  tax_id = "6204"
)

# Install the generated OrgDb package for Taenia solium
install.packages("annotation_blast2go/tsolium/org.Tsolium.eg.db", 
                 repos = NULL, 
                 type = "source")

