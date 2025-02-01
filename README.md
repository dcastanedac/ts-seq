# From cyst to tapeworm: Transcriptome analysis with R of *Taenia solium in-vitro* activation
The following repository contains the scripts used to validate and analyze *Taenia solium* transcriptome. These can be used as an starting point to explore the data and for hypothesis making.
The pipeline consists on the use of clustering analysis for validation of experimental groups, and filtering and visualization of differentially expressed genes (DEGs).
In addition, we generated an annotation package that can be used for over representation analysis (ORA) of functional categories in the dataset.

## Software and Packages used:
- R version 4.4.2
  - tidyverse: for data cleaning, manipulation (dplyr) and visualization (ggplot2)
  - DESeq2: for differential gene expression (DGE) analysis
  - tidyplots: for data visualization
  - eulerr: for visualization of DEGs across experimental conditions
  - annotationForge and annotationDbi: to generate an annotation package based on the data obtained with OmicsBox.
  - clusterProfiler: for ORA with Gene Set Enrichment Analysis (GSEA)
- RStudio 2024.12.0 Build 467
- OmicsBox version 3.3.2

## Samples and Conditions Description:
Each sample was labeled according to each individual's category and their treatment. **AT** and **SAT** makes reference to **+TA** and **-TA** respectively.
- **NC** makes reference to **non-cultured** cysts.
- **PRE_SAT** makes reference to **non-evaginated individuals cultured for 6h** in RPMI medium **-TA**. 
- **PRE_AT** makes reference to **non-evaginated individuals cultured for 6h** in RPMI medium **+TA**.
- **EV_SAT** makes reference to **recently evaginated indiviuals** (around the 12-24h range) in RPMI medium **-TA**. 
- **EV_AT** makes reference to **recently evaginated individuals** (around 12-24h range) in RPMI medium **+TA**.
- **POST_SAT** makes reference to **evaginated indiviuals kept in culture for 120h** in RPMI medium **-TA**. 
- **POST_AT** makes reference to **evaginated individuals kept in culture for 120h** in RPMI medium **+TA**.

## Scripts usage:
### For Technical Validation and Basic Data Exploration:
- **DGE_PCA_hierarchical_clustering.R:** Performs DGE analysis, PCA, and hierarchical clustering to validate and explore structure across replicates for each experimental conditions.
- **euler_plots.R:** Generates Euler diagrams to visualize comparisons of DEGs across each experimental condition, particularly -TA and +TA.

### As an example of further analysis for EV vs PRE (+TA):
- **volcano_plots.R:** Creates volcano plots for visualizing DGE results and to identify possible interesting genes.
- **blast2go_annotation.R:** Creates and install an annotation object of the data obtained from OmicsBox with assigned Gene Ontology (GO) terms for each gene.
- **gene_set_enrichment.R:** Performs GSEA using the previously created annotation object for all the GO terms belonging to the Biological Process (BP) ontology.

## GO Annotations and Annotations Package Folders:
The GO Annotations folder contain the following files:
- **GO_annotations.txt** is a text file obtained from OmicsBox used in the **blast2go_annotation.R** script to create an annotation object.
- **OmicsBoxProj_tsolium-annotation_details.pdf** is a file with details of the annotation project from the OmicsBox suite.
- **Tsolium_annotation_full-export.xlsx** is a spreadsheet containing the full results of the annotation in the OmicsBox suite.

The Annotations Package folder contain the files created by annotationForge, and can be used to install the annotation package used in this project.
Use the following code:
>> *install.packages("org.Tsolium.eg.db", repos = NULL, type = "source")*

## Contact Information:
Be free to contact David Castañeda to the following email in case of any question or concern: *david.castaneda.c@upch.pe*
