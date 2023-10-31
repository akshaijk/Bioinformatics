library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)
library(data.table)
library(tidyverse)
library(RCurl)
library(cowplot)
library(DelayedArray)
library(S4Arrays)
library(Biobase)
library(zlibbioc)
library(glmGamPoi)
library(XVector)
library(GenomicRanges)
library(GenomeInfoDb)
library(GenomeInfoDbData)
library(IRanges)
library(S4Vectors)
library(MatrixGenerics)
library(SummarizedExperiment)
library(BiocGenerics)
library(generics)
library(multtest)
library(glmGamPoi)
library(Matrix)



wt <- Read10X(data.dir = "C:/Users/aksha/Desktop/Soraya_snRNAseq/raw_data/WT")
il22 <- Read10X(data.dir = "C:/Users/aksha/Desktop/Soraya_snRNAseq/raw_data/il22")

wt_seurat <- CreateSeuratObject(counts = wt, project = "wt", min.cells = 3, min.features = 100)
il22_seurat <- CreateSeuratObject(counts = il22, project = "il22", min.cells = 3, min.features = 100)


merged_seurat <- merge(x = wt_seurat, y = il22_seurat,
                       add.cell.ids = c("wt_seurat", "il22_seurat")
)


View(merged_seurat@meta.data)      # SEE THE METADATA

merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100


# Compute percent ribo ratio
merged_seurat$riboRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^rp[sl][[:digit:]]|^rplp[[:digit:]]|^rpsa")
merged_seurat$riboRatio <- merged_seurat@meta.data$riboRatio / 100



# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Create a new column in the metadata dataframe called cells
metadata$cells <- rownames(metadata)


# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^wt_"))] <- "wt"
metadata$sample[which(str_detect(metadata$cells, "^il22_"))] <- "il22"
traceback()


# Rename columns
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                seq_folder = orig.ident,
                nGene = nFeature_RNA)



# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

View(merged_seurat@meta.data)


VlnPlot(merged_seurat, features = c("nUMI", "nGene", "mitoRatio", "riboRatio"), ncol = 4)



# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  geom_text(
    aes(label = ..count..),  # Display count values on bars
    stat = "count",          # Use the count statistic
    vjust = -0.5,            # Adjust the vertical position of labels
    size = 3                 # Adjust the font size of labels
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 400) +
  geom_vline(xintercept = 4000)


# Visualize the distribution of nGenes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.5)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=riboRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.1)


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)


# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 400) &
                            (nUMI <= 4000) &
                            (nGene >= 300) &
                            (log10GenesPerUMI > 0.75) & 
                            (mitoRatio < 0.10))



filtered_seurat <- subset(x = filtered_seurat, 
                          subset= (riboRatio < 0.10))  

filtered_seurat@meta.data

filtered_metadata <- filtered_seurat@meta.data 


VlnPlot(filtered_seurat, features = c("nUMI", "nGene", "mitoRatio", "riboRatio"), ncol = 4)




# Visualize the number of cell counts per sample
filtered_metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  geom_text(
    aes(label = ..count..),  # Display count values on bars
    stat = "count",          # Use the count statistic
    vjust = -0.5,            # Adjust the vertical position of labels
    size = 3                 # Adjust the font size of labels
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


filtered_metadata %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  geom_text(
    aes(label = ..count..),  # Display count values on bars
    stat = "count",          # Use the count statistic
    vjust = -0.5,            # Adjust the vertical position of labels
    size = 3                 # Adjust the font size of labels
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells")


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
filtered_metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

saveRDS(filtered_seurat, "C:/Users/aksha/Desktop/Soraya_snRNAseq/ribo_filtered_seurat.rds")

