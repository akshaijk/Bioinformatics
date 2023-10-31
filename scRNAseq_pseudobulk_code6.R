
library(DESeq2)
library(SingleCellExperiment)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(pheatmap)
library(apeglm)
library(png)
library(RColorBrewer)
library(data.table)

seurat_integrated <- readRDS("C:/Users/aksha/Desktop/Soraya_snRNAseq/seurat_integrated_rRNA_regressed.rds")


# Rename all identities. Do not use underscores in names of clusters
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "0" = "absorptive enterocytes 1",
                                  "1" = "ileal enterocytes 1",
                                  "2" = "absorptive enterocytes 2",
                                  "3" = "absorptive enterocytes 3",
                                  "4" = "absorptive enterocytes 4",
                                  "5" = "ileal enterocytes 2",
                                  "6" = "cluster 6",
                                  "7" = "goblet cells 1",
                                  "8" = "cluster 8",
                                  "9" = "meschenchymal cells",
                                  "10" = "epithelia",
                                  "11" = "absorptive enterocytes",
                                  "12" = "immune cells ",
                                  "13" = "ionocytes",
                                  "14" = "entero endocrine cells 1",
                                  "15" = "goblet cells 2",
                                  "16" = "vascular smooth muscles",
                                  "17" = "entero endocrine cells 2",
                                  "18" = "neurons",
                                  "19" = "cluster 19",
                                  "20" = "pancreas",
                                  "21" = "meschenchymal",
                                  "22" = "liver",
                                  "23" = "neurons",
                                  "24" = "cluster 24"
                                  )

counts <- seurat_integrated@assays$RNA@counts 

metadata <- seurat_integrated@meta.data


# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat_integrated@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)


## Check the assays present
assays(sce)


## Check the counts matrix
dim(counts(sce))
counts(sce)[1:6, 1:6]


# Extract unique names of clusters (= levels of cluster_id factor variable)
cluster_names <- levels(colData(sce)$cluster_id)
cluster_names

# Total number of clusters
length(cluster_names)

# Extract unique names of samples (= levels of sample_id factor variable)

sample_names <- c('wt', 'il22')


sample_names

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("cluster_id", "sample")]

head(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum") 

aggr_counts


# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]



# Using which() to look up tstrsplit() output
goblet_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "goblet cells 1")
goblet_idx

colnames(aggr_counts)[goblet_idx]
aggr_counts[1:10, goblet_idx]

##################################################### code works till this ########################################################################



# As a reminder, we stored our cell types in a vector called cluster_names
cluster_names


# Loop over all cell types to extract corresponding counts, and store information in a list

## Initiate empty list
counts_ls <- list()


for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)

##################################################### code works till this ########################################################################



# Reminder: explore structure of metadata
head(colData(sce))


# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(sample)

dim(metadata)
head(metadata)

##################################################### code works till this ########################################################################


# Number of cells per sample and cluster
t <- table(colData(sce)$sample,
           colData(sce)$cluster_id)

t


##################################################### code works till this ########################################################################


#Creating metadata list

## Initiate empty list
metadata_ls <- list()



for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster id with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, as.data.frame(metadata), 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}


# Explore the different components of the list
str(metadata_ls)


# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))

# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))



##################################################### code works till this ########################################################################






idx <- which(names(counts_ls) == "absorptive enterocytes 1")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]


cluster_counts
head(cluster_metadata)

# Check matching of matrix columns and metadata rows
all(colnames(cluster_counts) == rownames(cluster_metadata))


# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ sample_id)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample_id")


# Generate results object
res <- results(dds, 
               name = "group_id_stim_vs_ctrl",
               alpha = 0.05)




# Assuming you want to compare different groups
dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, design = ~ sample_id)
dds <- DESeq(dds)

resultsNames(dds)


