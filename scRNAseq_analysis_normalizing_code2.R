
################################# Gene-level filtering ################################################################################

# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 5

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]


# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)


#################################################### Normalize the data ################################################################


seurat_phase <- NormalizeData(filtered_seurat)

summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))

summary(seurat_phase@meta.data$riboRatio)


# Turn riboRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$riboFr <- cut(seurat_phase@meta.data$riboRatio, 
                                     breaks=c(-Inf, 0.01588, 0.0248, 0.039, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 30000, 
                                     verbose = FALSE)


top10 <- head(VariableFeatures(seurat_phase), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# Scale the counts
all.genes <- rownames(seurat_phase)
seurat_phase <- ScaleData(seurat_phase)

seurat_phase <- RunPCA(seurat_phase,  features = VariableFeatures(object = seurat_phase), ndims.print = 45)


# Plot the PCA colored by mitoFr
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")


# Plot the PCA colored by riboFr
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "riboFr",
        split.by = "riboFr")

#regress out the variation due to mitochondrial expression
#seurat_phase <- SCTransform(seurat_phase, vars.to.regress = c("mitoRatio"))

split_seurat <- SplitObject(seurat_phase, split.by = "sample")

split_seurat <- split_seurat[c("wt", "il22")]


options(future.globals.maxSize = 4000 * 1024^2)



for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio","riboRatio"), vst.flavor = "v2")
}





############################### donot run the code below ####################################################################


seurat_phase <- RunPCA(seurat_phase,  features = VariableFeatures(object = seurat_phase), ndims.print = 50)


# determine dimensionality of the data
ElbowPlot(seurat_phase, ndims = 50)


seurat_phase <- RunUMAP(seurat_phase, 
                        reduction = "pca", 
                        dims = 1:45)


8482 -> saved.seed
set.seed(saved.seed)


seurat_phase <- RunTSNE(
  seurat_phase,
  dims=1:45,
  seed.use = saved.seed, 
  perplexity=200,
  check_duplicates = FALSE
) 



# plot
p3 <- DimPlot(seurat_phase, reduction = 'tsne', group.by = 'sample')
p4 <- DimPlot(seurat_phase, reduction = 'umap', group.by = 'sample')

p3
p4



FindNeighbors(seurat_phase,dims=1:45) -> seurat_phase

FindClusters(seurat_phase,resolution = 1) -> seurat_phase

p5 <- DimPlot(seurat_phase,reduction = "tsne", pt.size = 0.5, split.by = "sample", label = TRUE) + ggtitle("tSNE with Perplexity 200")

p5


p5 <- DimPlot(seurat_phase,reduction = "umap", pt.size = 0.75, split.by = "sample", label = TRUE) + ggtitle("UMAP")

p5

p7 <- DimPlot(seurat_phase,reduction = "umap", pt.size = 0.75, group.by = "sample", label = TRUE) + ggtitle("UMAP")

p7





cluster_counts <- seurat_phase %>%
  group_by(sample, seurat_obj$ident) %>%
  summarize(CellCount = n())

# View the cluster counts
md <- seurat_phase@meta.data %>% as.data.table

cluster_cell_numbers <- md[, .N, by = c("sample", "seurat_clusters")]

cluster_cell_numbers <- as.data.frame(cluster_cell_numbers)

write.xlsx(cluster_cell_numbers, "C:/Users/aksha/Desktop/Soraya_snRNAseq/Outputs/PCA/cell_cluster_numbers_by_sample.xlsx")


