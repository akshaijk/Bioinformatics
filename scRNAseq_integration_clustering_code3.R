

split_seurat <- readRDS("C:/Users/aksha/Desktop/Soraya_snRNAseq/split_seurat.rds")




# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 



# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)




# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

## Don't run this during class
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")



# Scale data, run PCA and UMAP and visualize integrated data
seurat_integrated <- ScaleData(object = seurat_integrated)
seurat_integrated <- RunPCA(object = seurat_integrated)


ElbowPlot(seurat_phase, ndims = 50)


seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:45)

saved.seed  <-8400

seurat_integrated <- RunTSNE(
  seurat_integrated,dims=1:45,
  seed.use = saved.seed, 
  perplexity=200,
  check_duplicates = FALSE
) 


FindNeighbors(seurat_integrated,dims=1:45) -> seurat_integrated

FindClusters(seurat_integrated,resolution = 1) -> seurat_integrated

p9 <- DimPlot(seurat_integrated,reduction = "tsne", pt.size = 0.75, split.by = "sample", label = TRUE) + ggtitle("tSNE with Perplexity 200")

p9



p10 <- DimPlot(seurat_integrated,reduction = "umap", pt.size = 0.75, split.by = "sample", label = TRUE) + ggtitle("UMAP dims:45, resolution 1")

p10



p11 <- DimPlot(seurat_integrated,reduction = "tsne", pt.size = 0.75, label = TRUE) + ggtitle("tSNE with Perplexity 200")

p11



p12 <- DimPlot(seurat_integrated,reduction = "umap", pt.size = 0.75, label = TRUE) + ggtitle("umap")

p12




# View the cluster counts
md2 <- seurat_integrated@meta.data %>% as.data.table

integrated_cluster_cell_numbers <- md2[, .N, by = c("sample", "seurat_clusters")]

integrated_cluster_cell_numbers <- as.data.frame(integrated_cluster_cell_numbers)

write.csv(integrated_cluster_cell_numbers, "C:/Users/aksha/Desktop/Soraya_snRNAseq/rRNA_regressed_integrated_cluster_numbers.csv")




saveRDS(seurat_integrated, "C:/Users/aksha/Desktop/Soraya_snRNAseq/rRNA_seurat_integrated.rds")





################################################ dont run the code below ###############################################

# perform integration to correct for batch effects ------
obj.list <- SplitObject(seurat_phase, split.by = 'sample')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}



# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)


# find integration anchors (CCA)
integ_anchors <- FindIntegrationAnchors(object.list = obj.list,
                                        anchor.features = features)


# integrate data
seurat_integrated <- IntegrateData(anchorset = integ_anchors)

# Scale data, run PCA and UMAP and visualize integrated data
seurat_integrated <- ScaleData(object = seurat_integrated)
seurat_integrated <- RunPCA(object = seurat_integrated)


ElbowPlot(seurat_phase, ndims = 50)


seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:45)

saved.seed = 8624

seurat_integrated <- RunTSNE(
  seurat_integrated,dims=1:45,
  seed.use = saved.seed, 
  perplexity=200,
  check_duplicates = FALSE
) 


FindNeighbors(seurat_integrated,dims=1:45) -> seurat_integrated

FindClusters(seurat_integrated,resolution = 1) -> seurat_integrated

p9 <- DimPlot(seurat_integrated,reduction = "tsne", pt.size = 0.75, split.by = "sample", label = TRUE) + ggtitle("tSNE with Perplexity 200")

p9


p10 <- DimPlot(seurat_integrated,reduction = "umap", pt.size = 0.75, split.by = "sample", label = TRUE) + ggtitle("UMAP dims:45, resolution 1")

p10

p10 <- DimPlot(seurat_integrated,reduction = "umap", pt.size = 0.75, label = TRUE) + ggtitle("UMAP dims:45, resolution 1")

p10


p11 <- DimPlot(seurat_integrated,reduction = "tsne", pt.size = 0.75, label = TRUE) + ggtitle("tsne perplexity 200")

p11





# View the cluster counts
md2 <- seurat_integrated@meta.data %>% as.data.table

integrated_cluster_cell_numbers <- md2[, .N, by = c("sample", "seurat_clusters")]

integrated_cluster_cell_numbers <- as.data.frame(integrated_cluster_cell_numbers)

write_xlsx(integrated_cluster_cell_numbers, "C:/Users/aksha/Desktop/Soraya_snRNAseq/Outputs/PCA/integrated_cell_cluster_numbers_by_sample.xslx")

