

seurat_integrated <- readRDS("C:/Users/aksha/Desktop/Soraya_snRNAseq/seurat_integrated_rRNA_regressed.rds")

DefaultAssay(seurat_integrated) <- "RNA"

annotations <- gene_annotations_zebrafish_csv


# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

conserved_markers <- map_dfr(c(0:24), get_conserved)


# Extract top 20 markers per cluster
top20 <- conserved_markers %>% 
  mutate(avg_fc = (wt_avg_log2FC + il22_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)

top20 <- top10

write.csv(top20, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/top20_genes_rRNA_reg_2.csv", row.names = TRUE)



seurat_integrated %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend()



head(conserved_markers)

length(conserved_markers)

gene_list_pedro <- c('rbp2a', 'fabp1b.1', 'apoa1a', 'apoa4b.2.1', 'apobb.1',
                     'scarb1', 'chia.1', 'chia.2', 'mbl2', 'rida', 'gcshb',
                     'agr2', 'SPDEF', 'galnt7',
                     'pou2f3', 'spry2',
                     'tnfrsf11a', 'tnfrsf11b',
                     'neurod1', 'pax6b', 'nkx2.2a', 'pcsk1', 'gck', 'scg2b', 'scg3', 'scg5',
                     'mlnl',
                     'cftr', 'tmem51a', 'stap2b',
                     'best4', 'gucy2c', 'npr1a',
                     'phox2bb', 'numbl',
                     'desma',
                     'desmb',
                     'pdgfra', 'vim', 
                     'bcam', 'cxcl8b.1', 'rhag',
                     'aoc2', 'agtr2', 'emilin1a',
                     'mpeg1.1', 'marco',
                     'mpx',
                     'lck', 'nitr9', 'nitr4a', 'ccr9a', 'il7r',
                     'angptl6', 'agtr2'
)

# Remove duplicated genes
gene_list_pedro <- unique(gene_list_pedro)


my_genes <- c(
  "ucp1", "fabp6", "fabp2", "apoc2", "pck1", "slc13a2", "chia.2", 
  "cubn", "dab2", "amn", "lrp2b",
  "gck", "scg3", 
  "pax6b", "nkx2.2a", "sstr5", "cnfn", "icn", "agr2", "SPDEF", "stap2b", "cftr", "tmem51a", 
  "cldn15la", "vil1", "prox1a", "egr1", "fabp10a", "gc", "tfa", "apom", "prss1", "ela2l", 
  "gcga", "tagln", "acta2", "myl9a", "phox2bb", "numbl", "vgf", "nmu", "calb2a", 
  "otop2", "ca4b", "desmb", "mpeg1.1", "pfn1", "rgs2", "pou2f3", "trpm5", "vim", "bcam", 
  "aoc2", "esama", "lgr4", "lrig1"
)

my_dotplot <- DotPlot(seurat_integrated, features = my_genes, cols="RdBu") + RotatedAxis()
traceback()

my_dotplot

pedro_dotplot <- DotPlot(seurat_integrated, features = gene_list_pedro, cols="RdBu") + RotatedAxis()
traceback()

pedro_dotplot

# genes to test DE
genes_to_DE = c(
  "il22", "stat3", "tph1b", "trpa1b", "neurod1", "gcga", "pyyb",
   "vipb", "cdx1b", "ptger4c", "cldnb", "numbl", "phox2bb", "ntrk1", "neflb",
   "ache", "adrb2b", "desma", "desmb", "myh11a", "adora2b", "tagln",
  "pdgfaa", "tpm4b", "cd36", "fabp2",  "cldn2", "cftr", "dab2", "agt"
)


# violin plot in loop 
p <- list()  # Use list() to initialize an empty list

for (i in gene_list_pedro) {
  x <- VlnPlot(object = seurat_integrated, fhttp://127.0.0.1:9971/graphics/b0fe4900-8684-4a44-8c80-5101a2ecda5b.pngeatures = i)
  p[[i]] <- x  # Append the plot to the list using the gene name as the list element name
}


for (j in seq_along(p)) {
  print(p[[j]])
}

plot_dir <- ('C:/Users/aksha/Desktop/Soraya_snRNAseq/final_ppt/figs/expression_of_genes_vlnplot_pedrosuggested')

for (gene_name in names(p)) {
  # Print the plot
  print(p[[gene_name]])
  
  # Define the filename (use paste0 to concatenate the directory and gene name)
  filename <- file.path(plot_dir, paste0(gene_name, ".png"))
  
  # Save the plot as an image file (adjust file type if needed)
  png(filename)
  print(p[[gene_name]])
  dev.off()
}


# DE violin plot in loop 
p <- list()  # Use list() to initialize an empty list

for (i in genes_to_DE) {
  x <- VlnPlot(object = seurat_integrated, features = i, split.by = "sample")
  p[[i]] <- x  # Append the plot to the list using the gene name as the list element name
}


for (j in seq_along(p)) {
  print(p[[j]])
}

plot_dir <- ('C:/Users/aksha/Desktop/Soraya_snRNAseq/final_ppt/figs/DE_genes')

for (gene_name in names(p)) {
  # Print the plot
  print(p[[gene_name]])
  
  # Define the filename (use paste0 to concatenate the directory and gene name)
  filename <- file.path(plot_dir, paste0(gene_name, ".png"))
  
  # Save the plot as an image file (adjust file type if needed)
  png(filename)
  print(p[[gene_name]])
  dev.off()
}


# feature plot in loop 
f <- list()  # Use list() to initialize an empty list

for (i in gene_list_pedro) {
  y <- FeaturePlot(object = seurat_integrated,reduction= 'tsne', features = i, max.cutoff = 2.5,
                   cols = c("grey", "red"), label = TRUE)
  f[[i]] <- y  # Append the plot to the list using the gene name as the list element name
}


for (j in seq_along(f)) {
  print(f[[j]])
}




plot_dir <- ('C:/Users/aksha/Desktop/Soraya_snRNAseq/final_ppt/figs/pedro_genes_featureplot')

for (gene_name in names(f)) {
  # Print the plot
  print(f[[gene_name]])
  
  # Define the filename (use paste0 to concatenate the directory and gene name)
  filename <- file.path(plot_dir, paste0(gene_name, ".png"))
  
  # Save the plot as an image file (adjust file type if needed)
  png(filename,  width = 8 * 300, height = 8 * 300, res = 300)
  print(f[[gene_name]])
  dev.off()
}











FeaturePlot(seurat_integrated, reduction= 'tsne',  features = c("fabp6"), max.cutoff = 2.5,
            cols = c("grey", "red"), label = TRUE)










# distinguishing cluster markers
cluster0.markers <- FindMarkers(seurat_integrated, ident.1 = 0, ident.2 = c(2, 3), min.pct = 0.10)
head(cluster0.markers, n = 20)

all_markers <- FindAllMarkers(seurat_integrated, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)


write.csv(all_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/all_markers_rRNA_reg.csv", row.names = TRUE)






# find conserved markers


cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample")

cluster0_conserved_markers <- as.data.frame(cluster0_conserved_markers)
cluster0_conserved_markers$average_log2FC <- rowMeans(cluster0_conserved_markers[, c("il22_avg_log2FC", "wt_avg_log2FC")])
cluster0_conserved_markers <- cluster0_conserved_markers[order(-cluster0_conserved_markers$average_log2FC), ]
  
write_tsv(cluster0_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/integrated_find_conserved_markers/cluster0_conserved_markers.tsv")

cluster0_conserved_markers$average_log2FC <- rowMeans(cluster0_conserved_markers[, c("il22_avg_log2FC", "wt_avg_log2FC")])
cluster0_conserved_markers <- cluster0_conserved_markers[order(-cluster0_conserved_markers$average_log2FC), ]


cluster1_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                   ident.1 = 1,
                                                   grouping.var = "sample")

write.csv(cluster1_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster1_conserved_markers.csv", row.names = TRUE)




cluster2_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                   ident.1 = 2,
                                                   grouping.var = "sample")

write.csv(cluster2_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster2_conserved_markers.csv", row.names = TRUE)



cluster3_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                   ident.1 = 3,
                                                   grouping.var = "sample")

write.csv(cluster3_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster3_conserved_markers.csv", row.names = TRUE)


cluster4_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                   ident.1 = 4,
                                                   grouping.var = "sample")

write.csv(cluster4_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster4_conserved_markers.csv", row.names = TRUE)


cluster5_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                   ident.1 = 5,
                                                   grouping.var = "sample")

write.csv(cluster5_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster5_conserved_markers.csv", row.names = TRUE)



cluster6_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                   ident.1 = 6,
                                                   grouping.var = "sample")

write.csv(cluster6_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster6_conserved_markers.csv", row.names = TRUE)




cluster7_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                   ident.1 = 7,
                                                   grouping.var = "sample")

write.csv(cluster7_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster7_conserved_markers.csv", row.names = TRUE)




cluster8_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                   ident.1 = 8,
                                                   grouping.var = "sample")

write.csv(cluster8_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster8_conserved_markers.csv", row.names = TRUE)




cluster9_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                   ident.1 = 9,
                                                   grouping.var = "sample")

write.csv(cluster9_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster9_conserved_markers.csv", row.names = TRUE)




cluster10_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 10,
                                                    grouping.var = "sample")

write.csv(cluster10_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster10_conserved_markers.csv", row.names = TRUE)



cluster11_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 11,
                                                    grouping.var = "sample")

write.csv(cluster11_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster11_conserved_markers.csv", row.names = TRUE)


cluster12_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 12,
                                                    grouping.var = "sample")

write.csv(cluster12_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster12_conserved_markers.csv", row.names = TRUE)



cluster13_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 13,
                                                    grouping.var = "sample")

write.csv(cluster13_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster13_conserved_markers.csv", row.names = TRUE)




cluster14_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 14,
                                                    grouping.var = "sample")

write.csv(cluster14_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster14_conserved_markers.csv", row.names = TRUE)



cluster15_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 15,
                                                    grouping.var = "sample")

write.csv(cluster15_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster15_conserved_markers.csv", row.names = TRUE)


cluster15_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 15,
                                                    grouping.var = "sample")

write.csv(cluster15_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster15_conserved_markers.csv", row.names = TRUE)




cluster16_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 16,
                                                    grouping.var = "sample")

write.csv(cluster16_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster16_conserved_markers.csv", row.names = TRUE)



cluster17_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 17,
                                                    grouping.var = "sample")

write.csv(cluster17_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster17_conserved_markers.csv", row.names = TRUE)




cluster18_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 18,
                                                    grouping.var = "sample")

write.csv(cluster18_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster18_conserved_markers.csv", row.names = TRUE)




cluster19_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 19,
                                                    grouping.var = "sample")

write.csv(cluster19_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster19_conserved_markers.csv", row.names = TRUE)



cluster20_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 20,
                                                    grouping.var = "sample")

write.csv(cluster20_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster20_conserved_markers.csv", row.names = TRUE)




cluster21_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 21,
                                                    grouping.var = "sample")

write.csv(cluster21_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster21_conserved_markers.csv", row.names = TRUE)




cluster22_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 22,
                                                    grouping.var = "sample")

write.csv(cluster22_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster22_conserved_markers.csv", row.names = TRUE)



cluster23_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 23,
                                                    grouping.var = "sample")

write.csv(cluster23_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster23_conserved_markers.csv", row.names = TRUE)



cluster24_conserved_markers <- FindConservedMarkers(seurat_phase,
                                                    ident.1 = 24,
                                                    grouping.var = "sample")

write.csv(cluster24_conserved_markers, file = "C:/Users/aksha/Desktop/Soraya_snRNAseq/Find_Conserved_Markers/cluster24_conserved_markers.csv", row.names = TRUE)


