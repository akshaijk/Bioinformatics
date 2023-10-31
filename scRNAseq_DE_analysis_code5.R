seurat_integrated <- readRDS("C:/Users/aksha/Desktop/Soraya_snRNAseq/seurat_integrated_rRNA_regressed.rds")

DefaultAssay(seurat_integrated) <- "RNA"


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


DefaultAssay(seurat_integrated) <- "RNA"





seurat_integrated$celltype.condition <- paste(Idents(seurat_integrated), seurat_integrated$sample, sep="_")

seurat_integrated$celltype.condition 

seurat_integrated$celltype <- Idents(seurat_integrated)

Idents(seurat_integrated) <- "celltype.condition"



il22_response_goblet_cell_1 <- FindMarkers(seurat_integrated, ident.1 = "goblet cells 1_il22" , ident.2 = "goblet cells 1_wt", verbose = FALSE)
il22_response_goblet_cell_1
il22_response_goblet_cell_1 <- as.data.frame(il22_response_goblet_cell_1)
write.csv(il22_response_goblet_cell_1, "C:/Users/aksha/Desktop/Soraya_snRNAseq/il22_response_goblet_cell_1.csv")



EnhancedVolcano(il22_response_goblet_cell_1 , 
                rownames(il22_response_goblet_cell_1 ),
                x ="avg_log2FC", 
                y ="p_val_adj")



il22_response_absorptive_enterocytes_1 <- FindMarkers(seurat_integrated, ident.1 = "absorptive enterocytes 1_il22" , ident.2 = "absorptive enterocytes 1_wt", verbose = FALSE)
il22_response_absorptive_enterocytes_1
il22_response_absorptive_enterocytes_1 <- as.data.frame(il22_response_absorptive_enterocytes_1)
write.csv(il22_response_absorptive_enterocytes_1, "C:/Users/aksha/Desktop/Soraya_snRNAseq/il22_response_absorptive_enterocytes_1.csv")

EnhancedVolcano(il22_response_absorptive_enterocytes_1 , 
                rownames(il22_response_absorptive_enterocytes_1 ),
                x ="avg_log2FC", 
                y ="p_val_adj")


il22_response_absorptive_enterocytes_2 <- FindMarkers(seurat_integrated, ident.1 = "absorptive enterocytes 2_il22" , ident.2 = "absorptive enterocytes 2_wt", verbose = FALSE)
il22_response_absorptive_enterocytes_2
il22_response_absorptive_enterocytes_2 <- as.data.frame(il22_response_absorptive_enterocytes_2)
write.csv(il22_response_absorptive_enterocytes_2, "C:/Users/aksha/Desktop/Soraya_snRNAseq/il22_response_absorptive_enterocytes_2.csv")


EnhancedVolcano(il22_response_absorptive_enterocytes_2 , 
                rownames(il22_response_absorptive_enterocytes_2 ),
                x ="avg_log2FC", 
                y ="p_val_adj")

il22_response_absorptive_enterocytes_3 <- FindMarkers(seurat_integrated, ident.1 = "absorptive enterocytes 3_il22" , ident.2 = "absorptive enterocytes 3_wt", verbose = FALSE)
il22_response_absorptive_enterocytes_3
il22_response_absorptive_enterocytes_3 <- as.data.frame(il22_response_absorptive_enterocytes_2)
write.csv(il22_response_absorptive_enterocytes_3, "C:/Users/aksha/Desktop/Soraya_snRNAseq/il22_response_absorptive_enterocytes_3.csv")




EnhancedVolcano(il22_response_absorptive_enterocytes_3 , 
                rownames(il22_response_absorptive_enterocytes_3 ),
                x ="avg_log2FC", 
                y ="p_val_adj")




il22_response_absorptive_enterocytes_4 <- FindMarkers(seurat_integrated, ident.1 = "absorptive enterocytes 4_il22" , ident.2 = "absorptive enterocytes 4_wt", verbose = FALSE)
il22_response_absorptive_enterocytes_4
il22_response_absorptive_enterocytes_4 <- as.data.frame(il22_response_absorptive_enterocytes_4)
write.csv(il22_response_absorptive_enterocytes_4, "C:/Users/aksha/Desktop/Soraya_snRNAseq/il22_response_absorptive_enterocytes_4.csv")




EnhancedVolcano(il22_response_absorptive_enterocytes_4 , 
                rownames(il22_response_absorptive_enterocytes_4 ),
                x ="avg_log2FC", 
                y ="p_val_adj")




il22_response_absorptive_enterocytes_4 <- FindMarkers(seurat_integrated, ident.1 = "absorptive enterocytes 4_il22" , ident.2 = "absorptive enterocytes 4_wt", verbose = FALSE)
il22_response_absorptive_enterocytes_4
il22_response_absorptive_enterocytes_4 <- as.data.frame(il22_response_absorptive_enterocytes_4)
write.csv(il22_response_absorptive_enterocytes_4, "C:/Users/aksha/Desktop/Soraya_snRNAseq/il22_response_absorptive_enterocytes_4.csv")




EnhancedVolcano(il22_response_absorptive_enterocytes_4 , 
                rownames(il22_response_absorptive_enterocytes_4 ),
                x ="avg_log2FC", 
                y ="p_val_adj")




ileal_enterocytes_2 <- subset(seurat_integrated, idents = "absorptive enterocytes 1")
Idents(absorptive_enterocytes_1) <- "wt"

avg_absorptive_enterocytes_1 <- as.data.frame(log1p(AverageExpression(absorptive_enterocytes_1, verbose = FALSE)$RNA))
avg_absorptive_enterocytes_1$gene <- rownames(avg_absorptive_enterocytes_1)


# genes to test DE
genes.to.label = c(
  "il22", "il26", "stat3", "tph1b", "trpa1b", "neurod1", "gcga", "pyyb", "pyya",
  "nkx2.2", "vipb", "cdx1b", "ptger4c", "cldnb", "numbl", "phox2bb", "ntrk1", "neflb",
  "chrna", "ache", "cspr1a", "adrb2b", "desma", "desmb", "myh11a", "adora2b", "tagln",
  "pdgfaa", "tpm4b", "cd36", "fabp2", "slc66a1", "cldn2", "cftr", "dab2", "meis2b", "agt"
)



p2 <- ggplot(avg_absorptive_enterocytes_1, aes('wt', 'il22')) + geom_point() + ggtitle("absorptive enterocytes 1")
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)


p2


