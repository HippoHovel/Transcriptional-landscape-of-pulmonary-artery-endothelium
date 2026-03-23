#check available data
DimPlot(ECstory.j, cols=varibow_pal, label=T, split.by = "samplecondition")
DimPlot(ECstory.j, cols=varibow_pal, label=T, group.by = "cluster.identity")
DimPlot(ECstory.j, cols=varibow_pal, label=T, group.by = "seurat_clusters")

#confirm Endothelium
FeaturePlot_scCustom(ECstory.j, features=c("PECAM1"))

#Endosubset
ECstory.subset <- subset(ECstory.j, seurat_clusters %in% c("6","12"))

# mytmp <- "/home/isilon/users/o_lins/temp"
# dir.create(mytmp, showWarnings = FALSE, recursive = TRUE)
# # Override tempdir() for the session
# assignInNamespace("tempdir", function() mytmp, ns = "base")

DimPlot(ECstory.subset)

ECstory.subset <- NormalizeData(ECstory.subset)
ECstory.subset <- FindVariableFeatures(ECstory.subset)
ECstory.subset <- ScaleData(ECstory.subset)
ECstory.subset <- RunPCA(ECstory.subset)

ECstory.subset <- FindNeighbors(ECstory.subset, dims = 1:30, reduction = "pca")
ECstory.subset <- FindClusters(ECstory.subset, resolution = 0.04, cluster.name = "unintegrated_clusters")

ECstory.subset <- RunUMAP(ECstory.subset, reduction = "integrated.cca", dims = 1:40, reduction.name = "umap.cca")

DimPlot(ECstory.subset)
DimPlot(ECstory.subset, split.by = "samplecondition")

ECstory.subset$seurat_clusters <- Idents(ECstory.subset) 


#MAST for GOI
EC.cluster.recognition <- FindAllMarkers(ECstory.subset, test.use = "MAST")

top <- EC.cluster.recognition

top <- rownames_to_column(top, var="gene")

top$avg_log2FC <- as.numeric(top$avg_log2FC)



#check the top 20 cluster defining genes
top20 <- top
top20$valid = ""
top20$valid [(top20$p_val_adj < 0.05)] <- TRUE #check for significance
top20$valid[(top20$pct.1 < 0.2)] <- FALSE #pick genes that the cluster is expressing to a high %
top20$valid[(top20$pct.2 > 0.2)] <- FALSE #dismiss genes that the other populations are also expressing to a high %
top20 <- subset(top20, subset = top20$valid == TRUE )

top20 <- top20 %>% group_by(cluster) %>% top_n(20, rank) #take the top20 genes per cluster
top20 <- top20[order(top20$cluster),]
top20 <- na.omit(top20)

DotPlot_scCustom(ECstory.subset, features = c(rev(top20$gene)), flip_axes = T, x_lab_rotate = T, colors_use=c("#4f4ceb", "#807ee9", "#a8a7e6","#d8d8d8","#f38d8d","#f25656","#da0000"))+scale_size(range = c(1, 5))

plot_heatmap(dataset=ECstory.subset,
             markers = top20$gene,
             sort_var = c("seurat_clusters"),
             anno_var = c("seurat_clusters"),
             anno_colors= list(polychrome_pal),
             row_font_size = 7.5)

#EC QC
#check cluster 4
cluster4 <- subset(EC.cluster.recognition, cluster %in% c("4"))
cluster4

FeaturePlot_scCustom(ECstory.j, features=c("TYROBP", "FTL", "FCER1G", "C1QA"))
FeaturePlot_scCustom(ECstory.j, features=c("AIF1", "RBM47", "MSR1", "FTH1"))
FeaturePlot_scCustom(ECstory.j, features=c("LYZ", "C5AR1", "ACSL1", "C1orf162"))
FeaturePlot_scCustom(ECstory.j, features=c("IFI30", "FMN1", "SLC16A10", "MS4A7"))
FeaturePlot_scCustom(ECstory.j, features=c("CD163", "FCGR2A", "LST1", "KYNU"))
#cluster 4 shows heavy resemblance to Monocytes => filter out



#check cluster 3
cluster3 <- subset(EC.cluster.recognition, cluster %in% c("3"))
cluster3

FeaturePlot_scCustom(ECstory.j, features=c("CXCR4", "CCL5", "SYTL3", "STAT4"))
FeaturePlot_scCustom(ECstory.j, features=c("PTPRC", "PARP8", "CD69", "GZMA"))
FeaturePlot_scCustom(ECstory.j, features=c("ARHGAP15", "CD96", "CD3D", "RUNX3"))
FeaturePlot_scCustom(ECstory.j, features=c("RGS1", "CST7", "AOAH", "CD3E"))
FeaturePlot_scCustom(ECstory.j, features=c("GNG2", "SKAP1", "RHOH", "NKG7"))
FeaturePlot_scCustom(ECstory.j, features=c("TRBC2", "GZMH", "TRBC1", "KLRB1"))
FeaturePlot_scCustom(ECstory.j, features=c("TRBC1", "CTSW", "GZMK", "BCL11B"))
FeaturePlot_scCustom(ECstory.j, features=c("CD3E", "CD247", "BCL11B", "ICOS"))
FeaturePlot_scCustom(ECstory.subset, features="PTPRC", split.by = "seurat_clusters")#3 and 4 very high in CD45
#cluster 3 is overwhelmingly PTPRC positive and shows heavy resemblance to T-cells => filter out


#check cluster 2
cluster2 <- subset(EC.cluster.recognition, cluster %in% c("2"))
cluster2

FeaturePlot_scCustom(ECstory.j, features=c("FCN3", "ADGRL2", "AFF3", "STXBP6"))
FeaturePlot_scCustom(ECstory.j, features=c("BTNL9", "RGCC", "ADGRF5", "HECW2"))
FeaturePlot_scCustom(ECstory.j, features=c("HPGD", "SORBS1", "CA4", "PRICKLE2"))
FeaturePlot_scCustom(ECstory.j, features=c("PRKCE", "GALNT18", "TMEM100", "IL1RL1"))
#keep


#check cluster 1
cluster1 <- subset(EC.cluster.recognition, cluster %in% c("1"))
cluster1

FeaturePlot_scCustom(ECstory.j, features=c("PDE3A", "DKK2", "IGFBP3", "PTGIS"))
FeaturePlot_scCustom(ECstory.j, features=c("FN1", "BMX", "SULF1", "ABI3BP"))
FeaturePlot_scCustom(ECstory.j, features=c("NLGN1", "CXCL2", "HMCN1", "SOX5"))
FeaturePlot_scCustom(ECstory.j, features=c("GJA5", "SETBP1", "EDN1", "BGN"))
FeaturePlot_scCustom(ECstory.j, features=c("CPAMD8", "CP", "DEPP1", "SERPINE2"))
FeaturePlot_scCustom(ECstory.j, features=c("MEIS2", "CCN2", "CPE", "COL8A1"))
FeaturePlot_scCustom(ECstory.j, features=c("PLPP1", "EXT1", "OMD", "LTBP4"))
FeaturePlot_scCustom(ECstory.j, features=c("MECOM", "RBMS3", "CRIM1", "CLU"))
FeaturePlot_scCustom(ECstory.j, features=c("PTPRB", "MGP", "BMP6", "RHOB"))
#keep


#check cluster 0
cluster0 <- subset(EC.cluster.recognition, cluster %in% c("0"))
cluster0

FeaturePlot_scCustom(ECstory.j, features=c("ZNF385D", "TSHZ2", "ACKR1", "PLVAP"))
FeaturePlot_scCustom(ECstory.j, features=c("SPRY1", "ADAMTS9", "SNTG2", "POSTN"))
FeaturePlot_scCustom(ECstory.j, features=c("COL15A1", "CCL14", "DIPK2B", "EBF1"))
FeaturePlot_scCustom(ECstory.j, features=c("EMCN", "ADGRL4", "TLL1", "FABP5"))
FeaturePlot_scCustom(ECstory.j, features=c("PKP4", "AL109930.1", "MEOX1", "OLFM1"))
FeaturePlot_scCustom(ECstory.j, features=c("SPARCL1", "RPL10A", "RPL12", "RPS4X"))
FeaturePlot_scCustom(ECstory.j, features=c("RPS2", "RPL26", "RPS6", "RPS18"))
#keep

unique(Idents(ECstory.subset))
unique(ECstory.subset$seurat_clusters)
ECstory.subset$seurat_clusters <- Idents(ECstory.subset)
ECstory.subset.fil <- subset(ECstory.subset, seurat_clusters %!in% c("3", "4"))#filtering out non endothelial cells

DimPlot(ECstory.subset.fil, group.by = "seurat_clusters")

DimPlot(ECstory.subset.fil, group.by = "seurat_clusters") +
  scale_x_reverse() +
  scale_y_reverse()

#MAST DEGs with filtered ECs only
EC.cluster.recognition2 <- FindAllMarkers(ECstory.subset.fil, test.use = "MAST")

top2 <- EC.cluster.recognition2

#identify top genes from subclusters
#same filter as above
top2$valid = ""
top2$valid [(top2$p_val_adj < 0.05)] <- TRUE #check for significance
top2$valid[(top2$pct.1 < 0.2)] <- FALSE #pick genes that the cluster is expressing to a high %
top2$valid[(top2$pct.2 > 0.2)] <- FALSE #dismiss genes that the other populations are also expressing to a high %
top2 <- subset(top2, subset = top2$valid == TRUE )

top2 <- top2 %>% group_by(cluster) %>%top_n(1000, rank)
top2 <- top2[order(top2$cluster),]
top2 <- na.omit(top2)
#top.c1.2 <- top20

#top genes by cluster:
top.c0 <- subset(top2, cluster %in% c("0"))  
top.c1 <- subset(top2, cluster %in% c("1"))
top.c2 <- subset(top2, cluster %in% c("2"))

top.c0.2 <- subset(top.c0, avg_log2FC >= 0.5)
top.c1.2 <- subset(top.c1, avg_log2FC >= 0.5)
top.c2.2 <- subset(top.c2, avg_log2FC >= 0.5)

#sort only for genes that are significant in each population
top.sorted <- EC.cluster.recognition2 
top.sorted <- subset(top.sorted , gene %in% c(top.c0$gene,top.c1$gene,top.c2$gene))

top.sorted.c0 <- subset(top.sorted , cluster %in% c("0"))
top.sorted.c1 <- subset(top.sorted , cluster %in% c("1"))
top.sorted.c2 <- subset(top.sorted , cluster %in% c("2"))  



#set colours, ranks and labels  
top.sorted.c0$colour = ifelse((top.sorted.c0$avg_log2FC>= 0.5) & (top.sorted.c0$p_val_adj <= 0.05 & top.sorted.c0$p_val_adj >=0), "#b20505", 
                              ifelse((top.sorted.c0$avg_log2FC <= -0.5) & (top.sorted.c0$p_val_adj <= 0.05), "black", "lightgray"))
top.sorted.c0$rank <- abs(top.sorted.c0$avg_log2FC * -log10(top.sorted.c0$p_val_adj))
top.sorted.c0$labels = ifelse(top.sorted.c0$gene %in% head(top.sorted.c0[order(top.sorted.c0$rank, decreasing = T),], n = 30)$gene, ifelse(top.sorted.c0$colour != "gray",top.sorted.c0$gene, ""), "")

top.sorted.c1$colour = ifelse((top.sorted.c1$avg_log2FC>= 0.5) & (top.sorted.c1$p_val_adj <= 0.05 & top.sorted.c1$p_val_adj >=0), "#0312ad", 
                              ifelse((top.sorted.c1$avg_log2FC <= -0.5) & (top.sorted.c1$p_val_adj <= 0.05), "black", "lightgray"))
top.sorted.c1$rank <- abs(top.sorted.c1$avg_log2FC * -log10(top.sorted.c1$p_val_adj))
top.sorted.c1$labels = ifelse(top.sorted.c1$gene %in% head(top.sorted.c1[order(top.sorted.c1$rank, decreasing = T),], n = 30)$gene, ifelse(top.sorted.c1$colour != "gray",top.sorted.c1$gene, ""), "")

top.sorted.c2$colour = ifelse((top.sorted.c2$avg_log2FC>= 0.5) & (top.sorted.c2$p_val_adj <= 0.05 & top.sorted.c2$p_val_adj >=0), "#03ad46", 
                              ifelse((top.sorted.c2$avg_log2FC <= -0.5) & (top.sorted.c2$p_val_adj <= 0.05), "black", "lightgray"))
top.sorted.c2$rank <- abs(top.sorted.c2$avg_log2FC * -log10(top.sorted.c2$p_val_adj))
top.sorted.c2$labels = ifelse(top.sorted.c2$gene %in% head(top.sorted.c2[order(top.sorted.c2$rank, decreasing = T),], n = 30)$gene, ifelse(top.sorted.c2$colour != "gray",top.sorted.c2$gene, ""), "")

#Volcano Plots
ggplot(top.sorted.c0, aes(x = top.sorted.c0$avg_log2FC, y = -log10(top.sorted.c0$p_val_adj)))+
  geom_point(color = top.sorted.c0$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")  

ggplot(top.sorted.c1, aes(x = top.sorted.c1$avg_log2FC, y = -log10(top.sorted.c1$p_val_adj)))+
  geom_point(color = top.sorted.c1$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")

ggplot(top.sorted.c2, aes(x = top.sorted.c2$avg_log2FC, y = -log10(top.sorted.c2$p_val_adj)))+
  geom_point(color = top.sorted.c2$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")
#Volcanoplots have deceptive pvalues (overinfalted due to single cell test) -> go for MA plot  

#MA plot
#add average expression
geneLogSums <- log2(rowSums(GetAssayData(ECstory.subset.fil, "RNA")))
top.sorted.c0$logSum <- geneLogSums[top.sorted.c0$gene]
top.sorted.c1$logSum <- geneLogSums[top.sorted.c1$gene]
top.sorted.c2$logSum <- geneLogSums[top.sorted.c2$gene]

#isolate significant fingerprint gene for each subcluster
top.sorted.c0.FP <- subset(top.sorted.c0, colour == c("#b20505"))
top.sorted.c1.FP <- subset(top.sorted.c1, colour == c("#0312ad"))
top.sorted.c2.FP <- subset(top.sorted.c2, colour == c("#03ad46"))

ggplot(top.sorted.c0.FP, aes(x=logSum, y=avg_log2FC))+
  geom_text_repel(label = top.sorted.c0.FP$labels, nudge_y = 0.2, na.rm = FALSE, size=2.6, max.overlaps = 100) +
  geom_point(color = top.sorted.c0.FP$colour, alpha = 0.5, show.legend = TRUE)+
  theme_bw()+
  theme_linedraw() +
  xlab("Average expression") +
  ylab("Fold change (log2)")


ggplot(top.sorted.c1.FP, aes(x=logSum, y=avg_log2FC))+
  geom_text_repel(label = top.sorted.c1.FP$labels, nudge_y = 0.2, na.rm = FALSE, size=2.6, max.overlaps = 100) +
  geom_point(color = top.sorted.c1.FP$colour, alpha = 0.5, show.legend = TRUE)+
  theme_bw()+
  theme_linedraw() +
  xlab("Average expression") +
  ylab("Fold change (log2)")

ggplot(top.sorted.c2.FP, aes(x=logSum, y=avg_log2FC))+
  geom_text_repel(label = top.sorted.c2.FP$labels, nudge_y = 0.2, na.rm = FALSE, size=2.6, max.overlaps = 100) +
  geom_point(color = top.sorted.c2.FP$colour, alpha = 0.5, show.legend = TRUE)+
  theme_bw()+
  theme_linedraw() +
  xlab("Average expression") +
  ylab("Fold change (log2)")


#Enrichr for GSEA

library(devtools)
library(enrichR)
library(reshape2)
library(rrvgo)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

#define databases of interest
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022")

#rerun from here for different clusters
#run enrichr on ckuster sepcific genes
if (websiteLive) {
  enriched <- enrichr(top.sorted.c0.FP$gene, dbs)#switch with top.sorted.c1.FP$gene, top.sorted.c2.FP$gene for each cluster
}

#check with classic graph
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value")
}


enriched$GO_Biological_Process_2023$Count <- str_extract(enriched$GO_Biological_Process_2023$Overlap, ".+?(?=/)")
enriched$GO_Biological_Process_2023$Count <- as.numeric(enriched$GO_Biological_Process_2023$Count)

enrichedplot <- head(enriched$GO_Biological_Process_2023, 50)
enrichedplot2 <- enrichedplot

# enrichedplot <- head(enriched$GO_Biological_Process_2023, 50)
# enrichedplot2 <- enrichedplot

# ggplot(enrichedplot2, aes(x=reorder(Term, rev(P.value)), y=Count)) +
#   geom_bar(stat='identity', fill=("#c6841e")) +
#   geom_text(aes(label = format(P.value,scientific=T,digits=2)), hjust = 1.0,vjust = -0.6, colour = "#0d6d1d", size=3)+
#   geom_text(aes(label = Overlap), hjust = 1.0,vjust = 1.6, colour = "#3b3b3b", size=3)+
#   #geom_text(aes(label = Genes),hjust=1.0, colour="black", size=2.5)+
#   theme_bw()+
#   coord_flip()#+
# theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

#choose GOs curated to avoid redundancy
chosenGO.c0 <- c("Positive Regulation Of Protein Localization To Chromosome, Telomeric Region (GO:1904816)",
                 "Regulation Of Protein Localization To Cajal Body (GO:1904869)",
                 "Actin Filament Organization (GO:0007015)",
                 "Extracellular Matrix Organization (GO:0030198)",
                 "Protein Stabilization (GO:0050821)",
                 "Long-Chain Fatty Acid Transport (GO:0015909)",
                 "Positive Regulation Of Leukocyte Cell-Cell Adhesion (GO:1903039)",
                 "Positive Regulation Of T Cell Activation (GO:0050870)",
                 "MHC Class II Protein Complex Assembly (GO:0002399)",
                 "Negative Regulation Of Epithelial Cell Migration (GO:0010633)")

chosenGO.c1 <- c("Regulation Of Smooth Muscle Cell Proliferation (GO:0048660)",
                 "Regulation Of Cell Migration (GO:0030334)",
                 "Regulation Of Angiogenesis (GO:0045765)",
                 "Extracellular Matrix Assembly (GO:0085029)",
                 "Positive Regulation Of MAPK Cascade (GO:0043410)",
                 "Negative Regulation Of Cell Population Proliferation (GO:0008285)",
                 "Aortic Valve Morphogenesis (GO:0003180)",
                 "Positive Regulation Of Vasculature Development (GO:1904018)",
                 "Response To Hydrogen Peroxide (GO:0042542)",
                 "Negative Regulation Of Cellular Response To Growth Factor Stimulus (GO:0090288)")

chosenGO.c2 <- c("Positive Regulation Of Cyclin-Dependent Protein Kinase Activity (GO:1904031)",
                 "Positive Regulation Of DNA-binding Transcription Factor Activity (GO:0051091)",
                 "Regulation Of Small GTPase Mediated Signal Transduction (GO:0051056)",
                 "Positive Regulation Of Cell Cycle (GO:0045787)",
                 "Endothelin Receptor Signaling Pathway (GO:0086100)",
                 "Regulation Of Intracellular Signal Transduction (GO:1902531)",
                 "Negative Regulation Of Cell Adhesion (GO:0007162)",
                 "Positive Regulation Of NF-kappaB Transcription Factor Activity (GO:0051092)",
                 "Leukocyte Aggregation (GO:0070486)",
                 "Regulation Of Sodium Ion Transmembrane Transporter Activity (GO:2000649)")




enrichedplot3 <- subset(enrichedplot2, Term %in% chosenGO.c0)# switch with chosenGO.c1 and chosenGO.c2 for realted clusters


#Chrd diagram to show the genes of top curated ontologies
chordgenes <- subset(enrichedplot3, select=c("Term","Genes","Adjusted.P.value"))

chorddata <- na.omit(unlist(strsplit(chordgenes$Genes, ';')))
chorddata2 <-  chorddata[!duplicated(chorddata)]

chordgenes2 <- separate_rows(chordgenes, Genes, sep = ";")
chordgenes2$Adjusted.P.value <- -log10(chordgenes2$Adjusted.P.value)


gridcolors <- c("#000000","#8e4806","#8e7906","#6f8e06", "#488e06", 
                "#8e0619", "#068e71", "#06638e", "#50068e","#8e0686",
                ifelse(chorddata2 != 0, "#03ad46", "#5A61B7")) #switch True from ifelse statement with "#b20505" for cluster0, "#0312ad" for cluster 1 and "#03ad46" for cluster2

circos.par(start.degree = 270)
chordDiagram(chordgenes2, scale = FALSE, annotationTrack = c("name","grid"), link.lty = 1,annotationTrackHeight = c(0.03, 0.01))#, grid.col = gridcolors)
chordDiagram(chordgenes2, scale = FALSE, annotationTrack = c("grid"), link.lty = 1, grid.col = gridcolors)
circos.clear()


ECstory.subset.fil$subset.identity <- ""
ECstory.subset.fil$subset.identity[ECstory.subset.fil$seurat_clusters == "0"] <- "Immuno modulatory"
ECstory.subset.fil$subset.identity[ECstory.subset.fil$seurat_clusters == "1"] <- "Vascular tone modulatory"
ECstory.subset.fil$subset.identity[ECstory.subset.fil$seurat_clusters == "2"] <- "Vacular plasticity modulatory"

ECstory.subset.fil$namedclusters <- ECstory.subset.fil$subset.identity

#supp


DimPlot(ECstory.j, group.by = "cluster.identity", split.by = "samplecondition", cols= polychrome_pal)
DimPlot(ECstory, group.by = "cluster.identity", split.by = "samplecondition", cols= polychrome_pal)



p <- FeaturePlot_scCustom(
  ECstory.subset.fil,
  features = "ACTA2",
  colors_use = c("lightgrey", "red"),
  order = TRUE
)

p + scale_color_gradientn(
  colors = c("lightgrey", "red"),
  limits = c(0, 5),
  oob = scales::squish
)


p <- FeaturePlot_scCustom(
  ECstory.subset.fil,
  features = "LYVE1",
  colors_use = c("lightgrey", "red"),
  order = TRUE
)

p + scale_color_gradientn(
  colors = c("lightgrey", "red"),
  limits = c(0, 5),
  oob = scales::squish
)



