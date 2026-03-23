set.seed(1)
#create UMAP with stricter seperation for pseudotime
ECstory2 <- ECstory

ECstory2 <- FindVariableFeatures(ECstory2, nfeatures = 4000)
ECstory2 <- ScaleData(ECstory2)
ECstory2 <- RunPCA(ECstory2)


ECstory2 <- IntegrateLayers(
  object = ECstory2, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE,
  seed.use =1
)


ECstory2 <- FindNeighbors(ECstory2, reduction = "integrated.harmony", dims = 1:30)
ECstory2 <- FindClusters(ECstory2, resolution = 0.5, cluster.name = "harmony_clusters")

ECstory2 <- RunUMAP(ECstory2, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

DimPlot(ECstory2, reduction = "umap.harmony")

FeaturePlot_scCustom(ECstory2, reduction="umap.harmony", features="PECAM1")

DimPlot(ECstory2, reduction = "umap.harmony", label=T)

#subset endos
PAendos <- subset(ECstory2, harmony_clusters %in% c("10", "14", "9"))
DimPlot(PAendos, reduction = "umap.harmony")

PAendos <- FindNeighbors(PAendos, reduction = "integrated.harmony", dims = 1:30)
PAendos <- FindClusters(PAendos, resolution = 0.15)

PAendos <- RunUMAP(PAendos, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony", seed.use=1)

DimPlot(PAendos, reduction = "umap.harmony")

FeaturePlot_scCustom(PAendos, features="SULF1", reduction="umap.harmony")
FeaturePlot_scCustom(PAendos, features="PTPRC", reduction="umap.harmony")

PAendos <- JoinLayers(PAendos)
DimPlot(PAendos, reduction="umap.harmony", group.by = "seurat_clusters")
PAendos <- subset(PAendos, seurat_clusters %!in% c("3"))#reduce it dow to investigated 3 clusters


DimPlot(PAendos, reduction="umap.harmony", cols=c("0" = "#0312ad",  
                                                  "1" = "#b20505",  
                                                  "2" = "#03ad46" ))  

PAendos[["umap.harmony"]]@cell.embeddings[,2] <-
  -PAendos[["umap.harmony"]]@cell.embeddings[,2]

DimPlot(PAendos, reduction="umap.harmony")

# ---------------------------
# Step 1: Convert Seurat → Monocle CDS
# ---------------------------
cds <- as.cell_data_set(PAendos)

# Add UMAP embeddings
reducedDims(cds)$UMAP <- Embeddings(PAendos, "umap.harmony")

# Add cluster labels
cds@clusters$UMAP$clusters <- PAendos$seurat_clusters

# Learn trajectory
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)

# Optional: set root cells (or skip to let Monocle decide)
# root_cells <- WhichCells(PAendos, idents = "2")
cds <- order_cells(cds, root_cells = NULL)

# Add pseudotime to Seurat object
PAendos$pseudotime_monocle3 <- pseudotime(cds)


plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,          # don't label clusters
  label_leaves = TRUE,                # label terminal states
  label_branch_points = TRUE,         # show branch points
  graph_label_size = 3
)+
  scale_color_gradient(
    low = "lightgrey",  # orange
    high = "#ff8000"  # grey
  )




# ---------------------------
# Step 2: Identify genes that vary along pseudotime
# ---------------------------
var_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

top_genes <- rownames(
  head(arrange(var_genes, q_value), 500)
)

GOI <- c("EDNRB", "HPGD", "PRKG1", "TBX2", "RCSD1", "B3GALNT1", "NCALD", "S100A4",#acap
         "IL7R", "CD36", "SLC6A4", "ARHGAP18", "IL18R1",#gcap
         "DKK2", "IGFBP3", "FBLN5", "GJA5", "CXCL12", "LTBP4", "HEY1", "SEMA3G",#arterial
         "CPE", "CLU", "C7", "PTGS1", "EFEMP1", "MMRN1", "DKK3", "PLAT",#pulmonary veinous
         "COL15A1", "TSHZ2", "FLRT2", "OLFM1", "PLVAP", "VWA1", "SPRY1",#systemic veinous
         "PRX", "F2RL3", "CA4", "HPGD", "AFF3", "IL1RL1", #stress response
         "PTGIS", "PDE3A", "DKK2", "NLGN1", "CP", "SERPINE2", #vascular tone
         "ACKR1", "ADAM23", "TM4SF18", "TESC", "MEOX1") #immuno 

# ---------------------------
# Step 3: Prepare expression matrix
# ---------------------------
expr_matrix <- GetAssayData(PAendos, layer = "data")  # normalized data

# Keep only top variable genes
selected_genes <- intersect(rownames(expr_matrix), GOI)#GOI or top_genes
expr_subset <- expr_matrix[selected_genes, , drop = FALSE]

# Filter out NA pseudotime cells
valid_cells <- !is.na(PAendos$pseudotime_monocle3)
expr_subset <- expr_subset[, valid_cells]
pseudotime_vals <- PAendos$pseudotime_monocle3[valid_cells]

# Order cells by pseudotime
ordered_cells <- names(sort(pseudotime_vals))
expr_ordered <- expr_subset[, ordered_cells]

# Z-score scale the gene expression (row-wise)
expr_scaled <- t(scale(t(as.matrix(expr_ordered))))

# ---------------------------
# Step 4: Cluster & Pseudotime Annotations
# ---------------------------
# Cluster bar (Seurat clusters reordered)
clusters <- PAendos$seurat_clusters[ordered_cells]
clusters <- factor(clusters, levels = c("0", "1", "2"))

# Color mappings
cluster_colors <- c(
  "0" = "#0312ad",  # red"#0312ad"
  "1" = "#b20505",  # blue
  "2" = "#03ad46"   # green
)

# Pseudotime bar (lightgrey → orange)
pseudotime_ordered <- pseudotime_vals[ordered_cells]
pseudo_col_fun <- colorRamp2(
  range(pseudotime_ordered, na.rm = TRUE),
  c("lightgrey", "#ff8000")
)

# Annotations
col_anno <- HeatmapAnnotation(
  Cluster = clusters,
  Pseudotime = anno_simple(pseudotime_ordered, col = pseudo_col_fun),
  col = list(Cluster = cluster_colors),
  annotation_legend_param = list(Cluster = list(title = "Cluster"))
)



gene_clusters <- kmeans(expr_scaled, centers = 4)$cluster

# Turn into a factor to allow manual reordering
gene_clusters <- factor(gene_clusters, levels = c(3, 2, 4, 1))  # REORDER HERE

# Optional: Rename cluster labels
names(gene_clusters) <- rownames(expr_scaled)


# ---------------------------
# Step 5: Draw Heatmap
# ---------------------------

Heatmap(
  expr_scaled,
  name = "Z-score",
  cluster_rows = FALSE,
  row_split = gene_clusters,  # split vertically by k-means clusters
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  top_annotation = col_anno,
  column_title = "Top 500 genes across Monocle 3 pseudotime",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = gpar(fontsize = 8),  # <-- change font size for gene labels here
  col = colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
)










