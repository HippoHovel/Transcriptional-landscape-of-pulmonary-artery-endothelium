#load CA and AO datasets
counts <- read.table(gzfile("~/big sets important/coronary artery/GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt.gz"), header = TRUE, sep = "\t", row.names = 1)

CoronaryArtery <- CreateSeuratObject (counts=counts)


#CoronaryArtery[["RNA"]] <- as(CoronaryArtery[["RNA"]], Class = "Assay5")

#set up CA data
CoronaryArtery <- NormalizeData(CoronaryArtery)
CoronaryArtery <- FindVariableFeatures(CoronaryArtery)
CoronaryArtery <- ScaleData(CoronaryArtery)
CoronaryArtery <- RunPCA(CoronaryArtery)
#CoronaryArtery <- JoinLayers(CoronaryArtery)
CoronaryArtery$compartment <- "CA"

cell_names <- colnames(CoronaryArtery)

#data does not supply meta data but cell names contain sample numbers --> extract sample numbers and split data set for integration
sample_ids <- sub(".*\\.(\\d+)$", "\\1", cell_names)
CoronaryArtery$sample_id <- sample_ids


CoronaryArtery[["RNA"]] <- split(CoronaryArtery[["RNA"]], f = CoronaryArtery$sample_id)  

CoronaryArtery <- IntegrateLayers(
  object = CoronaryArtery, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE)

CoronaryArtery <- FindNeighbors(CoronaryArtery, reduction = "integrated.cca", dims = 1:30)
CoronaryArtery <- FindClusters(CoronaryArtery, resolution = 0.5, cluster.name = "cca_clusters")

CoronaryArtery <- RunUMAP(CoronaryArtery, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

DimPlot(CoronaryArtery)
DimPlot(CoronaryArtery, group.by = "sample_id")
FeaturePlot_scCustom(CoronaryArtery, features="CDH5")
FeaturePlot(CoronaryArtery, features="PECAM1")
FeaturePlot(CoronaryArtery, features="ACKR1")
FeaturePlot(CoronaryArtery, features="IL1RL1")
FeaturePlot(CoronaryArtery, features="PTGIS")
FeaturePlot(CoronaryArtery, features="SERPINE2")


#set up AO data
con4 <- read.table(gzfile("~/big sets important/aorta/GSM4704931_Con4.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
con4 <- CreateSeuratObject (counts=con4)
con4$sample.id <- "AOcon4"

con6 <- read.table(gzfile("~/big sets important/aorta/GSM4704932_Con6.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
con6 <- CreateSeuratObject (counts=con6)
con6$sample.id <- "AOcon6"

con9 <- read.table(gzfile("~/big sets important/aorta/GSM4704933_Con9.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
con9 <- CreateSeuratObject (counts=con9)
con9$sample.id <- "AOcon9"

TAA1 <- read.table(gzfile("~/big sets important/aorta/GSM4704934_TAA1.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
TAA1 <- CreateSeuratObject (counts=TAA1)
TAA1$sample.id <- "AOTAA1"

TAA2 <- read.table(gzfile("~/big sets important/aorta/GSM4704935_TAA2.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
TAA2 <- CreateSeuratObject (counts=TAA2)
TAA2$sample.id <- "AOTAA2"

TAA3 <- read.table(gzfile("~/big sets important/aorta/GSM4704936_TAA3.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
TAA3 <- CreateSeuratObject (counts=TAA3)
TAA3$sample.id <- "AOTAA3"

TAA4 <- read.table(gzfile("~/big sets important/aorta/GSM4704937_TAA4.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
TAA4 <- CreateSeuratObject (counts=TAA4)
TAA4$sample.id <- "AOTAA4"

TAA5 <- read.table(gzfile("~/big sets important/aorta/GSM4704938_TAA5.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
TAA5 <- CreateSeuratObject (counts=TAA5)
TAA5$sample.id <- "AOTAA5"

TAA6 <- read.table(gzfile("~/big sets important/aorta/GSM4704939_TAA6.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
TAA6 <- CreateSeuratObject (counts=TAA6)
TAA6$sample.id <- "AOTAA6"

TAA7 <- read.table(gzfile("~/big sets important/aorta/GSM4704940_TAA7.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
TAA7 <- CreateSeuratObject (counts=TAA7)
TAA7$sample.id <- "AOTAA7"

TAA8 <- read.table(gzfile("~/big sets important/aorta/GSM4704941_TAA8.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
TAA8 <- CreateSeuratObject (counts=TAA8)
TAA8$sample.id <- "AOTAA8"


Aorta <- merge(con4, y = c(con6, con9, TAA1, TAA2, TAA3, TAA4, TAA5, TAA6, TAA7, TAA8), 
               add.cell.ids = c("con4","con6", "con9", "TAA1", "TAA2", "TAA3", "TAA4", "TAA5", "TAA6", "TAA7", "TAA8"), 
               project = "Aorta")

Aorta <- NormalizeData(Aorta)
Aorta <- FindVariableFeatures(Aorta)
Aorta <- ScaleData(Aorta)
Aorta <- RunPCA(Aorta)
#Aorta <- JoinLayers(Aorta)
Aorta$compartment <- "AO"

Aorta <- IntegrateLayers(
  object = Aorta, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE)

Aorta <- FindNeighbors(Aorta, reduction = "integrated.cca", dims = 1:30)
Aorta <- FindClusters(Aorta, resolution = 0.5, cluster.name = "cca_clusters")
Aorta <- RunUMAP(Aorta, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

DimPlot(Aorta)
FeaturePlot(Aorta, features="PECAM1")
FeaturePlot(Aorta, features="ACKR1")
FeaturePlot(Aorta, features="IL1RL1")
FeaturePlot(Aorta, features="PTGIS")
FeaturePlot(Aorta, features="CA4")

ECstory.j <- subset(ECstory.j, sample.id %in% ("PAH4"))
ECstory.j$compartment <- "PA"
Comapartment.comp <- merge(ECstory.j, y = c(Aorta, CoronaryArtery), 
                           add.cell.ids = c("PA","AO", "CA"), 
                           project = "Compartments")


Comapartment.comp <- NormalizeData(Comapartment.comp)
Comapartment.comp <- FindVariableFeatures(Comapartment.comp)
Comapartment.comp <- ScaleData(Comapartment.comp)
Comapartment.comp <- RunPCA(Comapartment.comp)

Comapartment.comp <- IntegrateLayers(
  object = Comapartment.comp, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE, normalization.method = "LogNormalize") #takes quite a while 

DimPlot(Comapartment.comp)


#QC
Comapartment.comp[["percent.mt"]] <- PercentageFeatureSet(Comapartment.comp, pattern = "^MT-")

Comapartment.comp $QCcut <- ""
Comapartment.comp $QCcut <- "Viable"
Comapartment.comp $QCcut[Comapartment.comp $nFeature_RNA > 4000] <- "Multiplett"
Comapartment.comp $QCcut[Comapartment.comp $nFeature_RNA < 200] <- "Cell remnant"
Comapartment.comp $QCcut[Comapartment.comp $percent.mt > 10] <- "low quality read"

# temp_dir <- tempdir()
# files <- list.files(temp_dir, full.names = T)
# files
# file.remove(files)

memory.size(max=T) # gives the amount of memory obtained by the OS

memory.size(max=F) # gives the amount of memory being used

#cut inviable cells
Comapartment.comp <- subset(Comapartment.comp, QCcut %in% c("Viable"))

Comapartment.comp <- FindNeighbors(Comapartment.comp, reduction = "integrated.cca", dims = 1:30)
Comapartment.comp <- FindClusters(Comapartment.comp, resolution = 0.4, cluster.name = "cca_clusters")
Comapartment.comp <- RunUMAP(Comapartment.comp, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")

#Identify ECs
DimPlot(Comapartment.comp, cols = polychrome_pal, label=T)
FeaturePlot(Comapartment.comp, features="PECAM1")
FeaturePlot(Comapartment.comp, features = "PECAM1",order =T)+
  scale_colour_gradient2(name = "expression" , high = "#b800ef", mid =  "#d8d8d8")
FeaturePlot(Comapartment.comp, features="ACKR1")
FeaturePlot(Comapartment.comp, features="PTGIS")
FeaturePlot(Comapartment.comp, features="IL1RL1")
FeaturePlot(Comapartment.comp, features="LYVE1")
FeaturePlot(Comapartment.comp, features="CA4")

Comapartment.comp$seurat_clusters <- Idents(Comapartment.comp)
DimPlot(Comapartment.comp, cols = polychrome_pal, label=T, group.by= "seurat_clusters")
FeaturePlot_scCustom(Comapartment.comp, features="PECAM1")
Compartment.comp.endo <- subset(Comapartment.comp, seurat_clusters %in% c("7" ,"11"))

DimPlot(Compartment.comp.endo, cols = polychrome_pal)

Compartment.comp.endo$QCcut <- ""
Compartment.comp.endo$QCcut <- "Viable"
Compartment.comp.endo$QCcut[Compartment.comp.endo$nFeature_RNA > 4000] <- "Multiplett"
Compartment.comp.endo$QCcut[Compartment.comp.endo$nFeature_RNA < 200] <- "Cell remnant"
Compartment.comp.endo$QCcut[Compartment.comp.endo$percent.mt > 10] <- "low quality read"

DimPlot(Compartment.comp.endo, group.by = "QCcut")
DimPlot(Compartment.comp.endo, group.by = "samplecondition", split.by = "QCcut")


Compartment.comp.endo <- subset(Compartment.comp.endo, QCcut %in% c("Viable"))

Compartment.comp.endo <- JoinLayers(Compartment.comp.endo)

DimPlot(Compartment.comp.endo)

Compartment.comp.endo <- FindNeighbors(Compartment.comp.endo, reduction = "integrated.cca", dims = 1:30)
Compartment.comp.endo <- FindClusters(Compartment.comp.endo, resolution = 0.2, cluster.name = "cca_clusters")
Compartment.comp.endo <- RunUMAP(Compartment.comp.endo, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

DimPlot(Compartment.comp.endo)
FeaturePlot(Compartment.comp.endo, features="PECAM1")
FeaturePlot(Compartment.comp.endo, features="ACKR1")
FeaturePlot(Compartment.comp.endo, features="PTGIS")
FeaturePlot(Compartment.comp.endo, features="IL1RL1")
FeaturePlot(Compartment.comp.endo, features="LYVE1")
FeaturePlot(Compartment.comp.endo, features="CA4")



Compartment.comp.endo.deg <- FindAllMarkers(Compartment.comp.endo)
Compartment.comp.endo.deg3 <- subset(Compartment.comp.endo.deg, cluster %in% c("3"))
FeaturePlot(Compartment.comp.endo, features="PLN")
FeaturePlot(Compartment.comp.endo, features="ITGA8")
FeaturePlot(Compartment.comp.endo, features="NEXN")
FeaturePlot(Compartment.comp.endo, features="RCAN2")
FeaturePlot(Compartment.comp.endo, features="SYNPO2")
FeaturePlot(Compartment.comp.endo, features="AKAP6")
FeaturePlot(Compartment.comp.endo, features="LTBP1")
FeaturePlot(Compartment.comp.endo, features="ACTA2")

DimPlot(Compartment.comp.endo, split.by = "compartment", order = T)
DimPlot(Compartment.comp.endo, split.by = "samplecondition")
#SMC markers across the board
#endoMT? if so why are they again mostly found in Donor samples
#sort out



Compartment.comp.endo.deg4 <- subset(Compartment.comp.endo.deg, cluster %in% c("4"))
FeaturePlot(Compartment.comp.endo, features="SYTL3")
FeaturePlot(Compartment.comp.endo, features="CD3D")
FeaturePlot(Compartment.comp.endo, features="CD3E")
FeaturePlot(Compartment.comp.endo, features="CD2")
FeaturePlot(Compartment.comp.endo, features="IKZF1")
FeaturePlot(Compartment.comp.endo, features="BCL11B")
FeaturePlot(Compartment.comp.endo, features="PTPRC")
#immune cell markers
#sort out

Compartment.comp.endo$seurat_clusters <- Idents(Compartment.comp.endo)
Compartment.comp.endo.fil <- subset(Compartment.comp.endo, seurat_clusters %!in% c("3", "4"))
DimPlot(Compartment.comp.endo.fil, cols= c("#b20505","#0312ad","#03ad46"))
DimPlot(Compartment.comp.endo.fil, cols= c("#b20505","#0312ad","#03ad46"), split.by = "compartment")
Compartment.comp.endo.fil <- RunUMAP(Compartment.comp.endo.fil, reduction = "integrated.cca", dims = 1:30 , reduction.name = "umap.cca")
DimPlot(Compartment.comp.endo.fil, cols= c("#b20505","#0312ad","#03ad46"), split.by = "compartment")








DimPlot(Compartment.comp.endo.fil, cols = polychrome_pal, split.by = "compartment")
DimPlot(Compartment.comp.endo.fil, cols = c("#b20505","#0312ad","#03ad46"), split.by = "compartment")


FeaturePlot(Compartment.comp.endo.fil, features = "CA4")

FeaturePlot(Compartment.comp.endo.fil, features = "ACKR1", cols = "#b20505")
FeaturePlot_scCustom(Compartment.comp.endo.fil, features="ACKR1", colors_use = c("#e4b04e","#b20505"), split.by = "compartment")
FeaturePlot_scCustom(Compartment.comp.endo.fil, features="PTGIS", colors_use = c("#71cfd3","#142b9e"), split.by = "compartment")
FeaturePlot_scCustom(Compartment.comp.endo.fil, features="IL1RL1", colors_use = c("#98f12f","#03ad46"), split.by = "compartment")
FeaturePlot_scCustom(Compartment.comp.endo.fil, features="PRX", colors_use = c("#98f12f","#03ad46"), split.by = "compartment")

Compartment.comp.endo.fil$clusternames <- ""
Compartment.comp.endo.fil$clusternames[Compartment.comp.endo.fil$seurat_clusters == "0"] <- "Immunomodulatory"
Compartment.comp.endo.fil$clusternames[Compartment.comp.endo.fil$seurat_clusters == "1"] <- "Vascular tone modulating"
Compartment.comp.endo.fil$clusternames[Compartment.comp.endo.fil$seurat_clusters == "2"] <- "Vascular plasticity modulating"

Compartment.comp.endo.fil$clusternames <- factor(Compartment.comp.endo.fil$clusternames, levels = c("Immunomodulatory", "Vascular tone modulating", "Vascular plasticity modulating" ))

# DimPlot(Compartment.comp, split.by = "compartment", cols=polychrome_pal)
# FeaturePlot(Compartment.comp, features="LYVE1")
# DimPlot(Compartment.comp, split.by = "compartment", group.by = "samplecondition")

DimPlot(Compartment.comp.endo.fil, group.by = "clusternames", cols=c("#b20505","#0312ad","#03ad46"))


cluster0GOIhum <- c("TSHZ2", "ACKR1", "PLVAP","SPRY1", "ADAMTS9", "SNTG2", "POSTN", "COL15A1", "EBF1",
                    "FABP5", "NOSTRIN", "FABP4", "MEOX1", "OLFM1", "LRRC1", "CYTL1", "MEOX2",
                    "NEURL1B", "TPD52L1", "STOX2", "APOLD1", "ADAM23", "RBP7", "ANKRD29",
                    "PDE2A", "ELOVL7", "PLAT", "C7","NRP2","FLRT2")

Compartment.comp.endo.fil <- AddModuleScore(Compartment.comp.endo.fil, features=list(cluster0GOIhum), name="score0")
names(x = Compartment.comp.endo.fil[[]])

# FeaturePlot_scCustom(Compartment.comp.endo.fil, features="score01")
# VlnPlot_scCustom(Compartment.comp.endo.fil, features="score01", group.by="clusternames", colors_use = c("#0312ad","#b20505","#03ad46"), pt.size=0 )
# VlnPlot_scCustom(Compartment.comp.endo.fil, features="score21", group.by="seurat_clusters", colors_use = c("#5c109e","#ca8310","#0eb4b4"), pt.size=0 )
VlnPlot_scCustom(Compartment.comp.endo.fil, features="score01", group.by="seurat_clusters", colors_use = c("#b20505","#0312ad","#03ad46"), pt.size=0 )




cluster1GOIhum <- c("PDE3A", "DKK2", "IGFBP3", "PTGIS", "FN1", "SULF1", "BMX", "ABI3BP", "HMCN1", "SOX5",
                    "SETBP1", "BGN", "EDN1", "GJA5", "CP", "SERPINE2", "CPE", "COL8A1", "OMD",
                    "SERPINE1", "EFNA5", "EFEMP1", "CDH11", "TMOD1", "CELF2", "ST8SIA6", "ELN", "FBLIM1",
                    "PAPSS2", "CSMD1")


Compartment.comp.endo.fil <- AddModuleScore(Compartment.comp.endo.fil, features=list(cluster1GOIhum), name="score1")
names(x = Compartment.comp.endo.fil[[]])

FeaturePlot_scCustom(Compartment.comp.endo.fil, features="score11")
VlnPlot_scCustom(Compartment.comp.endo.fil, features="score11", group.by="seurat_clusters", colors_use = c("#b20505","#0312ad","#03ad46"), pt.size=0 )



cluster2GOIhum <- c("AFF3", "STXBP6", "BTNL9", "RGCC", "HECW2", "ADGRF5", "HPGD", "CAR4", "PRICKLE2",
                    "IL1RL1", "TMEM100", "PRX", "FENDRR", "IL7R","SLC6A4", "VIPR1", "NCALD", "EDNRB",
                    "ARHGAP6", "TBX3", "ARHGAP18", "PREX1", "PRKG1", "GMFG", "CX3CL1", "ISG20",
                    "CD109", "ROBO4", "ARHGEF10", "FOXF1")

Compartment.comp.endo.fil <- AddModuleScore(Compartment.comp.endo.fil, features=list(cluster2GOIhum), name="score2")
names(x = Compartment.comp.endo.fil[[]])

FeaturePlot_scCustom(Compartment.comp.endo.fil, features="score21")
VlnPlot_scCustom(Compartment.comp.endo.fil, features="score21", group.by="seurat_clusters", colors_use = c("#b20505","#0312ad","#03ad46"), pt.size=0 )



Compartment.comp.endo.fil$clusternames <- ""
Compartment.comp.endo.fil$clusternames[Compartment.comp.endo.fil$seurat_clusters == "0"] <- "Immunomodulatory"
Compartment.comp.endo.fil$clusternames[Compartment.comp.endo.fil$seurat_clusters == "1"] <- "Vascular tone modulating"
Compartment.comp.endo.fil$clusternames[Compartment.comp.endo.fil$seurat_clusters == "2"] <- "Vascular plasticity modulating"


Compartment.comp.endo.fil$clusternames <- factor(Compartment.comp.endo.fil$clusternames, levels=c("Immunomodulatory", "Vascular tone modulating", "Vascular plasticity modulating"))
DimPlot(Compartment.comp.endo.fil, group.by = "clusternames", cols =c("#b20505","#0312ad","#03ad46"), split.by = "compartment")

#Piechart
sum((Compartment.comp.endo.fil$clusternames == "Vascular tone modulating")) 
sum((Compartment.comp.endo.fil$clusternames == "Immunomodulatory")) 
sum((Compartment.comp.endo.fil$clusternames == "Vascular plasticity modulating"))

Compartment.comp.endo.AO <- subset(Compartment.comp.endo.fil, compartment %in% c("AO"))
DimPlot(Compartment.comp.endo.AO)

Compartment.comp.endo.CA <- subset(Compartment.comp.endo.fil, compartment %in% c("CA"))
DimPlot(Compartment.comp.endo.CA)

Compartment.comp.endo.PA <- subset(Compartment.comp.endo.fil, compartment %in% c("PA"))
DimPlot(Compartment.comp.endo.PA)


Compartment.comp.endo.df <- data.frame( # change subset based on the piechart to be plotted Compartment.comp.endo.PA/CA/AO
  group=c("Vascular tone modulating", "Immunomodulatory", "Vascular plasticity modulating"),
  value=c(sum((Compartment.comp.endo.PA$clusternames == "Vascular tone modulating")) ,
          sum((Compartment.comp.endo.PA$clusternames== "Immunomodulatory")),
          sum((Compartment.comp.endo.PA$clusternames == "Vascular plasticity modulating")))
)

#old piecharts
# Compute the position of labels
Compartment.comp.endo.df <- Compartment.comp.endo.df %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(Compartment.comp.endo.df$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) 
mutate(stdev = sd(value))


# Basic piechart
ggplot(Compartment.comp.endo.df, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=c("#03ad46","#b20505","#0312ad"))


#PA-----------------------------
meta_df_PA <- Compartment.comp.endo.PA@meta.data

#Count cells per clusternames × samplecondition × sample.ident
cell_counts_PA <- meta_df_PA %>%
  count(clusternames, samplecondition, sample.ident, name = "cell_count")

#Calculate percentage per sample.ident
cell_counts_PA <- cell_counts_PA %>%
  group_by(sample.ident) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

# View results
head(cell_counts_PA)


summary_stats_PA <- cell_counts_PA %>%
  group_by(clusternames) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )


#AO
meta_df_AO <- Compartment.comp.endo.AO@meta.data

#Count cells per clusternames × samplecondition × sample.ident
cell_counts_AO <- meta_df_AO %>%
  count(clusternames, samplecondition, sample.id, name = "cell_count")

#Calculate percentage per sample.ident
cell_counts_AO <- cell_counts_AO %>%
  group_by(sample.id) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

# View results
head(cell_counts_AO)


summary_stats_AO <- cell_counts_AO %>%
  group_by(clusternames) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )

summary_stats_AO$test <- ""
summary_stats_AO$test[summary_stats_AO$clusternames == "Vascular tone modulating"] <- "5"
summary_stats_AO$test[summary_stats_AO$clusternames == "Vascular plasticity modulating"] <- "50"
summary_stats_AO$test[summary_stats_AO$clusternames == "Immunomodulatory"] <- "100"


#CA
meta_df_CA <- Compartment.comp.endo.CA@meta.data

# Step 1: Count cells per clusternames × samplecondition × sample.ident
cell_counts_CA <- meta_df_CA %>%
  count(clusternames, samplecondition, sample_id, name = "cell_count")

# Step 2: Calculate percentage per sample.ident
cell_counts_CA <- cell_counts_CA %>%
  group_by(sample.id) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

# View results
print(cell_counts_CA, n = 35)


summary_stats_CA <- cell_counts_CA %>%
  group_by(clusternames) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )

print(summary_stats_CA)


#CA
meta_df_CA<- Compartment.comp.endo.CA@meta.data

#Count cells per clusternames × samplecondition × sample.ident
cell_counts_CA <- meta_df_CA %>%
  count(clusternames, samplecondition, sample_id, name = "cell_count")

#Calculate percentage per sample.ident
cell_counts_CA <- cell_counts_CA %>%
  group_by(sample_id) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

# View results
head(cell_counts_CA)


summary_stats_CA<- cell_counts_CA %>%
  group_by(clusternames) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )

summary_stats_CA$test <- ""
summary_stats_CA$test[summary_stats_AO$clusternames == "Vascular tone modulating"] <- "5"
summary_stats_CA$test[summary_stats_AO$clusternames == "Vascular plasticity modulating"] <- "50"
summary_stats_CA$test[summary_stats_AO$clusternames == "Immunomodulatory"] <- "100"



geom_text(aes(y = ypos, label = group), color = "white", size=6)
# geom_text(aes(y = ypos,x = xpos ,label=format(prop,digits=2)), color = "white", size=6)+
#scale_fill_manual(values=c("#0312ad","#b20505","#03ad46"))


#newpiecharts
ggplot(summary_stats_PA, aes(x="", y=median_percent, fill=cluster_names)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=c("#b20505","#0312ad","#03ad46"))

ggplot(summary_stats_AO, aes(x="", y=median_percent, fill=cluster_names)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=c("#b20505","#0312ad","#03ad46"))

ggplot(summary_stats_CA, aes(x="", y=median_percent, fill=cluster_names)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=c("#b20505","#0312ad","#03ad46"))



ggplot(cell_counts_PA, aes(x = clusternames, y = percent)) +
  geom_boxplot(fill = c("#b20505","#0312ad","#03ad46"), color = "black", outlier.shape = 21, outlier.fill = "red") +
  labs(
    title = "Distribution of Cluster Percentages Across Samples",
    x = "Cluster Name",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(cell_counts_AO, aes(x = clusternames, y = percent)) +
  geom_boxplot(fill = c("#b20505","#0312ad","#03ad46"), color = "black", outlier.shape = 21, outlier.fill = "red") +
  labs(
    title = "Distribution of Cluster Percentages Across Samples",
    x = "Cluster Name",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot(cell_counts_CA, aes(x = clusternames, y = percent)) +
  geom_boxplot(fill = c("#b20505","#0312ad","#03ad46"), color = "black", outlier.shape = 21, outlier.fill = "red") +
  labs(
    title = "Distribution of Cluster Percentages Across Samples",
    x = "Cluster Name",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot(cell_counts_CA, aes(x = clusternames, y = percent)) +
  geom_boxplot(
    fill = c("#b20505","#03ad46","#0312ad"),
    color = "black",
    outlier.shape = 21,
    outlier.fill = "red"
  ) +
  geom_jitter(
    aes(color = clusternames),
    width = 0.15,
    size = 2.5,
    alpha = 0.9
  ) +
  # geom_text(
  #   data = stats_PHPF,
  #   aes(
  #     x = clusternames,
  #     y = label_y,
  #     label = paste0(
  #       "Median:\n", round(median, 1), "%\n\n",
  #       "IQR:\n", round(q1, 1), "–", round(q3, 1), "%"
  #     )
  #   ),
  #   size = 2.5,
  #   fontface = "bold"
  # ) +
  scale_fill_manual(values = c("#CC4B4B","#57D487","#454FD1")) +
  scale_color_manual(values = c("#CC4B4B","#57D487","#454FD1")) +
  scale_y_continuous(limits = c(0, 100))+
  labs(
    title = "Cluster Proportions Across PHPF Samples",
    x = "Cluster",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  expand_limits(y = max(stats_PHPF$label_y) + 2)



ggplot(cell_counts_PA, aes(x = clusternames, y = percent)) +
  geom_boxplot(
    fill = c("#b20505","#03ad46","#0312ad"),
    color = "black",
    outlier.shape = 21#,
    # outlier.fill = "red"
  ) +
  geom_jitter(
    aes(color = clusternames),
    width = 0.15,
    size = 2.5,
    alpha = 0.9
  ) +
  # geom_text(
  #   data = stats_PHPF,
  #   aes(
  #     x = clusternames,
  #     y = label_y,
  #     label = paste0(
  #       "Median:\n", round(median, 1), "%\n\n",
  #       "IQR:\n", round(q1, 1), "–", round(q3, 1), "%"
  #     )
  #   ),
  #   size = 2.5,
  #   fontface = "bold"
  # ) +
  scale_fill_manual(values = c("#CC4B4B","#57D487","#454FD1")) +
  scale_color_manual(values = c("#CC4B4B","#57D487","#454FD1")) +
  scale_y_continuous(limits = c(0, 100))+
  labs(
    title = "Cluster Proportions Across PHPF Samples",
    x = "Cluster",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  expand_limits(y = max(stats_PHPF$label_y) + 2)



ggplot(cell_counts_AO, aes(x = clusternames, y = percent)) +
  geom_boxplot(
    fill = c("#b20505","#03ad46","#0312ad")#,
    #color = "black",
    #outlier.shape = 21#,
    # outlier.fill = "red"
  ) +
  geom_jitter(
    aes(color = clusternames),
    width = 0.15,
    size = 2.5,
    alpha = 0.9
  ) +
  # geom_text(
  #   data = stats_PHPF,
  #   aes(
  #     x = clusternames,
  #     y = label_y,
  #     label = paste0(
  #       "Median:\n", round(median, 1), "%\n\n",
  #       "IQR:\n", round(q1, 1), "–", round(q3, 1), "%"
  #     )
  #   ),
  #   size = 2.5,
  #   fontface = "bold"
  # ) +
  scale_fill_manual(values = c("#CC4B4B","#57D487","#454FD1")) +
  scale_color_manual(values = c("#CC4B4B","#57D487","#454FD1")) +
  scale_y_continuous(limits = c(0, 100))+
  labs(
    title = "Cluster Proportions Across PHPF Samples",
    x = "Cluster",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  expand_limits(y = max(stats_PHPF$label_y) + 2)


library(writexl)

cell_counts_AO$vascular_bed <- "AO"
cell_counts_CA$vascular_bed <- "CA"
cell_counts_PA$vascular_bed <- "PA"

combined_raw_vascular <- bind_rows(
  cell_counts_AO,
  cell_counts_CA,
  cell_counts_PA
)

combined_summary_vascular <- combined_raw_vascular %>%
  group_by(vascular_bed, clusternames) %>%
  summarise(
    median = median(percent),
    Q1 = quantile(percent,0.25),
    Q3 = quantile(percent,0.75),
    IQR = IQR(percent),
    mean = mean(percent),
    sd = sd(percent),
    n = n(),
    .groups="drop"
  )

combined_full_vascular <- left_join(
  combined_raw_vascular,
  combined_summary_vascular,
  by=c("vascular_bed","clusternames")
)

write_xlsx(combined_full_vascular,
           "EC_cluster_proportions_vascularbeds_SINGLE_TABLE.xlsx")

#mouse hoxnox

HoxNoxintegrated <- readRDS("HoxNoxintegratedtommy.rds")

DimPlot(HoxNoxintegrated, group.by = "Cell_annotation")
FeaturePlot_scCustom(HoxNoxintegrated, features="Prox1")

hoxnox_endo <- subset(HoxNoxintegrated, Cell_annotation %in% c("Endo 1", "Endo 2"))
DefaultAssay(hoxnox_endo) <- "RNA"
hoxnox_endo <- JoinLayers(hoxnox_endo)

hoxnox_endo <- ScaleData(hoxnox_endo, verbose = T)
hoxnox_endo <- RunPCA(hoxnox_endo, verbose = T)
hoxnox_endo <- RunUMAP(hoxnox_endo, dims = 1:30)

hoxnox_endo <- FindNeighbors(hoxnox_endo, reduction = "pca", dims = 1:30)
hoxnox_endo <- FindClusters(hoxnox_endo, resolution = 3.5)

DimPlot(hoxnox_endo, group.by = "seurat_clusters", cols = polychrome_pal)
hoxnox_endo$seurat_clusters_manual <- ""
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "0"] <- "V"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "1"] <- "II"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "2"] <- "III"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "3"] <- "III"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "4"] <- "III"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "5"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "6"] <- "V"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "7"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "8"] <- "I"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "9"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "10"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "11"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "12"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "13"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "14"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "15"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "16"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "17"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "18"] <- "IV"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "19"] <- "III"
hoxnox_endo$seurat_clusters_manual[hoxnox_endo$seurat_clusters == "20"] <- "IV"



DimPlot(hoxnox_endo, group.by = "seurat_clusters_manual",cols=c("#5c109e","#ca8310","#0eb4b4","#d535d7","#96c61e"))
DimPlot(hoxnox_endo, group.by = "seurat_clusters_manual",cols=c("#b20505","#0312ad","#025d11","#626fe8","#1acb5f"))

FeaturePlot_scCustom(hoxnox_endo, features= "Ackr1")
FeaturePlot_scCustom(hoxnox_endo, features= "Meox1")
FeaturePlot_scCustom(hoxnox_endo, features= "Plvap")
FeaturePlot_scCustom(hoxnox_endo, features= "Spry1")
FeaturePlot_scCustom(hoxnox_endo, features= "Adamts9")
FeaturePlot_scCustom(hoxnox_endo, features= "Sntg2")
FeaturePlot_scCustom(hoxnox_endo, features= "Postn")
FeaturePlot_scCustom(hoxnox_endo, features= "Col15a1")

cluster0GOI <- c("Tshz2", "Ackr1", "Plvap","Spry1", "Adamts9", "Sntg2", "Postn", "Col15a1", "Ebf1",
                 "Fabp5", "Nostrin", "Fabp4", "Meox1", "Olfm1", "Lrrc1", "Cytl1", "Meox2",
                 "Neurl1b", "Tpd52l1", "Stox2", "Apold1", "Adam23", "Rbp7", "Ankrd29",
                 "Pde2a", "Elovl7", "Plat", "C7","Nrp2","Flrt2" )
hoxnox_endo <- AddModuleScore(hoxnox_endo, features=list(cluster0GOI), name="score0")
names(x = hoxnox_endo[[]])

FeaturePlot_scCustom(hoxnox_endo, features="score01")
FeaturePlot_scCustom(hoxnox_endo, features="score01", split.by = "condition")
VlnPlot_scCustom(hoxnox_endo, features="score01", group.by="seurat_clusters_manual", colors_use = c("#5c109e","#ca8310","#0eb4b4","#d535d7","#96c61e"), pt.size=0.5)
VlnPlot_scCustom(hoxnox_endo, features="score01", group.by="seurat_clusters_manual", colors_use = c("#b20505","#0312ad","#025d11","#626fe8","#1acb5f"), pt.size=0)


cluster1GOI <- c("Pde3a", "Dkk2", "Igfbp3", "Ptgis", "Fn1", "Sulf1", "Bmx", "Abi3bp", "Hmcn1", "Sox5",
                 "Setbp1", "Bgn", "Edn1", "Gja5", "Cp", "Serpine2", "Cpe", "Col8a1", "Omd",
                 "Serpine1", "Efna5", "Efemp1", "Cdh11", "Tmod1", "Celf2", "St8sia6", "Eln", "Fblim1",
                 "Papss2", "Csmd1")
hoxnox_endo <- AddModuleScore(hoxnox_endo, features=list(cluster1GOI), name="score1")
names(x = hoxnox_endo[[]])

FeaturePlot_scCustom(hoxnox_endo, features="score11")
VlnPlot_scCustom(hoxnox_endo, features="score11", group.by="seurat_clusters_manual", colors_use = c("#5c109e","#ca8310","#0eb4b4","#d535d7","#96c61e"), pt.size=0.5)
VlnPlot_scCustom(hoxnox_endo, features="score11", group.by="seurat_clusters_manual", colors_use = c("#b20505","#0312ad","#025d11","#626fe8","#1acb5f"), pt.size=0)

cluster2GOI <- c("Aff3", "Stxbp6", "Btnl9", "Rgcc", "Hecw2", "Adgrf5", "Hpgd", "Car4", "Prickle2",
                 "Il1rl1", "Tmem100", "Prx", "Fendrr", "Il7r","Slc6a4", "Vipr1", "Ncald", "Ednrb",
                 "Arhgap6", "Tbx3", "Arhgap18", "Prex1", "Prkg1", "Gmfg", "Cx3cl1", "Isg20",
                 "Cd109", "Robo4", "Arhgef10", "Foxf1")
hoxnox_endo <- AddModuleScore(hoxnox_endo, features=list(cluster2GOI), name="score2")
names(x = hoxnox_endo[[]])

FeaturePlot_scCustom(hoxnox_endo, features="score21")
VlnPlot_scCustom(hoxnox_endo, features="score21", group.by="seurat_clusters_manual", colors_use = c("#5c109e","#ca8310","#0eb4b4","#d535d7","#96c61e"), pt.size=0.5)
VlnPlot_scCustom(hoxnox_endo, features="score21", group.by="seurat_clusters_manual", colors_use = c("#b20505","#0312ad","#025d11","#626fe8","#1acb5f"), pt.size=0)


hoxnox_endo$ec.clusters <- ""
hoxnox_endo$ec.clusters[hoxnox_endo$seurat_clusters_manual == "V"] <- "Defense"
hoxnox_endo$ec.clusters[hoxnox_endo$seurat_clusters_manual == "I"] <- "Immunomodulatory"
hoxnox_endo$ec.clusters[hoxnox_endo$seurat_clusters_manual == "III"] <- "Defense"
hoxnox_endo$ec.clusters[hoxnox_endo$seurat_clusters_manual == "II"] <- "Homeostatic"
hoxnox_endo$ec.clusters[hoxnox_endo$seurat_clusters_manual == "IV"] <- "Homeostatic"

DimPlot(hoxnox_endo, group.by = "ec.clusters", cols= c("#03ad46","#0312ad","#b20505"))


hoxnox_endo_nox <- subset(hoxnox_endo, condition %in% c("Nox"))

DimPlot(hoxnox_endo_nox, group.by = "seurat_clusters_manual")
hoxnox_endo_nox$seurat_clusters_manual[hoxnox_endo_nox$seurat_clusters_manual == "I"] <- "I"
hoxnox_endo_nox$seurat_clusters_manual[hoxnox_endo_nox$seurat_clusters_manual == "IV"] <- "II"
hoxnox_endo_nox$seurat_clusters_manual[hoxnox_endo_nox$seurat_clusters_manual == "V"] <- "III"
DimPlot(hoxnox_endo_nox, group.by = "seurat_clusters_manual", cols = c("#b20505","#0312ad","#03ad46"))

DimPlot(hoxnox_endo_nox, group.by = "ec.clusters")
DimPlot(hoxnox_endo_nox, group.by = "seurat_clusters_manual")

VlnPlot_scCustom(hoxnox_endo_nox, features="score01", group.by="seurat_clusters_manual", colors_use = c("#b20505","#0312ad","#025d11","#626fe8","#1acb5f"), pt.size=0)

cluster0GOI <- c("Tshz2", "Ackr1", "Plvap","Spry1", "Adamts9", "Sntg2", "Postn", "Col15a1", "Ebf1",
                 "Fabp5", "Nostrin", "Fabp4", "Meox1", "Olfm1", "Lrrc1", "Cytl1", "Meox2",
                 "Neurl1b", "Tpd52l1", "Stox2", "Apold1", "Adam23", "Rbp7", "Ankrd29",
                 "Pde2a", "Elovl7", "Plat", "C7","Nrp2","Flrt2" )
hoxnox_endo_nox <- AddModuleScore(hoxnox_endo_nox, features=list(cluster0GOI), name="score0")
names(x = hoxnox_endo_nox[[]])

FeaturePlot_scCustom(hoxnox_endo_nox, features="score01")
FeaturePlot_scCustom(hoxnox_endo_nox, features="score01", split.by = "condition")
VlnPlot_scCustom(hoxnox_endo_nox, features="score01", group.by="seurat_clusters_manual", colors_use = c("#5c109e","#ca8310","#0eb4b4","#d535d7","#96c61e"), pt.size=0.5)
VlnPlot_scCustom(hoxnox_endo_nox, features="score01", group.by="seurat_clusters_manual", colors_use = c("#b20505","#0312ad","#1acb5f"), pt.size=0)


cluster1GOI <- c("Pde3a", "Dkk2", "Igfbp3", "Ptgis", "Fn1", "Sulf1", "Bmx", "Abi3bp", "Hmcn1", "Sox5",
                 "Setbp1", "Bgn", "Edn1", "Gja5", "Cp", "Serpine2", "Cpe", "Col8a1", "Omd",
                 "Serpine1", "Efna5", "Efemp1", "Cdh11", "Tmod1", "Celf2", "St8sia6", "Eln", "Fblim1",
                 "Papss2", "Csmd1")
hoxnox_endo_nox <- AddModuleScore(hoxnox_endo_nox, features=list(cluster1GOI), name="score1")
names(x = hoxnox_endo_nox[[]])

FeaturePlot_scCustom(hoxnox_endo_nox, features="score11")
VlnPlot_scCustom(hoxnox_endo_nox, features="score11", group.by="seurat_clusters_manual", colors_use = c("#5c109e","#ca8310","#0eb4b4","#d535d7","#96c61e"), pt.size=0.5)
VlnPlot_scCustom(hoxnox_endo_nox, features="score11", group.by="seurat_clusters_manual", colors_use = c("#b20505","#0312ad","#1acb5f"), pt.size=0)

cluster2GOI <- c("Aff3", "Stxbp6", "Btnl9", "Rgcc", "Hecw2", "Adgrf5", "Hpgd", "Car4", "Prickle2",
                 "Il1rl1", "Tmem100", "Prx", "Fendrr", "Il7r","Slc6a4", "Vipr1", "Ncald", "Ednrb",
                 "Arhgap6", "Tbx3", "Arhgap18", "Prex1", "Prkg1", "Gmfg", "Cx3cl1", "Isg20",
                 "Cd109", "Robo4", "Arhgef10", "Foxf1")
hoxnox_endo_nox <- AddModuleScore(hoxnox_endo_nox, features=list(cluster2GOI), name="score2")
names(x = hoxnox_endo_nox[[]])

FeaturePlot_scCustom(hoxnox_endo_nox, features="score21")
VlnPlot_scCustom(hoxnox_endo_nox, features="score21", group.by="seurat_clusters_manual", colors_use = c("#5c109e","#ca8310","#0eb4b4","#d535d7","#96c61e"), pt.size=0.5)
VlnPlot_scCustom(hoxnox_endo_nox, features="score21", group.by="seurat_clusters_manual", colors_use = c("#b20505","#0312ad","#1acb5f"), pt.size=0)


#check hox vs nox deg
DimPlot(hoxnox_endo$, group.by = "condition")

Idents(hoxnox_endo) <- hoxnox_endo$condition
hoxnox.deg <- FindMarkers(hoxnox_endo, ident.1= "Hox", ident.2 = "Nox", test.use= "MAST")

top <- hoxnox.deg 

top <- rownames_to_column(top, var="gene")

top$avg_log2FC <- as.numeric(top$avg_log2FC)



#check the top 20 cluster defining genes
top20 <- top
top20$valid = ""
top20$valid [(top20$p_val_adj < 0.05)] <- TRUE #check for significance
top20$valid[(top20$pct.1 < 0.2)] <- FALSE #pick genes that the cluster is expressing to a high %
top20$valid[(top20$pct.2 > 0.2)] <- FALSE #dismiss genes that the other populations are also expressing to a high %
top20 <- subset(top20, subset = top20$valid == TRUE )

top20 <- top20[order(top20$cluster),]
top20 <- na.omit(top20)

DotPlot_scCustom(hoxnox_endo, features = c(rev(top20$gene)), flip_axes = T, x_lab_rotate = T, colors_use=c("#4f4ceb", "#807ee9", "#a8a7e6","#d8d8d8","#f38d8d","#f25656","#da0000"))+scale_size(range = c(1, 5))

hoxnox_endo <- ScaleData(hoxnox_endo, features = top20$gene)

plot_heatmap(dataset=hoxnox_endo,
             markers = top20$gene,
             sort_var = c("condition","seurat_clusters_manual"),
             anno_var = c("condition","seurat_clusters_manual"),
             anno_colors= list(polychrome_pal, c("#b20505","#0312ad","#025d11","#626fe8","#1acb5f")),
             row_font_size = 7.5)

top20$nr <- rownames(top20)

top20.2 <- top_n(top20, 100, avg_log2FC)
top20.2 <- top20.2[order(top20.2$avg_log2FC, decreasing=T ),]


top_n(x, n, wt)
top20.2 <- top20 %>% top_n(50)

top20.3 <- subset(top20, avg_log2FC <= 0)
top20.4 <- subset(top20, avg_log2FC >= 0)

top20.5 <- top_n(top20.3, 50, avg_log2FC)



plot_heatmap(dataset=hoxnox_endo,
             markers = top20.5$gene,
             sort_var = c("condition","seurat_clusters_manual"),
             anno_var = c("condition","seurat_clusters_manual"),
             anno_colors= list(polychrome_pal, c("#b20505","#0312ad","#025d11","#626fe8","#1acb5f")),
             row_font_size = 7.5)



DoHeatmap(hoxnox_endo, features=top20.4$gene, assay="RNA",slot="count")+scale_fill_gradientn(colors = c( "lightgrey", "#ff3829"))


DoHeatmap(hoxnox_endo, features=top20.4$gene, assay="RNA",slot="scale.data")+scale_fill_gradientn(colors = c( "lightgrey", "#ff3829"))


