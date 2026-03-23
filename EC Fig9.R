library(CellChat)




DimPlot(ECstory, cols=polychrome_pal)

FeaturePlot_scCustom(ECstory, features = c("PTPRC"))
FeaturePlot_scCustom(ECstory, features = c("GP1BB"))
FeaturePlot_scCustom(ECstory, features = c("ITGAM"))
FeaturePlot_scCustom(ECstory, features = c("ACTA2"))
FeaturePlot_scCustom(ECstory, features = c("VIM"))
FeaturePlot_scCustom(ECstory, features = c("DCN"))
FeaturePlot_scCustom(ECstory, features = c("ACTG2"))
FeaturePlot_scCustom(ECstory, features = c("VCAN"))
FeaturePlot_scCustom(ECstory, features = c("SOCS3"))
FeaturePlot_scCustom(ECstory, features = c("S100A10"))
FeaturePlot_scCustom(ECstory, features = c("LUM"))
FeaturePlot_scCustom(ECstory, features = c("TIMP1"))
FeaturePlot_scCustom(ECstory, features = c("FBLN1"))
FeaturePlot_scCustom(ECstory, features = c("APOE"))
FeaturePlot_scCustom(ECstory, features = c("RAMP1"))
FeaturePlot_scCustom(ECstory, features = c("APOD"))
FeaturePlot_scCustom(ECstory, features = c("GATA3"))
FeaturePlot_scCustom(ECstory, features = c("FOS"))
FeaturePlot_scCustom(ECstory, features = c("EIF2AK3"))
FeaturePlot_scCustom(ECstory, features = c("ABCC9"))
FeaturePlot_scCustom(ECstory, features = c("TFRC"))
FeaturePlot_scCustom(ECstory, features = c("CD36"))
FeaturePlot_scCustom(ECstory, features = c("GYPA"))
FeaturePlot_scCustom(ECstory, features = c("SELP"))
FeaturePlot_scCustom(ECstory, features = c("CALB2"))
FeaturePlot_scCustom(ECstory, features = c("KRT6A"))
FeaturePlot_scCustom(ECstory, features = c("KIT"))
FeaturePlot_scCustom(ECstory, features = c("KITLG"))
FeaturePlot_scCustom(ECstory, features = c("HBB"))
FeaturePlot_scCustom(ECstory, features = c("PTPRC"))
FeaturePlot_scCustom(ECstory, features = c("THY1"))
FeaturePlot_scCustom(ECstory, features = c("ITGAX"))
FeaturePlot_scCustom(ECstory, features = c("KLRB1"))
FeaturePlot_scCustom(ECstory, features = c("CD14"))
FeaturePlot_scCustom(ECstory, features = c("FCGR1A"))
FeaturePlot_scCustom(ECstory, features = c("CCR2"))
FeaturePlot_scCustom(ECstory, features = c("MARCO"))
FeaturePlot_scCustom(ECstory, features = c("SIGLEC1"))
FeaturePlot_scCustom(ECstory, features = c("SIGLEC1","ITGAX","MARCO"))
FeaturePlot_scCustom(ECstory, features = c("SIGLECL1"))
FeaturePlot_scCustom(ECstory, features = c("MRC1", "CD86"))
FeaturePlot_scCustom(ECstory, features = c("SIGLEC1", "ITGAM", "FCGR1A", "ITGAX", "HLA-DRA","ADGRE1"))
FeaturePlot_scCustom(ECstory, features = c("ADGRE1"))
FeaturePlot_scCustom(ECstory, features = c("MERTK"))
FeaturePlot_scCustom(ECstory, features = c("CD14", "FCGR3A"))
FeaturePlot_scCustom(ECstory, features = c("CD68"))
FeaturePlot_scCustom(ECstory, features = c("CD1A", "CD14", "CD19", "CD209"))
FeaturePlot_scCustom(ECstory, features = c("CD8A"))
FeaturePlot_scCustom(ECstory, features = c("GZMB"))
FeaturePlot_scCustom(ECstory, features = c("CADM1"))
FeaturePlot_scCustom(ECstory, features = c("FCER1A"))
FeaturePlot_scCustom(ECstory, features = c("SIGLECH"))
FeaturePlot_scCustom(ECstory, features = c("HLA-DRB1"))
FeaturePlot_scCustom(ECstory, features = c("IL3RA"))
FeaturePlot_scCustom(ECstory, features = c("MRC1"))
FeaturePlot_scCustom(ECstory, features = c("CD209"))
FeaturePlot_scCustom(ECstory, features = c("SIRPA"))
FeaturePlot_scCustom(ECstory, features = c("CLEC4C"))
FeaturePlot_scCustom(ECstory, features = c("ITGAX", "ITGAM", "CD1A", "CD1C"))
FeaturePlot_scCustom(ECstory, features = c("TNFRSF1A", "TNFRSF1B", "CD14", "FCGR3A","CCR2","CCR5"))
FeaturePlot_scCustom(ECstory, features = c("SELL"))
FeaturePlot_scCustom(ECstory, features = c("CCR2","CCR5","CXCR1","CX3CR1"))
FeaturePlot_scCustom(ECstory, features = c("CCR5"))
FeaturePlot_scCustom(ECstory, features = c("IL1RN"))



ECstory$clusters_tommy <- ""
ECstory$clusters_tommy[(ECstory$seurat_clusters == "6" )] <- "Endothelial 1"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "12")] <- "Endothelial 2"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "27")] <- "Endo MT"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "0" )] <- "Fibrobalsts"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "24")] <- "Myofibrobalsts"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "4" )] <- "Contractile SMC"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "14")] <- "Fibroblast like SMC"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "16")] <- "Oxygen sensing SMC" #pericytes seem to be a subset of these
ECstory$clusters_tommy[(ECstory$seurat_clusters == "9" )] <- "Synthetic SMC"    
ECstory$clusters_tommy[(ECstory$seurat_clusters == "21")] <- "B-cells" 
ECstory$clusters_tommy[(ECstory$seurat_clusters == "23")] <- "Platelets" 
ECstory$clusters_tommy[(ECstory$seurat_clusters == "20")] <- "Granulocytes"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "22")] <- "Erythrocytes"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "13")] <- "non Classical Monocytes"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "7" )] <- "Classical Monocytes"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "2" )] <- "Intermediate Monocytes"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "25")] <- "cDC"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "19")] <- "moDC"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "11")] <- "non-classical Monocytes"
ECstory$clusters_tommy[(ECstory$seurat_clusters == "15")] <- "intermediary Monocytes"
DimPlot(ECstory, group.by = "clusters_tommy", cols=polychrome_pal)


ECstory$clusters_tommy2 <- ""
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "6" )] <- "Endothelial cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "12")] <- "Endothelial cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "27")] <- "Endo MT"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "0" )] <- "Fibroblasts"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "24")] <- "Fibroblasts"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "4" )] <- "SMC"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "14")] <- "SMC"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "16")] <- "SMC" #pericytes seem to be a subset of these
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "9" )] <- "SMC"    
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "21")] <- "B-cells"  
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "23")] <- "Platelets" 
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "20")] <- "Granulocytes"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "22")] <- "Erythrocytes"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "13")] <- "Mono/Macs"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "7" )] <- "Mono/Macs"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "2" )] <- "Mono/Macs"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "25")] <- "Dendritic cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "19")] <- "Dendritic cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "18")] <- "Dendritic cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "11")] <- "Mono/Macs"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "5")] <- "Mono/Macs"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "17")] <- "Mono/Macs"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "15")] <- "Epithelial cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "26")] <- "Epithelial cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "8")] <- "NK cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "1")] <- "T-cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "3")] <- "T-cells"
ECstory$clusters_tommy2[(ECstory$seurat_clusters == "10")] <- "T-cells"

DimPlot(ECstory, group.by = "clusters_tommy2", cols=varibow_pal, split.by="samplecondition")

DimPlot(ECstory.subset.fil, group.by = "seurat_clusters")
DimPlot(ECstory.subset, group.by = "seurat_clusters")

DimPlot(ECstory.fil, group.by = "seurat_clusters", cols=polychrome_pal)

ECstory.fil$clusters_tommy3 <- ""
ECstory.fil$clusters_tommy4 <- ECstory.subset$seurat_clusters
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "27")] <- "Endo MT"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "0" )] <- "Fibroblasts"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "24")] <- "Fibroblasts"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "4" )] <- "SMC"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "14")] <- "SMC"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "16")] <- "SMC" #pericytes seem to be a subset of these
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "9" )] <- "SMC"    
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "21")] <- "B-cells"  
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "23")] <- "Lymphatic EC" 
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "20")] <- "Granulocytes"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "22")] <- "Erythrocytes"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "13")] <- "Mono/Macs"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "7" )] <- "Mono/Macs"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "2" )] <- "Mono/Macs"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "25")] <- "Dendritic cells"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "19")] <- "Dendritic cells"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "18")] <- "Dendritic cells"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "11")] <- "Mono/Macs"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "5")] <- "Mono/Macs"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "17")] <- "Mono/Macs"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "15")] <- "Epithelial cells"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "26")] <- "Epithelial cells"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "8")] <- "NK cells"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "1")] <- "T-cells"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "3")] <- "T-cells"
ECstory.fil$clusters_tommy3[(ECstory.fil$seurat_clusters == "10")] <- "T-cells"
ECstory.fil$clusters_tommy3[(ECstory.fil$clusters_tommy4 == "0")] <- "EC0"
ECstory.fil$clusters_tommy3[(ECstory.fil$clusters_tommy4 == "1")] <- "EC1"
ECstory.fil$clusters_tommy3[(ECstory.fil$clusters_tommy4 == "2")] <- "EC2"
ECstory.fil$clusters_tommy3[(ECstory.fil$clusters_tommy4 == "3")] <- "EC3"
ECstory.fil$clusters_tommy3[(ECstory.fil$clusters_tommy4 == "4")] <- "EC4"

DimPlot(ECstory.fil, group.by = "clusters_tommy3", cols=polychrome_pal)

DimPlot(ECstory.fil.subset.fil)
DimPlot(ECstory.fil.subset)

DimPlot(test, group.by = "clusters_tommy3")

DimPlot(ECstory.fil, group.by = "clusters_tommy2")
ECstory.cellchat <- subset(ECstory.fil, QCcut %in% c("Viable"))



Donor <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("Donor")) & (clusters_tommy3 %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))
PAH <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("IPAH")) & (clusters_tommy3 %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))
PF <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("IPF")) & (clusters_tommy3 %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))

PH <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("IPAH", "IPF")) & (clusters_tommy3 %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))


DimPlot(Donor, group.by = "clusters_tommy3")
DimPlot(PAH, group.by = "clusters_tommy3")
DimPlot(PF, group.by = "clusters_tommy3")
DimPlot(PH, group.by = "clusters_tommy3")



# DimPlot(ECstory.subset.fil, group.by = "namedclusters")
# Idents(ECstory.subset.fil) <- ECstory.subset.fil$namedclusters
# ECstory.subset.fil.Donor <- subset(ECstory.subset.fil, samplecondition %in% c("Donor"))
# ECstory.subset.fil.IPF <- subset(ECstory.subset.fil, samplecondition %in% c("IPF"))
# ECstory.subset.fil.PAH <- subset(ECstory.subset.fil, samplecondition %in% c("IPAH"))
# ECstory.subset.fil.PH <- subset(ECstory.subset.fil, samplecondition %in% c("IPAH","IPF"))


data.input <- GetAssayData(PH, assay = "RNA", slot = "data") # normalized data matrix
cellchat <- createCellChat(object = PH, group.by = "clusters_tommy3", assay = "RNA")
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data #PPI.mouse for protein prtoein interactions
showDatabaseCategory(CellChatDB)


cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, nboot = 20)
#saveRDS(cellchat, file = "cellchat_Fra2.rds")

#cellchat <- readRDS( file = "C:/Users/o_lins/Documents/cellchat_Fra2.rds")


cellchat <- filterCommunication(cellchat, min.cells = 5)

df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

#pal3 <- c("#7ED7D1FF", "#325A9BFF", "#F8A19FFF", "#F6222EFF", "#90AD1CFF", "#85660DFF", "#FA0087FF", "#2ED9FFFF", "#16FF32FF", "#FBE426FF", "#683B79FF","#038637","#f39d10")
# pal2 <- c("#7ED7D1FF", "#325A9BFF", "#F8A19FFF", "#F6222EFF", "#90AD1CFF", "#85660DFF", "#FA0087FF", "#2ED9FFFF", "#16FF32FF", "#FBE426FF")
# #pal3 <- c("#3640ae","#d64141","#038637")
# pal3 <- c("#038637","#3640ae","#d64141")
# #pal3 <- c("#d64141","#038637","#3640ae")
# pal3 <- c("#c7c7c8", "#6e79d6", "#8d02e1", "#6ed67c", "#c6d66e", "#d6bb6e","#6eb1d6" , "#d66ec5","#d66e6e", "#7b7b7e", "#683B79FF","#038637","#f39d10")
# pal3 <- c("#dddddd", "#b0b5de", "#8d02e1", "#99d9a2", "#dedf70", "#d6bb6e","#67c2df" , "#d66ec5","#e16659", "#7b7b7e", "#683B79FF","#038637","#f39d10")
# pal3 <- c("#dddddd", "#b0b5de", "#8d02e1", "#99d9a2", "#dedf70", "#d6bb6e","#67c2df" , "#d66ec5","#e16659", "#7b7b7e")
pal3 <- c("#7ED7D1FF", "#683B79FF", "#F8A19FFF", "#16FF32FF", "#FA0087FF", "#325A9BFF", "#2ED9FFFF","#f39d10","#038637", "#F6222EFF","#FBE426FF","#f39d10")



# pal3 <- polychrome_pal
# names(pal3) <- unique(PH$clusters_tommy2)
# names(pal3) <- unique(ECstory.subset.fil$namedclusters2)
# names(pal3) <- unique(ECstory.sc.Donor$clusters_tommy.sc)
# names(pal3) <- unique(EC.story@)
names(pal3) <- unique(Donor$clusters_tommy3)




groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count,top=0.3,color.use =pal3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=0.3,color.use =pal3,vertex.weight = groupSize, weight.scale = T, title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)

netVisual_circle(cellchat@netP$pathways,top=0.2,color.use =pal3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.7)

netVisual_circle(cellchat@net$count,top=0.5,color.use =pal3,sources.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=0.5,color.use =pal3,sources.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)

netVisual_circle(cellchat@net$count,top=1,color.use =pal3,targets.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=1,color.use =pal3,targets.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)


netVisual_bubble(cellchat, sources.use = 3, targets.use = ECstory.sc.Donor$clusters_tommy.sc, remove.isolate = TRUE, sort.by.source = T)
 netVisual_chord_gene(cellchat, sources.use = 3, targets.use = ECstory.sc.Donor$clusters_tommy.sc, lab.cex = 0.5, slot.name = "netP")
 netVisual_chord_gene(cellchat, sources.use = ECstory.sc.Donor$clusters_tommy.sc, targets.use =3, lab.cex = 0.5, slot.name = "netP")
 
print(p)

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], arrow.size = 0.1, layout = "chord")
}
netVisual_circle(cellchat@net$count,top=0.2,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.65, layout = "chord")
netVisual_aggregate(cellchat, layout = "chord")
par(mfrow = c(1,2))
netVisual_aggregate(cellchat, signaling = c("CXCL"), layout="chord")
#cellchat_PAH_allcells <- cellchat
#cellchat_PHPF_allcells <- cellchat
#cellchat_Donor_allcells <- cellchat
#cellchat_PH_allcells <- cellchat

##CellChat merge

object.list <- list(Donor = cellchat_Donor_allcells, PAH = cellchat_PAH_allcells)
object.list <- list(Donor = cellchat_Donor_allcells, PF = cellchat_PHPF_allcells)
object.list <- list(Donor = cellchat.subset.Donor, PF = cellchat.subset.PHPF)
object.list <- list(Donor = cellchat.subset.Donor, PAH = cellchat.subset.PAH)
object.list <- list(Donor = cellchat_Donor_allcells, PH = cellchat_PH_allcells)
cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

cellchat.merge@meta$clusters_tommy3
unique(cellchat.merge@meta$clusters_tommy3)

levels(cellchat.merge@meta$clusters_tommy3)
cellchat.merge@meta$clusters_tommy3 <- factor(cellchat.merge@meta$clusters_tommy3, 
                                              levels = c("EC0", "EC1", "EC2",
                                                         "Fibroblasts", "SMC",
                                                         "Mono/Macs", "NK cells", "T-cells", "B-cells", "Granulocytes", "Dendritic cells"))

levels(cellchat.merge@meta$clusters_tommy3)


names(pal3) <- unique(cellchat.merge@meta$clusters_tommy3)
# names(pal3) <- unique(cellchat.merge@meta$clusters_tommy.sc)
# names(pal3) <- unique(cellchat.merge@meta$namedclusters)
# pal3 <- c("#038637","#3640ae","#d64141")
# pal3 <- polychrome_pal
# names(pal3) <- unique(ECstory.subset.fil$namedclusters2)
# names(pal3) <- unique(ECstory.sc$clusters_tommy.sc)

par(mfrow=c(1,1))
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("count"), color.use = pal3)
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("weight"), color.use = pal3)
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("count"))
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("weight"))


pal3 <- c("#7ED7D1FF", "#683B79FF", "#F6222EFF", "#325A9BFF", "#038637", "#16FF32FF", "#2ED9FFFF","#f39d10","#FA0087FF", "#F8A19FFF","#FBE426FF","#f39d10")


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.merge, weight.scale = T, color.use = pal3)

netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, sources.use = 3, arrow.size=0.6)
netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, sources.use = 3, vertex.label.cex = 0.01, arrow.size=0.6)

netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, targets.use = 3, arrow.size=0.6)
netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, targets.use = 3, vertex.label.cex = 0.01, arrow.size=0.6)

netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, sources.use = 4, arrow.size=0.6)
netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, sources.use = 4, vertex.label.cex = 0.01, arrow.size=0.6)

netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, targets.use = 4, arrow.size=0.6)
netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, targets.use = 4, vertex.label.cex = 0.01, arrow.size=0.6)


netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, sources.use = 5, arrow.size=0.6)
netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, sources.use = 5, vertex.label.cex = 0.01, arrow.size=0.6)

netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, targets.use = 5, arrow.size=0.6)
netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight", color.use = pal3, targets.use = 5, vertex.label.cex = 0.01, arrow.size=0.6)

pal3 <- c("#7ED7D1FF", "#683B79FF", "#F8A19FFF", "#16FF32FF", "#FA0087FF", "#325A9BFF", "#2ED9FFFF","#f39d10","#038637", "#F6222EFF","#FBE426FF","#f39d10")

pal3 <- c("#7ED7D1FF", "#683B79FF", "#F8A19FFF", "#16FF32FF", "#FA0087FF", "#325A9BFF", "#2ED9FFFF","#f39d10","#038637", "#F6222EFF","#FBE426FF","#f39d10")



Mono/Macs         T-cells Dendritic cells             SMC     Fibroblasts             EC1         B-cells        NK cells             EC0 
"#7ED7D1FF"     "#683B79FF"     "#F8A19FFF"     "#16FF32FF"     "#FA0087FF"     "#325A9BFF"     "#2ED9FFFF"       "#f39d10"       "#038637" 
EC2    Granulocytes            <NA> 
  "#F6222EFF"     "#FBE426FF"       "#f39d10" 

netVisual_circle(cellchat.merge@net$count,top=0.3,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat.merge@net$weight,top=0.3,color.use =pal3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)


netVisual_bubble(cellchat.merge, sources.use = 3, targets.use = NULL ,  comparison = c(1, 2), title.name = "Increased signaling in LS", 
                 angle.x = 45, remove.isolate = F,color.text=c("black","#0491b3"), font.size = 8,color.heatmap = c("Spectral"), max.dataset = NULL )
  

netVisual_bubble(cellchat.merge, sources.use = NULL, targets.use = NULL,  comparison = NULL, max.dataset = NULL, title.name = "Increased signaling in LS", 
                 angle.x = 45, remove.isolate = F,color.text=c("black","#0491b3"), font.size = 8)

"#e06606","#0491b3"
targets <- c("Endothelial cells", "Epithelial cells", "Fibroblasts", "SMC")

netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("count"), remove.isolate = T)
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("weight"), remove.isolate = T)

ptm = Sys.time()
gg1 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1,2), color.use =c("grey","#e06606") )
gg2 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1,2), measure = "weight", color.use =c("grey","#e06606"))
gg1 + gg2

pos.dataset = "LS"
cellchat.merge <- identifyOverExpressedGenes(cellchat.merge, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
features.name = paste0(pos.dataset, ".merged")

net <- netMappingDEG(cellchat.merge, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat.merge, net = net, datasets = "LS",ligand.logFC = 0.05, receptor.logFC = NULL)

#nichenet----------------------

  
  
  seuratObj <- JoinLayers(ECstory)

  Idents(seuratObj) <- seuratObj$clusters_tommy3

organism = "human"

if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}


lr_network = lr_network %>% distinct(from, to)
head(lr_network)

weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network

###1. Define a “sender/niche” cell population and a “receiver/target”
## receiver
receiver = c("Club")
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10, assay_oi = "RNA")

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Basal","ATI","ATII", "Ciliated", "Club", "Goblet", "Mesothelial", "PNEC")

#"SMC", "T-cells", "NK cells", "Mono/Macs",
#"Granulocytes", "Fibroblasts", "Erythrocytes", "Epithelial cells",
#"Dendritic cells", "B-cells"

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10, assay_oi ="RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

###2. Define a gene set of interest:
seurat_obj_receiver= subset(seuratObj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["Manuscript_Identity_red", drop=TRUE]])

condition_oi = "IPAH"
condition_reference = "Donor"

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% dplyr::filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

###3. Define a set of potential ligands:
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

###4. Perform NicheNet ligand activity analysis: 
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-aupr_corrected) #%>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

#These ligands are expressed by one or more of the input sender cells. To see which cell population expresses which of these top-ranked ligands, you can run the following:
DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

DotPlot_scCustom(seuratObj, assay= "RNA", features = best_upstream_ligands %>% rev(), colors_use = c("#4f4ceb", "#807ee9","#a8a7e6","#d8d8d8","#f38d8d","#f25656","#da0000"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  #geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  RotatedAxis()+coord_flip()

###5) Infer receptors and top-predicted target genes of ligands  
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")+ theme(legend.text=element_text(angle=45,hjust=1.2, vjust=1.4),legend.title = element_text(vjust = 0.9), axis.text.x = element_text(size=6))

###D
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

###6) Add log fold change information of ligands from sender cells
# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratObj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratObj, condition_colname = "timetype2", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

# change colors a bit to make them more stand out
p_ligand_lfc = p_ligand_lfc + scale_fill_gradientn(colors = c("midnightblue","blue", "grey95", "grey99","firebrick1","red"),values = c(0,0.1,0.2,0.25, 0.40, 0.7,1), limits = c(vis_ligand_lfc %>% min() - 0.1, vis_ligand_lfc %>% max() + 0.1)) +theme(legend.text = element_text(size = 9, angle =45, hjust=1.2, vjust=1.4), legend.title = element_text(vjust = 0.9))
p_ligand_lfc


###7) Summary visualizations of the NicheNet analysis
# ligand activity heatmap
ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()

vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR") + theme(legend.text = element_text(size = 9, angle =45, hjust=1.2, vjust=1.4), legend.title = element_text(vjust = 0.9))

# ligand expression Seurat dotplot
order_ligands_adapted <- str_replace_all(order_ligands, "\\.", "-")
rotated_dotplot = DotPlot(seuratObj %>% subset(clusters_tommy3 %in% sender_celltypes), features = order_ligands_adapted) + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
figures_without_legend = cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  #p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  #ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

p_ligand_aupr +
  guide_legend()

###add aditional data continue here instead of "5)"
# By default, ligand_condition_specificty and receptor_condition_specificty are 0
prioritizing_weights = c("de_ligand" = 1,
                         "de_receptor" = 1,
                         "activity_scaled" = 2,
                         "exprs_ligand" = 1,
                         "exprs_receptor" = 1,
                         "ligand_condition_specificity" = 0.5,
                         "receptor_condition_specificity" = 0.5)

#lr_network <- lr_network[!is.na(lr_network)]
lr_network_renamed <- lr_network %>% dplyr::rename(ligand=from, receptor=to)
#colnames(lr_network_renamed) <- c("ligand","receptor")

# Only calculate DE for LCMV condition, with genes that are in the ligand-receptor network
DE_table <- calculate_de(seuratObj, celltype_colname = "clusters_tommy3",
                         condition_colname = "samplecondition", condition_oi = condition_oi,
                         features = union(expressed_ligands, expressed_receptors))

# Average expression information - only for LCMV condition
expression_info <- get_exprs_avg(seuratObj, "clusters_tommy3", condition_colname = "samplecondition", condition_oi = condition_oi)



# Calculate condition specificity - only for datasets with two conditions!
condition_markers <- FindMarkers(object = seuratObj, ident.1 = condition_oi, ident.2 = condition_reference,
                                 group.by = "samplecondition", min.pct = 0, logfc.threshold = 0,
                                 features = union(expressed_ligands, expressed_receptors)) %>% rownames_to_column("gene")

# Combine DE of senders and receivers -> used for prioritization
processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network_renamed,
                                          senders_oi = sender_celltypes, receivers_oi = receiver)

processed_expr_table <- process_table_to_ic(expression_info, table_type = "expression", lr_network_renamed)

processed_condition_markers <- process_table_to_ic(condition_markers, table_type = "group_DE", lr_network_renamed)

prior_table <- generate_prioritization_tables(processed_expr_table,
                                              processed_DE_table,
                                              ligand_activities,
                                              processed_condition_markers,
                                              prioritizing_weights = prioritizing_weights)

prior_table 

#only relevant columns:
prior_table %>% select(c('sender', 'receiver', 'ligand', 'receptor', 'scaled_lfc_ligand', 'scaled_lfc_receptor', 'scaled_p_val_ligand_adapted', 'scaled_p_val_receptor_adapted', 'scaled_avg_exprs_ligand', 'scaled_avg_exprs_receptor', 'scaled_lfc_ligand_group', 'scaled_lfc_receptor_group', 'scaled_activity'))


best_upstream_ligands = ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(desc(aupr_corrected)) %>% pull(test_ligand) %>% unique()

# DE analysis for each sender cell type
DE_table_all = Idents(seuratObj) %>% levels() %>% intersect(sender_celltypes) %>%
  lapply(get_lfc_celltype, seurat_obj = seuratObj, condition_colname = "samplecondition", condition_oi = condition_oi, condition_reference = condition_reference,
         expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) 
DE_table_all[is.na(DE_table_all)] = 0


order_ligands <- make.names(best_upstream_ligands) %>% rev()

# ligand activity heatmap
ligand_aupr_matrix <- ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand) %>%
  `rownames<-`(make.names(rownames(.))) %>% `colnames<-`(make.names(colnames(.)))

vis_ligand_aupr <- as.matrix(ligand_aupr_matrix[order_ligands, ], ncol=1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr, "Prioritized ligands","Ligand activity",
                                     color = "darkorange",legend_position = "top", x_axis_position = "top",
                                     legend_title = "AUPR") +
  theme(legend.text = element_text(size = 9, angle =45, hjust=1.2, vjust=1.4), legend.title = element_text(vjust = 0.9))



# LFC heatmap
# First combine ligand activities with DE information and make 
ligand_activities_de <- ligand_activities %>% select(test_ligand, aupr_corrected) %>% dplyr::rename(ligand = test_ligand) %>% left_join(DE_table_all %>% dplyr::rename(ligand = gene))
# ligand_activities_de <- ligand_activities
# ligand_activities_de <- ligand_activities_de %>% select(test_ligand, aupr_corrected)
# colnames(ligand_activities_de) <- c("ligand", "aupr_corrected")
# colnames(DE_table_all) <- c("gene" ,"early fibrotic gCap" ,"aCap",  "gCap" ,"fibrotic aCap" ,"CD45+/CD31+"  ,"Vein" ,"Artery","late fibrotic gCap")
# ligand_activities_de %>% left_join(DE_table_all)



ligand_activities_de[is.na(ligand_activities_de)] <- 0
lfc_matrix <- ligand_activities_de  %>% select(-ligand, -aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand) %>%
  `rownames<-`(make.names(rownames(.))) %>% `colnames<-`(make.names(colnames(.)))
vis_ligand_lfc <- lfc_matrix[order_ligands,]

p_ligand_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc, "Prioritized ligands","LFC in Sender",
                                               low_color = "midnightblue", mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",
                                               legend_position = "top", x_axis_position = "top", legend_title = "LFC") +
  theme(axis.text.y = element_text(face = "italic"))


# ligand expression Seurat dotplot
order_ligands_adapted <- str_replace_all(order_ligands, "\\.", "-")


rotated_dotplot <- DotPlot(seuratObj %>% subset(clusters_tommy3 %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") +
  # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
  coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 10), axis.text.x = element_text(size=12, angle=45,hjust=1, vjust=1))+
  scale_x_discrete(position = "top") 

rotated_dotplot #+ theme(axis.text.x = element_text(angle=90))

x_axis_position = "top"

# Combine figures and legend separately
figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12),
                          axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) +
    ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text(size = 9, angle =0)) + ylab(""),
  align = "hv",
  nrow = 1)
# rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8))


legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1))

combined_plot <- cowplot::plot_grid(figures_without_legend, legends, nrow = 2, align = "hv")
print(combined_plot)

mylevels <- c("gCap","early fibrotic gCap","late fibrotic gCap", "aCap", "fibrotic aCap", "Vein", "Artery", "CD45+/CD31+","dESC", "Lymphatic")


prior_table2 <- prior_table
prior_table2 <- subset(prior_table, sender %in% mylevels)
prior_table2$ID <- ""
prior_table2$ID <- prior_table2$sender
prior_table2$ID[prior_table2$ID == "gCap"] <- "1"
prior_table2$ID[prior_table2$ID == "early fibrotic gCap"] <- "2"
prior_table2$ID[prior_table2$ID == "late fibrotic gCap"] <- "3"
prior_table2$ID[prior_table2$ID == "aCap"] <- "4"
prior_table2$ID[prior_table2$ID == "fibrotic aCap"] <- "5"
prior_table2$ID[prior_table2$ID == "Vein"] <- "6"
prior_table2$ID[prior_table2$ID == "Artery"] <- "7"
#prior_table2$ID[prior_table2$ID == "CD45+/CD31+"] <- "8"
#prior_table2$ID[prior_table2$ID == "dESC"] <- "9"
#prior_table2$ID[prior_table2$ID == "Lymphatic"] <- "10"


prior_table2$ID <- ""
prior_table2$ID <- prior_table2$sender
prior_table2$ID[prior_table2$ID == "T-cell"] <- "1"
prior_table2$ID[prior_table2$ID == "B-cell"] <- "2"
prior_table2$ID[prior_table2$ID == "Dendritic"] <- "3"
prior_table2$ID[prior_table2$ID == "Granulocytes"] <- "4"
prior_table2$ID[prior_table2$ID == "alv. Mac"] <- "5"
prior_table2$ID[prior_table2$ID == "int. Mac"] <- "6"
prior_table2$ID[prior_table2$ID == "Mesenchymal"] <- "7"



prior_table2$sender <- factor(prior_table2$sender)

make_mushroom_plot(prior_table, top_n = 30) + theme(axis.title.x = element_text(size = 0, angle =0), axis.text.x = element_text(angle =90))+theme(
  legend.text = element_text(size = 5,hjust = 0.5, lineheight = 121), 
  legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))#+NoLegend()




active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")+ theme(legend.text=element_text(angle=45,hjust=1.2, vjust=1.4),legend.title = element_text(vjust = 0.9))

#signaling network
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final_mouse.rds"))

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
sig_network <- readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_mouse_21122021.rds"))
gr_network <- readRDS(url("https://zenodo.org/record/7074291/files/gr_network_mouse_21122021.rds"))


ligands_oi <- c("Il1b", "Tgfb1","Vgfa") # this can be a list of multiple ligands if required
targets_oi <- c("Adamts1","Dusp6", "Edn1", "Egr1", "Emp1", "Ets1","Flt1", "Mef2c", "Mgp", "Ndrg1", "Plat")

ligands_oi <- c("App", "Tgfb1","Vgfa") # this can be a list of multiple ligands if required
targets_oi <- c("Adamts1","Dusp6", "Edn1", "Egr1", "Emp1", "Ets1","Flt1", "Mef2c", "Mgp", "Ndrg1", "Plat")

active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,
                                                      targets_all = targets_oi,
                                                      weighted_networks = weighted_networks,
                                                      ligand_tf_matrix = ligand_tf_matrix,
                                                      top_n_regulators = 4,
                                                      minmax_scaling = TRUE) 


graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                   ligands_all = ligands_oi, targets_all = targets_oi,
                                                   sig_color = "indianred", gr_color = "steelblue")

# To render the graph in RStudio Viewer, uncomment following line of code
# DiagrammeR::render_graph(graph_min_max, layout = "tree")

# To export/draw the svg, you need to install DiagrammeRsvg
graph_svg <- DiagrammeRsvg::export_svg(DiagrammeR::render_graph(graph_min_max, layout = "tree", output = "graph"))
cowplot::ggdraw() + cowplot::draw_image(charToRaw(graph_svg))



# Read in expression data and update
seuratObj <- readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
seuratObj <- UpdateSeuratObject(seuratObj)
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

# Run Liana
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse") # convert gene names

# Load LIANA mouse consensus ligand-receptor network
lr_network_liana <- select_resource("MouseConsensus")[[1]] %>% decomplexify()
lr_network_liana <- lr_network_liana %>% rename(from = source_genesymbol, to = target_genesymbol) %>% select(from, to) %>% mutate(source = "liana")

# Define receiver cell type
receiver = "CD8 T"
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)

# Define sender cell types
sender_celltypes <- c("CD4 T","Treg", "Mono", "NK", "B", "DC")

# Get expressed genes in the sender
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

# Get gene set of interest through DE analysis
seurat_obj_receiver <- subset(seuratObj, idents = receiver)
condition_oi <- "LCMV"; condition_reference <- "SS"
DE_table_receiver <- FindMarkers(object = seurat_obj_receiver,
                                 ident.1 = condition_oi, ident.2 = condition_reference,
                                 group.by = "aggregate", 
                                 min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)

# Get potential ligands
ligands <- lr_network_liana %>% pull(from) %>% unique()
receptors <- lr_network_liana %>% pull(to) %>% unique()
expressed_ligands <- intersect(ligands,expressed_genes_sender)
expressed_receptors <- intersect(receptors,expressed_genes_receiver)
potential_ligands <- lr_network_liana %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
potential_ligands
##  [1] "Adam10" "Adam23" "Icosl"  "Thy1"   "H2-K1"  "App"    "Itgb2"  "Icam1"  "Lck"    "Fcer2a" "Cd48"   "Icam2"  "H2-M3"  "Vcam1"  "Cd22"   "Cd86"  
## [17] "B2m"    "Sirpa"  "Adam17" "Selplg" "Btla"   "F11r"   "Fn1"    "Ebi3"   "Hp"     "Tgfb1"  "Lgals1" "Vcan"   "Mif"    "Il15"   "Copa"   "Ybx1"  
## [33] "Fam3c"  "Ccl2"   "Crlf2"  "Ccl22"  "Cxcl10" "Cxcl16" "Cd72"   "Txlna"  "Psen1"  "F13a1"  "C1qb"   "Gnai2"  "Tgm2"   "Gnas"   "Ccl5"   "Mia"   
## [49] "H2-T23" "Lyz1"   "Lyz2"

# Constructing the ligand-target matrix
gr_network_mouse <- readRDS(url("https://zenodo.org/records/7074291/file/gr_network_mouse_21122021.rds"))

# Define optimum weight as the one calculated for the NicheNet LR network
optim_weight <- optimized_source_weights_df %>% filter(source == "nichenet_verschueren") %>% pull(avg_weight) # 0.2788104
source_weights <- optimized_source_weights_df %>%
  add_row(source = "liana", avg_weight = optim_weight, .before = 1) %>% rename(weight=avg_weight)

# Aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks_liana <- construct_weighted_networks(lr_network = lr_network_liana, sig_network = sig_network_mouse, gr_network = gr_network_mouse, source_weights_df = source_weights)

# Downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks_liana <- apply_hub_corrections(weighted_networks = weighted_networks_liana,
                                                 lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                                 gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight))

# Construct ligand-target matrix using only the potential ligands to save time
ligand_target_matrix_liana <- construct_ligand_target_matrix(weighted_networks = weighted_networks_liana, ligands = as.list(potential_ligands), algorithm = "PPR",
                                                             damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                             ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))

# Filter background genes and the gene set of interest to only the ones in the ligand-target matrix
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix_liana)]
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix_liana)]

# Perform ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix_liana, potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
## # A tibble: 51 × 6
##    test_ligand auroc  aupr aupr_corrected pearson  rank
##    <chr>       <dbl> <dbl>          <dbl>   <dbl> <dbl>
##  1 Ebi3        0.663 0.389          0.244  0.300      1
##  2 H2-M3       0.608 0.292          0.147  0.179      2
##  3 H2-T23      0.611 0.278          0.133  0.153      3
##  4 Lck         0.621 0.268          0.122  0.103      4
##  5 H2-K1       0.605 0.268          0.122  0.141      5
##  6 Sirpa       0.613 0.263          0.117  0.143      6
##  7 Cd48        0.612 0.257          0.112  0.0860     7
##  8 Tgfb1       0.597 0.254          0.108  0.202      8
##  9 Ccl22       0.605 0.250          0.104  0.125      9
## 10 Ccl5        0.606 0.248          0.102  0.0343    10
## # ℹ 41 more rows

###circle plot visualization
library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(circlize)

ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
head(lr_network)



seuratObj <- fra2dataDFr_integrated
seuratObj <- JoinLayers(seuratObj)
#seuratObj[["RNA"]] <- as(object = seuratObj[["RNA"]], Class = "Assay")
Idents(seuratObj) <- seuratObj$nichenet


# Convert gene names
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

sender_celltypes <- "Endothelial"
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = c("Dendritic","Granulocytes","B-cell", "T-cell", "Monocytes", "Mesenchymal", "Epithelial"), 
  condition_colname = "timetype2",
  condition_oi = "TG8w", condition_reference = "WT8w", 
  sender = sender_celltypes, 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks,
  top_n_ligands = 20)
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"

nichenet_output$ligand_activities
nichenet_output$ligand_expression_dotplot
nichenet_output$ligand_differential_expression_heatmap
VlnPlot(seuratObj, features = c("Ptprc", "H2-M3", "Cxcl10"), split.by = "timetype2", pt.size = 0, combine = TRUE)
nichenet_output$ligand_target_heatmap

ligand_type_indication_df <- assign_ligands_to_celltype(seuratObj,
                                                        nichenet_output$top_ligands,
                                                        celltype_col = "nichenet") 

ligand_type_indication_df %>% head()

ligand_type_indication_df$ligand_type %>% table()


head(nichenet_output$ligand_target_df)

active_ligand_target_links_df <- nichenet_output$ligand_target_df
active_ligand_target_links_df$target_type <- "Endothelial" # needed for joining tables
circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
                                           active_ligand_target_links_df,
                                           cutoff = 0.90) 

head(circos_links)



ligand_colors <- c("General" = "#377EB8", "B-cell" = "#4DAF4A", "T-cell" = "#984EA3",
                   "Monocytes" = "#FF7F00", "Mesenchymal" = "#c7ea38", "Granulocytes" = "#F781BF",
                   "Dendritic"= "#E41A1C", "Epithelial" = "#38d2ea", "Endothelial" ="#999999") 
target_colors <- c("Endothelial" = "#999999") 

ligand_colors <- c("Endothelial" = "#999999") 

target_colors <- c("General" = "#377EB8", "B-cell" = "#4DAF4A", "T-cell" = "#984EA3",
                   "Monocytes" = "#FF7F00", "Mesenchymal" = "#c7ea38", "Granulocytes" = "#F781BF",
                   "Dendritic"= "#E41A1C", "Epithelial" = "#38d2ea", "Endothelial" ="#999999") 

vis_circos_obj <- prepare_circos_visualization(circos_links,
                                               ligand_colors = ligand_colors,
                                               target_colors = target_colors,
                                               celltype_order = NULL) 


make_circos_plot(vis_circos_obj, transparency = FALSE, args.circos.text = list(cex = 0.5)) 


make_circos_plot(vis_circos_obj, transparency = TRUE, args.circos.text = list(cex = 0.5))


par(bg = "transparent")

# Default celltype order
celltype_order <- unique(circos_links$ligand_type) %>% sort() %>% .[. != "General"] %>% c(., "General")

# Create legend
circos_legend <- ComplexHeatmap::Legend(
  labels = celltype_order,
  background = ligand_colors[celltype_order],
  type = "point",
  grid_height = unit(3, "mm"),
  grid_width = unit(3, "mm"),
  labels_gp = grid::gpar(fontsize = 8)
)

circos_legend_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(circos_legend))

make_circos_plot(vis_circos_obj, transparency = TRUE, args.circos.text = list(cex = 0.5))
p_circos_no_legend <- recordPlot()


cowplot::plot_grid(p_circos_no_legend, circos_legend_grob, rel_widths = c(1, 0.1))


svg("ligand_target_circos.svg", width = 10, height = 10)
cowplot::plot_grid(p_circos_no_legend, circos_legend_grob, rel_widths = c(1, 0.1))
dev.off()


lr_network_top_df <- nichenet_output$ligand_receptor_df %>%
  mutate(target_type = "LCMV_CD8T_receptor") %>%
  rename(target=receptor) %>%
  inner_join(ligand_type_indication_df)

receptor_colors <- c("LCMV_CD8T_receptor" = "#E41A1C")

vis_circos_receptor_obj <- prepare_circos_visualization(lr_network_top_df,
                                                        ligand_colors = ligand_colors,
                                                        target_colors = receptor_colors) 


make_circos_plot(vis_circos_receptor_obj, transparency = FALSE,
                 link.visible = TRUE,  args.circos.text = list(cex = 0.8))











ECstory.subset.fil$namedclusters2 <- ""
ECstory.subset.fil$namedclusters2[ECstory.subset.fil$namedclusters == "Immune modulatory"] <- "Immune modulatory"
ECstory.subset.fil$namedclusters2[ECstory.subset.fil$namedclusters == "Homeostatic"] <- "Homeostatic"
ECstory.subset.fil$namedclusters2[ECstory.subset.fil$namedclusters == "Vaso active"] <- "Defense"

DimPlot(ECstory.subset.fil, group.by = "namedclusters2")


netVisual_circle(cellchat@net$weight,top=0.5,color.use =pal3,vertex.weight = groupSize, weight.scale = T, title.name = "Interaction weights/strength", vertex.label.cex=0.8, idents.use = )

unique(ECstory$clusters_tommy2)

#sctalk--------




ECstory.sc <- JoinLayers(ECstory) 

ECstory.sc $clusters_tommy.sc<- ""
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "Endo1"] <- "Endothelial"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "Endo2"] <- "Endothelial"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "Erythrocytes"] <- "Erythrocytes"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "B-cells"] <- "B-cells"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "Dendritic cells"] <- "Dendritic Cells"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "Epithelial Cells"] <- "Epithelial Cells"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "Granulocytes"] <- "Granulocytes"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "Lymphatic Endo"] <- "Lymphatic Endo"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "Fibroblasts"] <- "Fibroblasts"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "Mono/Macs"] <- "Mono/Macs"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "NK-cells"] <- "NK-cells"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "SMC1"] <- "SMC"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "SMC2"] <- "SMC"
ECstory.sc $clusters_tommy.sc[ECstory.sc $clusters_tommy5 == "T-cells"] <- "T-cells"




Idents(ECstory.sc) <- ECstory.sc$clusters_tommy.sc
DimPlot(ECstory.sc, group.by = "clusters_tommy.sc")
FeaturePlot(ECstory.sc, features = "ACKR1" )
populations.use <- ECstory.sc$clusters_tommy.sc

ECstory.sc.Donor <- subset(ECstory.sc, samplecondition %in% c("Donor"))



file.lab <- "TIP"

GenerateEdgeWeights(seurat.object = ECstory.sc,
                    file.label = file.lab,
                    species = "human",
                    verbose=T)


edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)


GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)
##
ECstory.subset.fil.sc <- JoinLayers(ECstory.subset.fil) 
Idents(ECstory.subset.fil.sc) <- ECstory.subset.fil.sc$namedclusters
DimPlot(ECstory.subset.fil.sc, group.by = "namedclusters")
populations.use <- ECstory.subset.fil.sc$namedclusters

file.lab <- "TIP2"

GenerateEdgeWeights(seurat.object = ECstory.subset.fil.sc,
                    file.label = file.lab,
                    species = "human",
                    verbose=T)


edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)


GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)

paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)


results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4,
                                     return.results = TRUE)
head(results.table)

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")
col.set <- polychrome_pal
col.set <- c("#c10023","#0099cc","#008e17", "#fde800","#aa00ff","#5f777f","#f2a287")
col.set <- c("#c10023", "#5f777f", "#00cc99", "#f60000", "#c10023", "#4ffc00","#fb8500", "#0053c8", "#0099cc", "#aa00ff", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

#p.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO","EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc", "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")

p.labels <- unique(ECstory.sc$clusters_tommy.sc)

names(col.set) <- p.labels

## Generate a circle plot
cell.network.file = paste0("Permutation_tests_", file.lab, "_network.csv")
CellCirclePlot(input.file = cell.network.file, 
               adj_pval_thresh = 0.01, 
               col.set = col.set)

knitr::include_graphics('TIP_small_circle_plot.png')

## Generate an inbound vs outbound plot
InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05)
knitr::include_graphics('TIP_small_inbound_outbound.png')

InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05,
                    lab.weight.thresh = 10000)

## Plot top ligands for a cell population
path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
cell.type <- "Endothelial"
PlotTopLigands(input.file = path.file, 
               col.use = "#0053c8",
               cell.identity = cell.type)
knitr::include_graphics('TopLigands.png')

## First, select ligands upregulated in source population
source.marker.table <- Seurat::FindMarkers(ECstory.sc, ident.1 = "Endothelial", 
                                           logfc.threshold = 1.5, 
                                           only.pos = TRUE)
source.marker.genes <- rownames(source.marker.table)

path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
edge.score.file <- paste0(file.lab, "_all_ligand_receptor_network_edges.csv")
source.population <- "Endothelial"
target.populations <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")

source.population <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
target.populations <- "Endothelial"
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")


source.population <- c("Endothelial")
target.populations <- c("Fibroblasts", "SMC")
population.cols <- c("#0053c8", "#fb8500", "#c10023")

source.population <- "Endothelial"
target.populations <- c("Lymphatic Endo","Mono/Macs")
population.cols <- c("#aa00ff", "#5f777f","#0053c8")

source.population <- c("Endothelial")
target.populations <- c("Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f")

  NetworkTreePlot(path.file = path.file, 
                edge.score.file = edge.score.file,
                source.population = source.population,
                target.populations = target.populations,
                source.marker.genes = source.marker.genes,
                population.cols = population.cols)

knitr::include_graphics('tree_plot.png')


#----------
file.lab <- "TIP3"

GenerateEdgeWeights(seurat.object = ECstory.sc.Donor,
                    file.label = file.lab,
                    species = "human",
                    verbose=T)


edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)


GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)

paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)


results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4,
                                     return.results = TRUE)
head(results.table)

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")
col.set <- polychrome_pal
col.set <- c("#c10023","#0099cc","#008e17", "#fde800","#aa00ff","#5f777f","#f2a287")
col.set <- c("#c10023", "#5f777f", "#00cc99", "#f60000", "#c10023", "#4ffc00","#fb8500", "#0053c8", "#0099cc", "#aa00ff", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

#p.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO","EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc", "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")

p.labels <- unique(ECstory.sc$clusters_tommy.sc)

names(col.set) <- p.labels

## Generate a circle plot
cell.network.file = paste0("Permutation_tests_", file.lab, "_network.csv")
CellCirclePlot(input.file = cell.network.file, 
               adj_pval_thresh = 0.01, 
               col.set = col.set)

knitr::include_graphics('TIP_small_circle_plot.png')

## Generate an inbound vs outbound plot
InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05)
knitr::include_graphics('TIP_small_inbound_outbound.png')

InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05,
                    lab.weight.thresh = 10000)

## Plot top ligands for a cell population
path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
cell.type <- "Endothelial"
PlotTopLigands(input.file = path.file, 
               col.use = "#0053c8",
               cell.identity = cell.type)
knitr::include_graphics('TopLigands.png')

## First, select ligands upregulated in source population
source.marker.table <- Seurat::FindMarkers(ECstory.sc, ident.1 = "Endothelial", 
                                           logfc.threshold = 1.5, 
                                           only.pos = TRUE)
source.marker.genes <- rownames(source.marker.table)

path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
edge.score.file <- paste0(file.lab, "_all_ligand_receptor_network_edges.csv")
source.population <- "Endothelial"
target.populations <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")

source.population <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
target.populations <- "Endothelial"
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")


source.population <- c("Endothelial")
target.populations <- c("Fibroblasts", "SMC")
population.cols <- c("#0053c8", "#fb8500", "#c10023")

source.population <- "Endothelial"
target.populations <- c("Lymphatic Endo","Mono/Macs")
population.cols <- c("#aa00ff", "#5f777f","#0053c8")

source.population <- c("Endothelial")
target.populations <- c("Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f")

NetworkTreePlot(path.file = path.file, 
                edge.score.file = edge.score.file,
                source.population = source.population,
                target.populations = target.populations,
                source.marker.genes = source.marker.genes,
                population.cols = population.cols)

knitr::include_graphics('tree_plot.png')

#----------
file.lab <- "TIP4"

ECstory.sc.PHPF <- subset(ECstory.sc, samplecondition %in% c("IPF"))

GenerateEdgeWeights(seurat.object = ECstory.sc.PHPF,
                    file.label = file.lab,
                    species = "human",
                    verbose=T)


edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)


GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)

paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)


results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4,
                                     return.results = TRUE)
head(results.table)

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")
col.set <- polychrome_pal
col.set <- c("#c10023","#0099cc","#008e17", "#fde800","#aa00ff","#5f777f","#f2a287")
col.set <- c("#00896e", "#5f777f", "#00cc99", "#f2a287", "#c10023", "#4ffc00","#fb8500", "#0053c8", "#0099cc", "#aa00ff", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

#p.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO","EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc", "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")
#00896e
p.labels <- unique(ECstory.sc$clusters_tommy.sc)

names(col.set) <- p.labels

## Generate a circle plot
cell.network.file = paste0("Permutation_tests_", file.lab, "_network.csv")
CellCirclePlot(input.file = cell.network.file, 
               adj_pval_thresh = 0.01, 
               col.set = col.set)

knitr::include_graphics('TIP_small_circle_plot.png')

## Generate an inbound vs outbound plot
InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05)
knitr::include_graphics('TIP_small_inbound_outbound.png')

InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05,
                    lab.weight.thresh = 10000)

## Plot top ligands for a cell population
path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
cell.type <- "Endothelial"
PlotTopLigands(input.file = path.file, 
               col.use = "#0053c8",
               cell.identity = cell.type)
knitr::include_graphics('TopLigands.png')

## First, select ligands upregulated in source population
source.marker.table <- Seurat::FindMarkers(ECstory.sc, ident.1 = "Endothelial", 
                                           logfc.threshold = 1.5, 
                                           only.pos = TRUE)
source.marker.genes <- rownames(source.marker.table)

path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
edge.score.file <- paste0(file.lab, "_all_ligand_receptor_network_edges.csv")
source.population <- "Endothelial"
target.populations <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")

source.population <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
target.populations <- "Endothelial"
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")


source.population <- c("Endothelial")
target.populations <- c("Fibroblasts", "SMC")
population.cols <- c("#0053c8", "#fb8500", "#c10023")

source.population <- "Endothelial"
target.populations <- c("Lymphatic Endo","Mono/Macs")
population.cols <- c("#aa00ff", "#5f777f","#0053c8")

source.population <- c("Endothelial")
target.populations <- c("Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f")

NetworkTreePlot(path.file = path.file, 
                edge.score.file = edge.score.file,
                source.population = source.population,
                target.populations = target.populations,
                source.marker.genes = source.marker.genes,
                population.cols = population.cols)

knitr::include_graphics('tree_plot.png')

#---------------

file.lab <- "TIP5"

ECstory.sc.PAH <- subset(ECstory.sc, samplecondition %in% c("IPAH"))

GenerateEdgeWeights(seurat.object = ECstory.sc.PAH,
                    file.label = file.lab,
                    species = "human",
                    verbose=T)


edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)


GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)

paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)


results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4,
                                     return.results = TRUE)
head(results.table)

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")
col.set <- polychrome_pal
col.set <- c("#c10023","#0099cc","#008e17", "#fde800","#aa00ff","#5f777f","#f2a287")
col.set <- c("#00896e", "#5f777f", "#00cc99", "#f2a287", "#c10023", "#4ffc00","#fb8500", "#0053c8", "#0099cc", "#aa00ff", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

#p.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO","EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc", "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")
#00896e
p.labels <- unique(ECstory.sc$clusters_tommy.sc)

names(col.set) <- p.labels

## Generate a circle plot
cell.network.file = paste0("Permutation_tests_", file.lab, "_network.csv")
CellCirclePlot(input.file = cell.network.file, 
               adj_pval_thresh = 0.01, 
               col.set = col.set)

knitr::include_graphics('TIP_small_circle_plot.png')

## Generate an inbound vs outbound plot
InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05)
knitr::include_graphics('TIP_small_inbound_outbound.png')

InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05,
                    lab.weight.thresh = 10000)

## Plot top ligands for a cell population
path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
cell.type <- "Endothelial"
PlotTopLigands(input.file = path.file, 
               col.use = "#0053c8",
               cell.identity = cell.type)
knitr::include_graphics('TopLigands.png')

## First, select ligands upregulated in source population
source.marker.table <- Seurat::FindMarkers(ECstory.sc, ident.1 = "Endothelial", 
                                           logfc.threshold = 1.5, 
                                           only.pos = TRUE)
source.marker.genes <- rownames(source.marker.table)

path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
edge.score.file <- paste0(file.lab, "_all_ligand_receptor_network_edges.csv")
source.population <- "Endothelial"
target.populations <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")

source.population <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
target.populations <- "Endothelial"
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")


source.population <- c("Endothelial")
target.populations <- c("Fibroblasts", "SMC")
population.cols <- c("#0053c8", "#fb8500", "#c10023")

source.population <- "Endothelial"
target.populations <- c("Lymphatic Endo","Mono/Macs")
population.cols <- c("#aa00ff", "#5f777f","#0053c8")

source.population <- c("Endothelial")
target.populations <- c("Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f")

NetworkTreePlot(path.file = path.file, 
                edge.score.file = edge.score.file,
                source.population = source.population,
                target.populations = target.populations,
                source.marker.genes = source.marker.genes,
                population.cols = population.cols)

knitr::include_graphics('tree_plot.png')


#-----------
ECstory.subset.fil.sc$namedclusters2 <- ""
ECstory.subset.fil.sc$namedclusters2[ECstory.subset.fil.sc$namedclusters == "Immune modulatory"] <- "Immunemodulatory"
ECstory.subset.fil.sc$namedclusters2[ECstory.subset.fil.sc$namedclusters == "Homeostatic"] <- "Homeostatic"
ECstory.subset.fil.sc$namedclusters2[ECstory.subset.fil.sc$namedclusters == "Tone modulatory"] <- "Defense"

ECstory.subset.fil$namedclusters2 <- ""
ECstory.subset.fil$namedclusters2[ECstory.subset.fil$namedclusters == "Immune modulatory"] <- "Defense"
ECstory.subset.fil$namedclusters2[ECstory.subset.fil$namedclusters == "Homeostatic"] <- "Homeostatic"
ECstory.subset.fil$namedclusters2[ECstory.subset.fil$namedclusters == "Tone modulatory"] <- "Immunomodulatory"

file.lab <- "TIP6"

ECstory.subset.fil.sc.Donor <- subset(ECstory.subset.fil.sc, samplecondition %in% c("Donor"))

unique(Idents(ECstory.subset.fil.sc.Donor))


Idents(ECstory.subset.fil.sc.Donor) <- ECstory.subset.fil.sc.Donor$namedclusters2

GenerateEdgeWeights(seurat.object = ECstory.subset.fil.sc.Donor,
                    file.label = file.lab,
                    species = "human",
                    verbose=T)


edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)


GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)

paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)


results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4,
                                     return.results = TRUE)
head(results.table)

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")
col.set <- polychrome_pal
col.set <- c("#c10023","#0099cc","#008e17", "#fde800","#aa00ff","#5f777f","#f2a287")
col.set <- c("#00896e", "#5f777f", "#00cc99", "#f2a287", "#c10023", "#4ffc00","#fb8500", "#0053c8", "#0099cc", "#aa00ff", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

col.set <- c("#0053c8", "#008e17","#c10023")


#p.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO","EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc", "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")
#00896e
p.labels <- unique(ECstory.subset.fil.sc.Donor$namedclusters2)

names(col.set) <- p.labels

## Generate a circle plot
cell.network.file = paste0("Permutation_tests_", file.lab, "_network.csv")
CellCirclePlot(input.file = cell.network.file, 
               adj_pval_thresh = 0.05, 
               col.set = col.set)

knitr::include_graphics('TIP_small_circle_plot.png')

## Generate an inbound vs outbound plot
InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05)
knitr::include_graphics('TIP_small_inbound_outbound.png')

InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05,
                    lab.weight.thresh = 10000)

## Plot top ligands for a cell population
path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
cell.type <- "Endothelial"
PlotTopLigands(input.file = path.file, 
               col.use = "#0053c8",
               cell.identity = cell.type)
knitr::include_graphics('TopLigands.png')

## First, select ligands upregulated in source population
source.marker.table <- Seurat::FindMarkers(ECstory.sc, ident.1 = "Endothelial", 
                                           logfc.threshold = 1.5, 
                                           only.pos = TRUE)
source.marker.genes <- rownames(source.marker.table)

path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
edge.score.file <- paste0(file.lab, "_all_ligand_receptor_network_edges.csv")
source.population <- "Endothelial"
target.populations <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")

source.population <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
target.populations <- "Endothelial"
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")


source.population <- c("Endothelial")
target.populations <- c("Fibroblasts", "SMC")
population.cols <- c("#0053c8", "#fb8500", "#c10023")

source.population <- "Endothelial"
target.populations <- c("Lymphatic Endo","Mono/Macs")
population.cols <- c("#aa00ff", "#5f777f","#0053c8")

source.population <- c("Endothelial")
target.populations <- c("Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f")

NetworkTreePlot(path.file = path.file, 
                edge.score.file = edge.score.file,
                source.population = source.population,
                target.populations = target.populations,
                source.marker.genes = source.marker.genes,
                population.cols = population.cols)

knitr::include_graphics('tree_plot.png')

#-----------------

file.lab <- "TIP7"

ECstory.subset.fil.sc.PHPF <- subset(ECstory.subset.fil.sc, samplecondition %in% c("IPF"))

unique(Idents(ECstory.subset.fil.sc.PHPF))


Idents(ECstory.subset.fil.sc.PHPF) <- ECstory.subset.fil.sc.PHPF$namedclusters2

GenerateEdgeWeights(seurat.object = ECstory.subset.fil.sc.PHPF,
                    file.label = file.lab,
                    species = "human",
                    verbose=T)


edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)


GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)

paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)


results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4,
                                     return.results = TRUE)
head(results.table)

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")
col.set <- polychrome_pal
col.set <- c("#c10023","#0099cc","#008e17", "#fde800","#aa00ff","#5f777f","#f2a287")
col.set <- c("#00896e", "#5f777f", "#00cc99", "#f2a287", "#c10023", "#4ffc00","#fb8500", "#0053c8", "#0099cc", "#aa00ff", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

col.set <- c("#0053c8", "#008e17","#c10023")


#p.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO","EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc", "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")
#00896e
p.labels <- unique(ECstory.subset.fil.sc.Donor$namedclusters2)

names(col.set) <- p.labels

## Generate a circle plot
cell.network.file = paste0("Permutation_tests_", file.lab, "_network.csv")
CellCirclePlot(input.file = cell.network.file, 
               adj_pval_thresh = 0.05, 
               col.set = col.set)

knitr::include_graphics('TIP_small_circle_plot.png')

## Generate an inbound vs outbound plot
InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05)
knitr::include_graphics('TIP_small_inbound_outbound.png')

InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05,
                    lab.weight.thresh = 10000)

## Plot top ligands for a cell population
path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
cell.type <- "Endothelial"
PlotTopLigands(input.file = path.file, 
               col.use = "#0053c8",
               cell.identity = cell.type)
knitr::include_graphics('TopLigands.png')

## First, select ligands upregulated in source population
source.marker.table <- Seurat::FindMarkers(ECstory.sc, ident.1 = "Endothelial", 
                                           logfc.threshold = 1.5, 
                                           only.pos = TRUE)
source.marker.genes <- rownames(source.marker.table)

path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
edge.score.file <- paste0(file.lab, "_all_ligand_receptor_network_edges.csv")
source.population <- "Endothelial"
target.populations <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")

source.population <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
target.populations <- "Endothelial"
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")


source.population <- c("Endothelial")
target.populations <- c("Fibroblasts", "SMC")
population.cols <- c("#0053c8", "#fb8500", "#c10023")

source.population <- "Endothelial"
target.populations <- c("Lymphatic Endo","Mono/Macs")
population.cols <- c("#aa00ff", "#5f777f","#0053c8")

source.population <- c("Endothelial")
target.populations <- c("Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f")

NetworkTreePlot(path.file = path.file, 
                edge.score.file = edge.score.file,
                source.population = source.population,
                target.populations = target.populations,
                source.marker.genes = source.marker.genes,
                population.cols = population.cols)

knitr::include_graphics('tree_plot.png')


#---------------

file.lab <- "TIP8"

ECstory.subset.fil.sc.PAH <- subset(ECstory.subset.fil.sc, samplecondition %in% c("IPAH"))

unique(Idents(ECstory.subset.fil.sc.PAH))


Idents(ECstory.subset.fil.sc.PAH) <- ECstory.subset.fil.sc.PAH$namedclusters2

GenerateEdgeWeights(seurat.object = ECstory.subset.fil.sc.PAH,
                    file.label = file.lab,
                    species = "human",
                    verbose=T)


edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)


GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)

paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)


results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4,
                                     return.results = TRUE)
head(results.table)

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")
col.set <- polychrome_pal
col.set <- c("#c10023","#0099cc","#008e17", "#fde800","#aa00ff","#5f777f","#f2a287")
col.set <- c("#00896e", "#5f777f", "#00cc99", "#f2a287", "#c10023", "#4ffc00","#fb8500", "#0053c8", "#0099cc", "#aa00ff", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8","#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

col.set <- c("#0053c8", "#008e17","#c10023")


#p.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO","EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc", "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")
#00896e
p.labels <- unique(ECstory.subset.fil.sc.Donor$namedclusters2)

names(col.set) <- p.labels

## Generate a circle plot
cell.network.file = paste0("Permutation_tests_", file.lab, "_network.csv")
CellCirclePlot(input.file = cell.network.file, 
               adj_pval_thresh = 0.05, 
               col.set = col.set)

knitr::include_graphics('TIP_small_circle_plot.png')

## Generate an inbound vs outbound plot
InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05)
knitr::include_graphics('TIP_small_inbound_outbound.png')

InboundOutboundPlot(input.file = cell.network.file, 
                    col.set = col.set,
                    p.adj.thresh = 0.05,
                    lab.weight.thresh = 10000)

## Plot top ligands for a cell population
path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
cell.type <- "Homeostatic"
PlotTopLigands(input.file = path.file, 
               col.use = "#0053c8",
               cell.identity = cell.type)


cell.type <- "Immune modulatory"
PlotTopLigands(input.file = path.file, 
               col.use = "#c10023",
               cell.identity = cell.type)

cell.type <- "Defense"
PlotTopLigands(input.file = path.file, 
               col.use = "#008e17",
               cell.identity = cell.type)

knitr::include_graphics('TopLigands.png')

## First, select ligands upregulated in source population
source.marker.table <- Seurat::FindMarkers(ECstory.sc, ident.1 = "Endothelial", 
                                           logfc.threshold = 1.5, 
                                           only.pos = TRUE)
source.marker.genes <- rownames(source.marker.table)

path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
edge.score.file <- paste0(file.lab, "_all_ligand_receptor_network_edges.csv")
source.population <- "Endothelial"
target.populations <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")

source.population <- c("Endothelial", "Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
target.populations <- "Endothelial"
population.cols <- c("#aa00ff", "#fb8500", "#c10023", "#5f777f","#0053c8")


source.population <- c("Endothelial")
target.populations <- c("Fibroblasts", "SMC")
population.cols <- c("#0053c8", "#fb8500", "#c10023")

source.population <- "Endothelial"
target.populations <- c("Lymphatic Endo","Mono/Macs")
population.cols <- c("#aa00ff", "#5f777f","#0053c8")

source.population <- c("Endothelial")
target.populations <- c("Lymphatic Endo", "Fibroblasts", "SMC", "Mono/Macs")
population.cols <- c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f")

source.population <- c("Homeostatic")
target.populations <- c("Homeostatic", "Immune modulatory", "Defense")
population.cols <- c("#0053c8", "#008e17","#c10023")


"Immune modulatory"
ECstory.subset.fil.sc$namedclusters2[ECstory.subset.fil.sc$namedclusters == "Homeostatic"] <- "Homeostatic"
ECstory.subset.fil.sc$namedclusters2[ECstory.subset.fil.sc$namedclusters == "Tone modulatory"] <- "Defense"


NetworkTreePlot(path.file = path.file, 
                edge.score.file = edge.score.file,
                source.population = source.population,
                target.populations = target.populations,
                source.marker.genes = source.marker.genes,
                population.cols = population.cols)

knitr::include_graphics('tree_plot.png')


#ECstorysupp
DimPlot(ECstory, group.by = "clusters_tommy5", cols=polychrome_pal)

ECstory$majorpop <- ""
ECstory$majorpop[(ECstory$clusters_tommy5 == "B-cells" )] <- "Immune cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "Dendritic cells")] <- "Immune cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "Endo1")] <- "Endothelial cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "Endo2")] <- "Endothelial cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "Epithelial Cells")] <- "Epithelial"
ECstory$majorpop[(ECstory$clusters_tommy5 == "Erythrocytes" )] <- "Erythrocytes"
ECstory$majorpop[(ECstory$clusters_tommy5 == "Fibroblasts" )] <- "Mesenchymal cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "Granulocytes" )] <- "Immune cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "Lymphatic Endo" )] <- "Endothelial cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "Mono/Macs" )] <- "Immune cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "NK-cells" )] <- "Immune cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "SMC1" )] <- "Mesenchymal cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "SMC2" )] <- "Mesenchymal cells"
ECstory$majorpop[(ECstory$clusters_tommy5 == "T-cells" )] <- "Immune cells"

levelsmajorpop <- c("Immune cells","Endothelial cells", "Mesenchymal cells", "Epithelial","Erythrocytes")
ECstory$majorpop <- factor(ECstory$majorpop, levels= levelsmajorpop)

DimPlot(ECstory, group.by = "majorpop", cols=c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f"))
DimPlot(ECstory, group.by = "majorpop", cols=c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f"), split.by = "samplecondition")

ECstoryDonor <- subset(ECstory, samplecondition %in% c("Donor"))
ECstoryPAH <- subset(ECstory, samplecondition %in% c("IPAH"))
ECstoryPHPF <- subset(ECstory, samplecondition %in% c("IPF"))

DimPlot(ECstoryDonor, group.by = "majorpop", cols=c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f"), pt.size = 0.05)
DimPlot(ECstoryPAH, group.by = "majorpop", cols=c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f"), pt.size = 0.05)
DimPlot(ECstoryPHPF, group.by = "majorpop", cols=c("#0053c8", "#aa00ff", "#fb8500", "#c10023", "#5f777f"), pt.size = 0.05)



#venn
library(VennDiagram)

sharedligands <- c("SPP1",	"PTPRM",	"PPIA",	"POSTN", 	"PECAM1",	"MIF",	"JAM2"	"IGFBP3",	"ICAM2"	, "HLA-E	HLA-DRB5	HLA-DRB1	HLA-DRA	HLA-DQB1	HLA-DPB1	HLA-DPA1	HLA-DMA	GJA1	FN1	EFNB2	COL4A2	CHOLESTEROL -LIPA	CD99	CD55	CD34	CCL14	BMP6	APP	ANXA1	ITGAV+ITGB1
)
