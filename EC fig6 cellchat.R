library(CellChat)

ECstory.fil$clusters3 <- ""
ECstory.fil$clusters4 <- ECstory.subset$seurat_clusters
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "27")] <- "Endo MT"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "0" )] <- "Fibroblasts"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "24")] <- "Fibroblasts"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "4" )] <- "SMC"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "14")] <- "SMC"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "16")] <- "SMC" #pericytes seem to be a subset of these
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "9" )] <- "SMC"    
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "21")] <- "B-cells"  
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "23")] <- "Lymphatic EC" 
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "20")] <- "Granulocytes"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "22")] <- "Erythrocytes"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "13")] <- "Mono/Macs"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "7" )] <- "Mono/Macs"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "2" )] <- "Mono/Macs"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "25")] <- "Dendritic cells"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "19")] <- "Dendritic cells"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "18")] <- "Dendritic cells"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "11")] <- "Mono/Macs"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "5")] <- "Mono/Macs"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "17")] <- "Mono/Macs"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "15")] <- "Epithelial cells"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "26")] <- "Epithelial cells"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "8")] <- "NK cells"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "1")] <- "T-cells"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "3")] <- "T-cells"
ECstory.fil$clusters3[(ECstory.fil$seurat_clusters == "10")] <- "T-cells"
ECstory.fil$clusters3[(ECstory.fil$clusters4 == "0")] <- "EC0"
ECstory.fil$clusters3[(ECstory.fil$clusters4 == "1")] <- "EC1"
ECstory.fil$clusters3[(ECstory.fil$clusters4 == "2")] <- "EC2"
ECstory.fil$clusters3[(ECstory.fil$clusters4 == "3")] <- "EC3"
ECstory.fil$clusters3[(ECstory.fil$clusters4 == "4")] <- "EC4"


DimPlot(ECstory.fil, group.by = "clusters3", cols=polychrome_pal)

DimPlot(ECstory.fil.subset.fil)
DimPlot(ECstory.fil.subset)

DimPlot(test, group.by = "clusters3")

DimPlot(ECstory.fil, group.by = "clusters2")
ECstory.cellchat <- subset(ECstory.fil, QCcut %in% c("Viable"))


Donor <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("Donor")) & (clusters3 %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))
PH <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("IPAH", "IPF")) & (clusters3 %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))


#init cell chat
data.input <- GetAssayData(Donor, assay = "RNA", slot = "data") # normalized data matrix #switchdataset to PH or Donor
cellchat <- createCellChat(object = Donor, group.by = "clusters3", assay = "RNA") #switchdataset to PH or Donor
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


pal3 <- c("#7ED7D1FF", "#683B79FF", "#F8A19FFF", "#16FF32FF", "#FA0087FF", "#325A9BFF", "#2ED9FFFF","#f39d10","#038637", "#F6222EFF","#FBE426FF","#f39d10")

names(pal3) <- unique(Donor$clusters3)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count,top=0.3,color.use =pal3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=0.3,color.use =pal3,vertex.weight = groupSize, weight.scale = T, title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)

netVisual_circle(cellchat@net$count,top=0.5,color.use =pal3,sources.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=0.5,color.use =pal3,sources.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)

netVisual_circle(cellchat@net$count,top=1,color.use =pal3,targets.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=1,color.use =pal3,targets.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)

#run only one of these commented lines depending on the current loaded set
#cellchat_PH_allcells <- cellchat
#cellchat_Donor_allcells <- cellchat

#continue after computing both parts of the merged object: cellchat_PH_allcells and cellchat_Donor_allcells


object.list <- list(Donor = cellchat_Donor_allcells, PH = cellchat_PH_allcells)
cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

cellchat.merge@meta$clusters3
unique(cellchat.merge@meta$clusters3)

levels(cellchat.merge@meta$clusters3)
cellchat.merge@meta$clusters3 <- factor(cellchat.merge@meta$clusters3, 
                                              levels = c("EC0", "EC1", "EC2",
                                                         "Fibroblasts", "SMC",
                                                         "Mono/Macs", "NK cells", "T-cells", "B-cells", "Granulocytes", "Dendritic cells"))

levels(cellchat.merge@meta$clusters3)


names(pal3) <- unique(cellchat.merge@meta$clusters3)

par(mfrow=c(1,1))
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("count"), color.use = pal3)
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("weight"), color.use = pal3)
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("count"))
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("weight"))



pal3 <- c("#7ED7D1FF", "#683B79FF", "#F6222EFF", "#325A9BFF", "#038637", "#16FF32FF", "#2ED9FFFF","#f39d10","#FA0087FF", "#F8A19FFF","#FBE426FF","#f39d10")


par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat.merge, weight.scale = T, color.use = pal3)

par(mfrow = c(1,2), xpd=TRUE)
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










#Figure8


Donor2 <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("Donor")) & (cluster.identity %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))
PAH <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("IPAH")) & (cluster.identity %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))
PF <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("IPF")) & (cluster.identity %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))

PH <- JoinLayers(subset(ECstory.cellchat, (samplecondition %in% c("IPAH", "IPF")) & (cluster.identity %!in% c("Platelets", "Endo MT", "Erythrocytes", "Epithelial cells", "EC3", "EC4"))))

#init cell chat
data.input <- GetAssayData(PF, assay = "RNA", slot = "data") # normalized data matrix #switchdataset to PH or Donor
cellchat <- createCellChat(object = PF, group.by = "cluster.identity", assay = "RNA") #switchdataset to PH or Donor
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


#cellchatDonor2 <- cellchat
#cellchatPAH <- cellchat
#cellchatPF <- cellchat

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count,top=0.3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=0.3,vertex.weight = groupSize, weight.scale = T, title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)

object.list <- list(Donor = cellchatDonor2, PAH = cellchatPAH) 
object.list <- list(Donor = cellchatDonor2, PAH = cellchatPF)
cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)


netVisual_bubble(cellchat.merge, sources.use = 3, targets.use = NULL ,  comparison = c(1, 2), title.name = "Increased signaling in LS", 
                 angle.x = 45, remove.isolate = F,color.text=c("black","#0491b3"), font.size = 8,color.heatmap = c("Spectral"), max.dataset = NULL )


netVisual_bubble(cellchat.merge, sources.use = NULL, targets.use = 3,  comparison = c(1, 2), max.dataset = NULL, title.name = "Increased signaling in LS", 
                 angle.x = 45, remove.isolate = F,color.text=c("black","#0491b3"), font.size = 8)







#PHPF
enrichrinput <- c("PLXNA2", "NLGN1", "ITGA6", "ITGB4", "ITGA6", "ITGB1", "IGF1R", "HRH1", "ADGRL2", "CD66")#receptors PHPF

enrichrinput2 <- c("TNFRSF1A", "CD47", "NRP1", "PLXNA2", "ADGRG6", "NCL", "CXCR4", "ITGA10", "ITGB1", "CDH5", "CD99", "CALCRL")#receptors PAH

#shared
enrichrinput <- c("SPP1", "PTPRM", "PPIA", "POSTN", "PECAM1", "MIF", "JAM2", "IGFBP3", "ICAM2", "HLA-E", "HLA-DRB5", "HLA-DRB1",
                  "HLA-DRA", "HLA-DQB1", "HLA-DPB1", "HLA-DPA1", "HLA-DMA", "GJA1", "FN1", "EFNB2", "COL4A2", "LIPA", "CD99",
                  "CD55", "CD34", "CCL14", "BMP6", "APP", "ANXA1", "ITGAV", "ITGB1")#shared ligands

enrichrinput2 <- c("VEGFR1", "ITGAV", "ITGB1", "PLXND1", "PTPRM", "BSG", "PECAM1", "ITGA5", "P4HB", "JAM2", "PARD3", "GJA1", "EPHA4", "ACKR1", "RORA", "SELP", "CD74")#shared receptors


#PAH
enrichrinput <- c("TGFB1", "NRG1", "MDK", "LAMB2", "HLA-DQA1", "GRN", "GAS6", "CDH5", "APOE")#ligand PAH

enrichrinput2 <- c("SEMA6A", "PGE2-PTGES3", "NAMPT", "LAMC1", "LAMB1", "ICAM1", 
                   "HLA-F", "HLA-C", "HLA-B", "HLA-A", "CXCL2", "CXCL12", 
                   "COL4A1", "CCL2")#ligand PH-PF


#enrichr
"#e06606","#0491b3"

enrichrpicks <- mergedPH$gene[mergedPH$colour == "#0491b3"]
enrichrpicks2 <- mergedPH$gene[mergedPH$colour == "#e06606"]


enriched <- enrichr(enrichrpicks, dbs) #custom.gene.list.up is a simple vector with up regulated gene names 

enriched$GO_Biological_Process_2023$Count <- str_extract(enriched$GO_Biological_Process_2023$Overlap, ".+?(?=/)")
enriched$GO_Biological_Process_2023$Count <- as.numeric(enriched$GO_Biological_Process_2023$Count)

enrichedplot <- head(enriched$GO_Biological_Process_2023, 100)
enrichedplot2 <- head(enrichedplot,  n=50) 

enriched2 <- enrichr(enrichrpicks2, dbs) #custom.gene.list.down is a simple vector with down regulated gene names

enriched2$GO_Biological_Process_2023$Count <- str_extract(enriched2$GO_Biological_Process_2023$Overlap, ".+?(?=/)")
enriched2$GO_Biological_Process_2023$Count <- as.numeric(enriched2$GO_Biological_Process_2023$Count)

enrichedplot3 <- head(enriched2$GO_Biological_Process_2023, 5000)

enrichedplot4 <- head(enrichedplot3,  n=50)
#enrichedplot4$Count <- enrichedplot4$Count *-1 #to visualize down regulation




enrichedplot4 <- enrichedplot4[order(enrichedplot4$Adjusted.P.value),]

enrichedplot2 <- enrichedplot2[order(enrichedplot2$Adjusted.P.value),]



enrichedplot2 <- enrichedplot2 %>%
  mutate(
    Overlap_Count = as.numeric(str_extract(Overlap, "^[0-9]+")),
    Overlap_Size = as.numeric(str_extract(Overlap, "[0-9]+$")),
    Overlap_Percentage = (Overlap_Count / Overlap_Size) * 100
  )




enrichedplot4 <- enrichedplot4 %>%
  mutate(
    Overlap_Count = as.numeric(str_extract(Overlap, "^[0-9]+")),
    Overlap_Size = as.numeric(str_extract(Overlap, "[0-9]+$")),
    Overlap_Percentage = (Overlap_Count / Overlap_Size) * 100
  )



ggplot() +
  geom_bar(data=enrichedplot4, stat='identity',
           aes(x=reorder(Term, rev(P.value)), y=Overlap_Percentage),
           fill=("#e06606"))+
  geom_text(data=enrichedplot4,
            aes(x=reorder(Term, rev(P.value)), y=0.1, label = format(Adjusted.P.value,scientific=T,digits=2)),
            hjust = 0,vjust = 0,
            colour = "#404040",
            size=3)+
  #geom_text(data=enrichedplot4,
  #aes(x=reorder(Term, rev(Adjusted.P.value)), y=0.1, label =Genes),
  #hjust = 0,vjust = 0.5, 
  #colour = "black", 
  #size=2.8)+
  
  geom_bar(data=enrichedplot2, 
           stat='identity',
           aes(x=reorder(Term, rev(P.value)), 
               y=Overlap_Percentage),
           #fill=factor(ifelse(Term =="Negative Regulation Of Cellular Process(GO:0048523)","purple","red"))))+ 
           #ifelse(enrichedplot2$Term=="Negative Regulation Of Cellular Process (GO:0048523)",fill="yellow", fill="red"))+
           #+
           fill=("#0491b3"))+
  #geom_text(data=enrichedplot2,
  #aes(x=reorder(Term, rev(Adjusted.P.value)), y=-0.1, label = Genes),
  #hjust = 1,vjust = 0.5, 
  #colour = "black", 
  #size=2.8)+
  geom_text(data=enrichedplot2,
            aes(x=reorder(Term, rev(P.value)), y=0.1, label = format(Adjusted.P.value,scientific=T,digits=2)),
            hjust = 0,vjust = 0, 
            colour = "#404040", 
            size=3)+
  
  #ylim(-5,8)+ #change to form fit your data
  theme_test()+
  #theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+#vanish standard axis names
  coord_flip()

