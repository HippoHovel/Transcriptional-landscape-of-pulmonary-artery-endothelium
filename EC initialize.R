set.seed(1)
library(hrbrthemes)
library(RColorBrewer)
library(Seurat)
library(cowplot)
library(ggplot2)
library(scran)
library(SeuratData)
library(patchwork)
library(Matrix)
library(dplyr)
library(umap)
library(reticulate)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(SeuratDisk)
library(SeuratData)
library(ggsignif)
library(BiocManager)
library(scater)
library(patchwork)
library(cowplot)
library(viridis)
library(tidyverse)
library(reshape2)
library(parallel)
library(DoubletFinder)
library(scCustomize)
library(devtools)
library(enrichR)
library(reshape2)
library(plotly)
library(Scillus)
library(circlize)
library(STRINGdb)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(monocle3)
library(SeuratWrappers)
library("scProportionTest")

#define function for excluding certain identifyers when subsetting
'%!in%' <- function(x,y)!('%in%'(x,y))

#initialize some colour palettes 
polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome",shuffle_pal = T)
varibow_pal <- DiscretePalette_scCustomize(num_colors = 50, palette = "varibow", shuffle_pal = T)
mypastelcol = c("#C881AE","#81C881","#818AC8", "#C7C881", "#cb90f4", "#86cdcd", "#a0a0a0", "#cfad64")
myphasecol = c("#804daa", "#c7c881", "#323227")


ECstory <- readRDS("~/ECstory_publication.rds")

#QC cut
# ECstory$QCcut <- ""
# ECstory$QCcut <- "Viable"
# ECstory$QCcut[ECstory$nFeature_RNA > 4000] <- "Multiplett"
# ECstory$QCcut[ECstory$nFeature_RNA < 200] <- "Cell remnant"
# ECstory$QCcut[ECstory$percent.mt > 10] <- "low quality read"
DimPlot(ECstory, group.by = "QCcut", split.by = "samplecondition")
ECstory.fil <- subset(ECstory, QCcut %in% c("Viable")) #skip this step for custom filtering

#newer seurat versions require joining of layers 
ECstory.j <- JoinLayers(ECstory.fil)








