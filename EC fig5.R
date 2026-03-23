
ECstory.subset.fil$hvd <- ""
ECstory.subset.fil$hvd[ECstory.subset.fil$samplecondition == "Donor"] <- "Donor"
ECstory.subset.fil$hvd[ECstory.subset.fil$samplecondition == "IPAH"] <- "PH"
ECstory.subset.fil$hvd[ECstory.subset.fil$samplecondition == "IPF"] <- "PH"

DimPlot(ECstory.subset.fil, group.by = "samplecondition", cols= c("#17a4ad","#ad179d","#97ab17") )+
DimPlot(ECstory.subset.fil,group.by = "hvd", cols= c("#17a4ad","#ad179d") )


hvd2.deg <- FindAllMarkers(ECstory.subset.fil)
hvd2.deg  <- hvd2.deg[order(hvd2.deg $avg_log2FC),]
top <- hvd2.deg

#top <- rownames_to_column(top, var="gene")

top$avg_log2FC <- as.numeric(top$avg_log2FC)

top$colour = ifelse((top$avg_log2FC>= 0.5) & (top$p_val_adj <= 0.05 & top$p_val_adj >=0), "#ad179d", 
                    ifelse((top$avg_log2FC <= -0.5) & (top$p_val_adj <= 0.05), "#17a4ad", "gray"))

top$rank <- abs(top$avg_log2FC * -log10(top$p_val_adj))

top$labels = ifelse(top$gene %in% head(top[order(top$rank, decreasing = T),], n = 100)$gene, ifelse(top$colour != "gray",top$gene, ""), "")  
top$labels = ifelse(top$gene %in% head(top[order(top$rank, decreasing = T),], n = 30)$gene, ifelse(top$colour != "gray",top$gene, ""), "")  




top20 <- top
top20$valid = ""
top20$valid [(top20$p_val_adj < 0.05)] <- TRUE
top20$valid[(top20$pct.1 < 0.2)] <- FALSE
top20$valid[(top20$pct.2 > 0.2)] <- FALSE
top20 <- subset(top20, subset = top20$valid == TRUE )

top20 <- top20 %>% group_by(cluster) %>% top_n(300, avg_log2FC)
top20 <- top20[order(top20$cluster),]
top20 <- na.omit(top20)

top20$avg_log2FC[top20$cluster == "Donor"] <- top20$avg_log2FC *-1
top20$colour[top20$cluster == "Donor" & top20$colour != "gray"] <- "#17a4ad"
top20$labels = ifelse(top20$gene %in% head(top20[order(top20$rank, decreasing = T),], n = 30)$gene, ifelse(top20$colour != "gray",top20$gene, ""), "")  

top20 <- subset(top20, colour %!in% c("gray"))

ggplot(top20, aes(x = top20$avg_log2FC, y = -log10(top20$p_val_adj)))+
  geom_text_repel(label = top20$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_point(color = top20$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")


#Heatmap code
top20 <- top
top20$valid = ""
top20$valid [(top20$p_val_adj < 0.05)] <- TRUE
top20$valid[(top20$pct.1 < 0.2)] <- FALSE
top20$valid[(top20$pct.2 > 0.2)] <- FALSE
top20 <- subset(top20, subset = top20$valid == TRUE )

top20 <- top20 %>% group_by(cluster) %>% top_n(15, avg_log2FC)
top20 <- top20[order(top20$cluster),]
top20 <- na.omit(top20)





ECstory.subset.fil <- ScaleData(ECstory.subset.fil, features=top20$gene)

plot_heatmap(ECstory.subset.fil,
             markers = top20$gene,
             sort_var = c("hvd","seurat_clusters"),
             anno_var = c("hvd", "seurat_clusters"),
             anno_colors= list(c("#17a4ad","#ad179d"), c("#b20505","#0312ad","#03ad46")),
             row_font_size = 6.5)







#seperate genes for clusters
#red
red_cluster <- subset(ECstory.subset.fil, seurat_clusters %in% "0")


unique(Idents(red_cluster))
Idents(red_cluster) <- red_cluster$hvd
red_cluster.deg <- FindAllMarkers(red_cluster)


red_cluster.deg  <- red_cluster.deg[order(red_cluster.deg $avg_log2FC),]
top.red <- red_cluster.deg

#blue
blue_cluster <- subset(ECstory.subset.fil, seurat_clusters %in% "1")


Idents(blue_cluster) <- blue_cluster$hvd
blue_cluster.deg <- FindAllMarkers(blue_cluster)


blue_cluster.deg  <- blue_cluster.deg[order(blue_cluster.deg $avg_log2FC),]
top.blue <- blue_cluster.deg

#green
green_cluster <- subset(ECstory.subset.fil, seurat_clusters %in% "2")


Idents(green_cluster) <- green_cluster$hvd
green_cluster.deg <- FindAllMarkers(green_cluster)


green_cluster.deg  <- green_cluster.deg[order(green_cluster.deg $avg_log2FC),]
top.green <- green_cluster.deg


#top <- rownames_to_column(top, var="gene")

top <- top.red #only run 1 of these three lines depending on current investigated
top <- top.blue
top <- top.green

top$avg_log2FC <- as.numeric(top$avg_log2FC)

top$colour = ifelse((top$avg_log2FC>= 0.5) & (top$p_val_adj <= 0.05 & top$p_val_adj >=0), "#ad179d", 
                    ifelse((top$avg_log2FC <= -0.5) & (top$p_val_adj <= 0.05), "#17a4ad", "gray"))

top$rank <- abs(top$avg_log2FC * -log10(top$p_val_adj))

top$labels = ifelse(top$gene %in% head(top[order(top$rank, decreasing = T),], n = 100)$gene, ifelse(top$colour != "gray",top$gene, ""), "")  
top$labels = ifelse(top$gene %in% head(top[order(top$rank, decreasing = T),], n = 30)$gene, ifelse(top$colour != "gray",top$gene, ""), "")  




top20 <- top
top20$valid = ""
top20$valid [(top20$p_val_adj < 0.05)] <- TRUE
top20$valid [(top20$p_val_adj < 0.05)] <- TRUE
top20$valid[(top20$pct.1 < 0.2)] <- FALSE
top20$valid[(top20$pct.2 > 0.2)] <- FALSE
top20 <- subset(top20, subset = top20$valid == TRUE )

top20 <- top20 %>% group_by(cluster) %>% top_n(3000, avg_log2FC)
top20 <- top20[order(top20$cluster),]
top20 <- na.omit(top20)
#top20 <- subset(top20, colour %!in% c("gray"))

#top20$avg_log2FC[top20$cluster == "Donor"] <- top20$avg_log2FC *-1
top20$colour[top20$cluster == "Donor" & top20$colour != "grey"] <- "#17a4ad"
top20$labels = ifelse(top20$gene %in% head(top20[order(top20$rank, decreasing = T),], n = 20)$gene, ifelse(top20$colour != "gray",top20$gene, ""), "")  


red_cluster <- ScaleData(red_cluster, features=top20$gene)

plot_heatmap(red_cluster,
             markers = top20$gene,
             sort_var = c("hvd","seurat_clusters"),
             anno_var = c("hvd", "seurat_clusters", "sample.ident"),
             anno_colors= list(c("#17a4ad","#ad179d"), c("#b20505"), polychrome_pal),
             row_font_size = 6.5)



blue_cluster <- ScaleData(blue_cluster, features=top20$gene)

plot_heatmap(blue_cluster,
             markers = top20$gene,
             sort_var = c("hvd","seurat_clusters"),
             anno_var = c("hvd", "seurat_clusters", "sample.ident"),
             anno_colors= list(c("#17a4ad","#ad179d"), c("#0312ad"), polychrome_pal),
             row_font_size = 6.5)


green_cluster <- ScaleData(green_cluster, features=top20$gene)

plot_heatmap(green_cluster,
             markers = top20$gene,
             sort_var = c("hvd","seurat_clusters"),
             anno_var = c("hvd", "seurat_clusters", "sample.ident"),
             anno_colors= list(c("#17a4ad","#ad179d"), c("#03ad46"), polychrome_pal),
             row_font_size = 6.5)


top20 <- subset(top20, gene %!in% c("SFTPC"))
top20$avg_log2FC[top20$cluster == "Donor"] <- top20$avg_log2FC *-1
ggplot(top20, aes(x = top20$avg_log2FC, y = -log10(top20$p_val_adj)))+
  geom_text_repel(label = top20$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_point(color = top20$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")



top.red <- subset(top.red, gene %!in% c("SFTPC"))
top.red$avg_log2FC[top.red$cluster == "Donor"] <- top20$avg_log2FC *-1

top.blue <- subset(top.blue, gene %!in% c("SFTPC"))
top.blue$avg_log2FC[top.blue$cluster == "Donor"] <- top20$avg_log2FC *-1

top.green <- subset(top.green, gene %!in% c("SFTPC"))
top.green$avg_log2FC[top.green$cluster == "Donor"] <- top20$avg_log2FC *-1

top.red <- subset(top.red, colour %!in% "gray")
top.blue <- subset(top.blue, colour %!in% "gray")
top.green <- subset(top.green, colour %!in% "gray")

top20 <- subset(top20, gene %!in% c("SFTPC"))
top20.up <- subset(top20, colour %in% "#ad179d")
top20.down <- subset(top20, colour %in% "#17a4ad")

top.red.up <- subset(top.red, colour %in% "#ad179d")
top.red.down <- subset(top.red, colour %in% "#17a4ad")

top.blue.up <- subset(top.blue, colour %in% "#ad179d")
top.blue.down <- subset(top.blue, colour %in% "#17a4ad")

top.green.up <- subset(top.green, colour %in% "#ad179d")
top.green.down <- subset(top.green, colour %in% "#17a4ad")

ggplot(top.red.up, aes(x = top.red.up$avg_log2FC, y = -log10(top.red.up$p_val_adj)))+
  geom_text_repel(label = top.red.up$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_point(color = top.red.up$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")
