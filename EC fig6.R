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

top <- top.red #only run 1 of these three lines depending on current cluster investigated 
top <- top.blue #only run 1 of these three lines depending on current cluster investigated 
top <- top.green #only run 1 of these three lines depending on current cluster investigated 

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

top20 <- top20 %>% group_by(cluster) %>% top_n(30, avg_log2FC)
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







#Enrichr for GSEA

library(devtools)
library(enrichR)
library(reshape2)
library(rrvgo)

top20 <- top
top20$valid = ""
top20$valid [(top20$p_val_adj < 0.05)] <- TRUE
top20$valid[(top20$pct.1 < 0.2)] <- FALSE
top20$valid[(top20$pct.2 > 0.2)] <- FALSE
top20 <- subset(top20, subset = top20$valid == TRUE )

top20 <- top20 %>% group_by(cluster) %>% top_n(3000, avg_log2FC)
top20 <- top20[order(top20$cluster),]
top20 <- na.omit(top20)
top20 <- subset(top20, colour %!in% c("gray"))


top20$colour[top20$cluster == "Donor" & top20$colour != "grey"] <- "#17a4ad"
top20$labels = ifelse(top20$gene %in% head(top20[order(top20$rank, decreasing = T),], n = 20)$gene, ifelse(top20$colour != "gray",top20$gene, ""), "")  


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
  enriched <- enrichr(top20$gene, dbs)#switch with top.sorted.c1.FP$gene, top.sorted.c2.FP$gene for each cluster
}

#check with classic graph
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value")
}
  
  
enriched$GO_Biological_Process_2023$Count <- str_extract(enriched$GO_Biological_Process_2023$Overlap, ".+?(?=/)")
enriched$GO_Biological_Process_2023$Count <- as.numeric(enriched$GO_Biological_Process_2023$Count)

enrichedplot <- head(enriched$GO_Biological_Process_2023, 100)
enrichedplot2 <- head(enriched$GO_Biological_Process_2023, 10)

# chosen.red  <- c("Regulation Of Cell Adhesion (GO:0030155)",
#                  "Positive Regulation Of Leukocyte Cell-Cell Adhesion (GO:1903039)",
#                  "Negative Regulation Of Metabolic Process (GO:0009892)",
#                  "Extracellular Matrix Organization (GO:0030198)",
#                  "Protein Stabilization (GO:0050821)",
#                  "Long-Chain Fatty Acid Transport (GO:0015909)",
#                  "Positive Regulation Of Leukocyte Cell-Cell Adhesion (GO:1903039)",
#                  "Positive Regulation Of T Cell Activation (GO:0050870)",
#                  "MHC Class II Protein Complex Assembly (GO:0002399)",
#                  "Negative Regulation Of Epithelial Cell Migration (GO:0010633)")


#Chrd diagram to show the genes of top curated ontologies
chordgenes <- subset(enrichedplot2, select=c("Term","Genes","Adjusted.P.value"))

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






#revigo
enrichr.terms <- as.data.frame(enriched$GO_Biological_Process_2023)
enrichr.terms <- as.data.frame(enrichedplot2)

enrichr.GO <- enrichr.terms
enrichr.GO <- extract(enrichr.GO, Term, into=c("Term","GO"), regex = "(.*?) \\((.*?)\\)")

simMatrix <- calculateSimMatrix(enrichr.GO$GO,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")



#enrichr.GO.red <- enrichr.GO
reduced <- rownames(simMatrix)

enrichr.GO.red <- enrichr.GO[enrichr.GO$GO %in% reduced, ]

scores <- -log(enrichr.GO.red$Adjusted.P.value)
names(scores) <- enrichr.GO.red$GO

#run revigo
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold = 0.1,
                                orgdb="org.Hs.eg.db")

#scatterPlot(simMatrix, reducedTerms)
#treemapPlot(reducedTerms)

#filter out parent ontologiesa
reducedTerms2 <- reducedTerms %>% arrange(desc(score))
reducedTerms2$dupes <- duplicated(reducedTerms2$parentTerm)

reducedTerms3 <- subset(reducedTerms2, dupes %in% c("FALSE"))

#revigoresult <- subset(reducedTerms2, select=c("parent", "parentSimScore", "parentTerm"))
revigoresult<- reducedTerms3 %>% arrange(desc(score))


#isolate genes from Gos
GOgenes <- enrichr.GO.red %>%
  mutate(Genes = strsplit(as.character(Genes), ";")) %>%
  unnest(Genes) %>%
  filter(Genes != "") 

top.enrichr <- top20 #change to current cluster
top.enrichr$Genes <- toupper(top.enrichr$gene)

GOgenes$l2fc <- ""


#GOgenes$l2fc <- ifelse(top$Genes == GOgenes$Genes, top$avg_log2FC, "FALSE")

# iic <- c(1:nrow(GOgenes))
# ic <- c(1:nrow(top.enrichr))
# 
# top.filter <- top20
# 
# for (ii in iic){
#   for (i in ic){
#     GOgenes$l2fc[ii] = ifelse(top.filter$gene[i] == GOgenes$Genes[ii], top.filter$avg_log2FC[i], GOgenes$l2fc[ii])
#     GOgenes$Adjusted.P.value[ii] = ifelse(top.filter$gene[i] == GOgenes$Genes[ii], top.filter$p_val_adj[i], GOgenes$Adjusted.P.value[ii]) 
#   }
# }

# Align GOgenes$Genes to top.filter$gene
idx <- match(GOgenes$Genes, top.filter$gene)

# Replace only matched rows
GOgenes$l2fc[!is.na(idx)] <- top.filter$avg_log2FC[idx[!is.na(idx)]]
GOgenes$Adjusted.P.value[!is.na(idx)] <- top.filter$p_val_adj[idx[!is.na(idx)]]


# GOgenes <- GOgenes %>%
#   left_join(
#     top.enrichr %>% 
#       select(Genes, avg_log2FC, p_val_adj),
#     by = "Genes"
#   ) %>%
#   mutate(
#     l2fc = if_else(!is.na(avg_log2FC), avg_log2FC, l2fc),
#     Adjusted.P.value = if_else(!is.na(p_val_adj), p_val_adj, Adjusted.P.value)
#   ) %>%
#   select(-avg_log2FC, -p_val_adj)

revigoresult <- revigoresult[order(revigoresult$score, decreasing=T),]
revigoresult$parentTerm <- tolower(revigoresult$parentTerm)#parent? 

#only take top 10 ontologies
TopGOI <- head(revigoresult, 10)
GOgenes$Term <- tolower(GOgenes$Term)
TopGOI$term <- tolower(TopGOI$term) 

#
TopGOI.genes <- subset(GOgenes, subset= Term %in% TopGOI$term)#parent?

TopGOI.genes$dupes <- duplicated(TopGOI.genes$Genes)

TopGOI.genes$l2fc <- as.numeric(TopGOI.genes$l2fc)


TopGOI.genes2 <- subset(TopGOI.genes, l2fc >= 1) #change to "<=-1" or ">=1" accordingly

TopGOI.genes2 <- TopGOI.genes

TopGOI.genes3 <- subset(TopGOI.genes2, dupes %in% c("FALSE"))

TopGOI.genes3 <- subset(TopGOI.genes2, select = c("Term", "Genes", "l2fc"))
TopGOI.genes5 <- subset(TopGOI.genes2, select = c("Term", "Genes", "l2fc", "Adjusted.P.value"))
TopGOI.genes5$l2fc <- as.numeric(TopGOI.genes5$l2fc)

test <- TopGOI.genes3
test$l2fc <- as.numeric(test$l2fc)

#"#17a4ad","#ad179d"

gridcolors <- c("#000000","#8e4806","#8e7906","#6f8e06", "#488e06", 
                "#8e0619", "#068e71", "#06638e", #"#50068e",#"#8e0686",
                ifelse(test %>% distinct(l2fc)>=0, "#ad179d", "#17a4ad"))

# gridcolors <- c("#000000","#8e4806","#8e7906","#6f8e06", "#488e06", 
#                 "#8e0619", "#068e71", "#06638e", "#50068e","#8e0686",
#                 ifelse(test %>% distinct(avg_log2FC)>=0, "#ad179d", "#17a4ad"))


circos.par(start.degree = 270)
chordDiagram(test, scale = FALSE, annotationTrack = "grid", link.lty = 1, grid.col = gridcolors)
circos.clear()









