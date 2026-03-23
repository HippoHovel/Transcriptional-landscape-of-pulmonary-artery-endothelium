














EC.story.0 <- subset(ECstory.subset.fil, ECstory.j_clusters %in% c("0"))
Idents(EC.story.0)
Idents(EC.story.0) <- EC.story.0$samplecondition
c0.PAHvPF.deg <- FindMarkers(EC.story.0, ident.1 ="IPAH", ident.2 = "IPF", test.use="MAST" )
c0.PAHvD.deg <- FindMarkers(EC.story.0, ident.1 ="IPAH", ident.2 = "Donor", test.use="MAST" )
c0.PFvD.deg <- FindMarkers(EC.story.0, ident.1 ="IPF", ident.2 = "Donor", test.use="MAST" )

c0.PAHvD.deg$disease <- "PAH"
c0.PFvD.deg$disease <- "PFPH"

#---

c0.PAHvD.deg2 <- c0.PAHvD.deg

c0.PAHvD.deg2 <- rownames_to_column(c0.PAHvD.deg2, var="gene")

c0.PAHvD.deg2$avg_log2FC <- as.numeric(c0.PAHvD.deg2$avg_log2FC)


c0.PAHvD.deg2$colour <- "#e06606" 



c0.PAHvD.deg2$rank <- abs(c0.PAHvD.deg2$avg_log2FC * -log10(c0.PAHvD.deg2$p_val_adj))

c0.PAHvD.deg2$labels = ifelse(c0.PAHvD.deg2$gene %in% head(c0.PAHvD.deg2[order(c0.PAHvD.deg2$rank, decreasing = T),], n = 100)$gene, ifelse(c0.PAHvD.deg2$colour != "gray",c0.PAHvD.deg2$gene, ""), "")  
c0.PAHvD.deg2$labels = ifelse(c0.PAHvD.deg2$gene %in% head(c0.PAHvD.deg2[order(c0.PAHvD.deg2$rank, decreasing = T),], n = 30)$gene, ifelse(c0.PAHvD.deg2$colour != "gray",c0.PAHvD.deg2$gene, ""), "")  

ggplot(c0.PAHvD.deg, aes(x = c0.PAHvD.deg2$avg_log2FC, y = -log10(c0.PAHvD.deg2$p_val_adj)))+
  geom_text_repel(label = c0.PAHvD.deg2$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_point(color = c0.PAHvD.deg2$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")

#---

c0.PFvD.deg2 <- c0.PFvD.deg

c0.PFvD.deg2 <- rownames_to_column(c0.PFvD.deg2, var="gene")

c0.PFvD.deg2$avg_log2FC <- as.numeric(c0.PFvD.deg2$avg_log2FC)


c0.PFvD.deg2$colour <- "#0491b3" 



c0.PFvD.deg2$rank <- abs(c0.PFvD.deg2$avg_log2FC * -log10(c0.PFvD.deg2$p_val_adj))

c0.PFvD.deg2$labels = ifelse(c0.PFvD.deg2$gene %in% head(c0.PFvD.deg2[order(c0.PFvD.deg2$rank, decreasing = T),], n = 100)$gene, ifelse(c0.PFvD.deg2$colour != "gray",c0.PFvD.deg2$gene, ""), "")  
c0.PFvD.deg2$labels = ifelse(c0.PFvD.deg2$gene %in% head(c0.PFvD.deg2[order(c0.PFvD.deg2$rank, decreasing = T),], n = 30)$gene, ifelse(c0.PFvD.deg2$colour != "gray",c0.PFvD.deg2$gene, ""), "")  

ggplot(c0.PFvD.deg, aes(x = c0.PFvD.deg2$avg_log2FC, y = -log10(c0.PFvD.deg2$p_val_adj)))+
  geom_text_repel(label = c0.PFvD.deg2$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_point(color = c0.PFvD.deg2$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")

#----------

pick <- c0.PFvD.deg2

pick$check <- ifelse(c0.PFvD.deg2$gene %in% c0.PAHvD.deg2$gene, T, F)
pick <- subset(pick, check == T)

reduce.c0.PFvD.deg2 <- subset(c0.PFvD.deg2, gene %in% pick$gene)
reduce.c0.PAHvD.deg2 <- subset(c0.PAHvD.deg2, gene %in% pick$gene)

reduce.c0.PFvD.deg2$colour <- "#e06606"
reduce.c0.PAHvD.deg2$colour <- "#0491b3"

# merge <- rbind(reduce.c0.PFvD.deg2, reduce.c0.PAHvD.deg2)
# merge <- subset(merge, p_val_adj != 1)

# merge1 <- subset(reduce.c0.PFvD.deg2, select = c("gene", "avg_log2FC", "disease"))
# merge2 <- subset(reduce.c0.PAHvD.deg2, select = c("gene", "avg_log2FC", "disease"))
# 
# merge <- rbind(merge1, merge2)

merge1 <- subset(reduce.c0.PFvD.deg2, select = c("gene", "avg_log2FC", "disease","p_val_adj", "colour"))
names(merge1)[names(merge1) == 'avg_log2FC'] <- 'PHPF_L2FC'
#merge1 <- subset(merge1, p_val_adj <= 0.05)

merge2 <- subset(reduce.c0.PAHvD.deg2, select = c("gene", "avg_log2FC", "disease","p_val_adj", "colour"))
names(merge2)[names(merge2) == 'avg_log2FC'] <- 'PAH_L2FC'
#merge2 <- subset(merge2, p_val_adj <= 0.05)

pick <- merge1

pick$check <- ifelse(merge1$gene %in% merge2$gene, T, F)
pick <- subset(pick, check == T)

merge1 <- subset(merge1, gene %in% pick$gene)
merge2 <- subset(merge2, gene %in% pick$gene)

# merge$PFL2FC <- merge1$avg_log2FC
# merge1$PAHL2FC <- merge2$PAHL2FC

names(merge1)[names(merge1) == 'p_val_adj'] <- 'PF_pval_adj'
names(merge2)[names(merge2) == 'p_val_adj'] <- 'PAH_pval_adj'

merge1 <- merge1[order(merge1$gene),]
merge2 <- merge2[order(merge2$gene),]

merge <- subset(merge1, select = c("gene", "PHPF_L2FC","PF_pval_adj"))  
merge$PAH_L2FC <- merge2$PAH_L2FC
merge$PAH_pval_adj  <- merge2$PAH_pval_adj 

merge$rankPF <- abs(merge$PHPF_L2FC * -log10(merge$PF_pval_adj))
merge$rankPAH <- abs(merge$PAH_L2FC * -log10(merge$PAH_pval_adj))

merge$labels = ifelse(merge$gene %in% head(merge[order(merge$rankPF, decreasing = T),], n = 15)$gene, merge$gene, "")  
merge$topPF = ifelse(merge$gene %in% head(merge[order(merge$rankPF, decreasing = T),], n = 15)$gene, T, "") 
merge$topPAH = ifelse(merge$gene %in% head(merge[order(merge$rankPAH, decreasing = T),], n = 15)$gene, T, "")

merge$labels <- ""
merge$labels <- ifelse((merge$topPAH == T | merge$topPF == T), merge$gene, "")

# pick <- c(c0.PFvD.deg2, c0.PAHvD.deg2)
# pick <- as.data.frame(pick)
# pick <- as.data.frame(c0.PFvD.deg2, c0.PAHvD.deg2)

# ggplot(pick, aes(x = reduce.c0.PFvD.deg2$avg_log2FC[reduce.c0.PFvD.deg2$gene %in%reduce.c0.PAHvD.deg2$gene], 
#                  y = reduce.c0.PAHvD.deg2$avg_log2FC[reduce.c0.PFvD.deg2$gene %in%reduce.c0.PAHvD.deg2$gene]))+
#   geom_point(color = merge$colour, alpha = 0.5, show.legend = TRUE)+
#   geom_hline(yintercept = 0.5, linetype = "dashed") +
#   geom_vline(xintercept = 0.5, linetype = "dashed") +
#   geom_hline(yintercept = -0.5, linetype = "dashed") +
#   geom_vline(xintercept = -0.5, linetype = "dashed") +
#   theme_bw()+
#   ylim(-7.5, 7.5)+
#   xlim(-7.5, 7.5)+
#   xlab("l2fc PAH") +
#   ylab("l2fc PHPF")


ggplot(merge, aes(x = merge$PAH_L2FC, 
                  y = merge$PHPF_L2FC))+
  geom_point(color = "red", alpha = 0.5, show.legend = TRUE)+
  geom_text_repel(label = merge$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_bw()+
  ylim(-7.5, 7.5)+
  xlim(-7.5, 7.5)+
  xlab("l2fc PAH") +
  ylab("l2fc PHPF")


mergedPH <- subset(merge,(merge$PF_pval_adj <= 0.05) | (merge$PAH_pval_adj <= 0.05))
mergedPH <- merge[(merge$PF_pval_adj <= 0.05) | (merge$PAH_pval_adj <= 0.05)]

ggplot(mergedPH, aes(x = mergedPH$PAH_L2FC, 
                   y = mergedPH$PHPF_L2FC))+
  geom_point(color = mergedPH$colour, alpha = 0.5, show.legend = TRUE)+
  #geom_text_repel(label = mergedPH$labels2, nudge_y = 3, na.rm = FALSE, size=2.5, max.overlaps = 100, color=mergedPH$colour) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_bw()+
  ylim(-7.5, 7.5)+
  xlim(-7.5, 7.5)+
  xlab("l2fc PAH") +
  ylab("l2fc PHPF")

mergedPH$colour <- ""

#mergedPH$colour[(mergedPH$PF_pval_adj <= 0.05) | (mergedPH$PAH_pval_adj <= 0.05)] <- "black"

mergedPH$colour[(mergedPH$PF_pval_adj < mergedPH$PAH_pval_adj)] <- "#0491b3"
mergedPH$colour[(mergedPH$PAH_pval_adj < mergedPH$PF_pval_adj)] <- "#e06606"
mergedPH$colour[(mergedPH$PF_pval_adj <= 0.05) & (mergedPH$PAH_pval_adj <= 0.05)] <- "#109e09"
#mergedPH$colour[(mergedPH$PF_pval_adj >= 0.05) & (mergedPH$PAH_pval_adj >= 0.05)] <- "grey"


mergedPH$topPF = ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPF, decreasing = T),], n = 15)$gene, T, "") 
mergedPH$topPAH = ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPAH, decreasing = T),], n = 15)$gene, T, "")

mergedPH$labels <- ""
mergedPH$labels <- ifelse((mergedPH$topPAH == T | mergedPH$topPF == T), mergedPH$gene, "")

mergedPH$labels2 <- ""
mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPAH, decreasing = T),], n = 50)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), "")
mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPF, decreasing = T),], n = 50)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), "")

mergedPH$labels2 <- ""
mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPAH, decreasing = T),], n = 60)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPF, decreasing = T),], n = 60)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), ""))

mergedPH2 <- subset(mergedPH, colour %!in% "#109e09")                 
mergedPH2$labels2 <- ifelse(mergedPH2$gene %in% head(mergedPH2[order(mergedPH2$rankPAH, decreasing = T),], n = 30)$gene, ifelse(mergedPH2$colour != "#109e09",mergedPH2$gene, ""), mergedPH2$labels2 <- ifelse(mergedPH2$gene %in% head(mergedPH2[order(mergedPH2$rankPF, decreasing = T),], n = 30)$gene, ifelse(mergedPH2$colour != "#109e09",mergedPH2$gene, ""), ""))

mergedPH3 <- subset(mergedPH2, labels2 %!in% "") 

ggplot(mergedPH3, aes(x = mergedPH3$PAH_L2FC, 
                    y = mergedPH3$PHPF_L2FC))+
  geom_point(color = mergedPH3$colour, alpha = 0.5, show.legend = TRUE)+
  geom_text_repel(label = mergedPH3$labels2, nudge_y = 0 , na.rm = FALSE, size=2, max.overlaps = 100, color=mergedPH3$colour) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_bw()+
  ylim(-7.5, 7.5)+
  xlim(-7.5, 7.5)+
  xlab("l2fc PAH") +
  ylab("l2fc PHPF")

#---------------

#mergedPH2$PHPF_L2FC[(mergedPH2$colour == "#0491b3") & (mergedPH2$PHPF_L2FC <= 0.5)]

mergedPHfil <- subset(mergedPH2, colour == "#0491b3") 

mergedPHfil2 <- subset(mergedPH2, colour == "#e06606") 


mergedPH.curated <- c("IGHGP", "IGHM","LERFS","RGS2","RGS1", "APOC1", "SPP1","APOE", "FOS", "FOSB", "HSPA1A", "HSPA1B", "JUN", "SRPX",
                    "RHOB", "TSC22D3", "C1QB", "LINC00910", "DEPP1", "DNAJB1","KLF6","CSRNP1", "CTSB","CTSD" ,
                    "DDIT4", "FTL","KLF6","MAP1LC3B","CXCL8", "CXCR4",#,"SPARCL1", "SPDYE2")
                    
                    "")


#--------------



mergedPH4 <- subset(mergedPH3, select = c("gene","PAH_L2FC","PHPF_L2FC"))

mergedPH5 <- mergedPH4 %>% pivot_longer(!gene, names_to="disease", values_to = "l2fc")

ggplot(mergedPH5, aes(disease, l2fc, fill=disease)) +
  geom_point()+ 
  geom_line(aes(group=gene)) 


mergedPH5$colour <- ""
PAHGOI <- mergedPH3$gene[mergedPH3$colour == "#e06606"]
mergedPH5$colour[(mergedPH5$gene %in% PAHGOI) & (mergedPH5$disease == "PAH_L2FC")] <- "#e06606"
mergedPH5$colour[(mergedPH5$gene %in% PAHGOI) & (mergedPH5$disease != "PAH_L2FC")] <- "grey"

PHPFGOI <- mergedPH3$gene[mergedPH3$colour == "#0491b3"]
mergedPH5$colour[(mergedPH5$gene %in% PHPFGOI) & (mergedPH5$disease == "PHPF_L2FC")] <- "#0491b3"
mergedPH5$colour[(mergedPH5$gene %in% PHPFGOI) & (mergedPH5$disease != "PHPF_L2FC")] <- "grey"  

mergedPH5$label <- ""
mergedPH5$label <- ifelse(mergedPH5$colour != "grey", mergedPH5$gene, "")


ggplot() +
  geom_point(aes(mergedPH5$disease, mergedPH5$l2fc, fill=mergedPH5$disease),
             color =mergedPH5$colour, size=2)+
  geom_point(aes(mergedPH5$disease, mergedPH5$l2fc, fill=mergedPH5$disease),
             color ="black", shape=1, size=2, alpha=0.6)+ 
  geom_line(aes(mergedPH5$disease, mergedPH5$l2fc, group=mergedPH5$gene), 
            color ="black", size=0.4, alpha=0.4, linetype="dashed")+
  #geom_text_repel(aes(mergedPH5$disease, mergedPH5$l2fc),
  #label = mergedPH5$label, nudge_x = 0, na.rm = T, size=2.5, max.overlaps = 100, color=mergedPH5$colour) +
  theme_bw()+
  xlab("diseases") +
  ylab("logfoldchange")+
  #ylim(-4,4.5)+
  NoLegend()








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



#cluster 0 pick
enrichedplot2 <- subset(enrichedplot2, Term %in% c("Positive Regulation Of Cell Migration (GO:0030335)",
                                                   "Cytoplasmic Translation (GO:0002181)",
                                                   "protein-RNA Complex Assembly (GO:0022618)",
                                                   "Gene Expression (GO:0010467)",
                                                   "Positive Regulation Of Cell-Substrate Adhesion (GO:0010811)"))


enrichedplot4 <- subset(enrichedplot4, Term %in% c("Negative Regulation Of ERBB Signaling Pathway (GO:1901185)",
                                                   "Positive Regulation Of MAPK Cascade (GO:0043410)",
                                                   "Negative Regulation Of Epidermal Growth Factor Receptor Signaling Pathway (GO:0042059)",
                                                   "Long-Chain Fatty Acid Transport (GO:0015909)",
                                                   "Eosinophil Chemotaxis (GO:0048245)"))








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


















#retry only cluster 1
EC.story.1 <- subset(ECstory.subset.fil, ECstory.j_clusters %in% c("1"))
Idents(EC.story.1)
Idents(EC.story.1) <- EC.story.1$samplecondition
#c1.PAHvPF.deg <- FindMarkers(EC.story.1, ident.1 ="IPAH", ident.2 = "IPF", test.use="MAST" )
c1.PAHvD.deg <- FindMarkers(EC.story.1, ident.1 ="IPAH", ident.2 = "Donor", test.use="MAST" )
c1.PFvD.deg <- FindMarkers(EC.story.1, ident.1 ="IPF", ident.2 = "Donor", test.use="MAST" )

c1.PAHvD.deg$disease <- "PAH"
c1.PFvD.deg$disease <- "PFPH"

#---

c1.PAHvD.deg2 <- c1.PAHvD.deg

c1.PAHvD.deg2 <- rownames_to_column(c1.PAHvD.deg2, var="gene")

c1.PAHvD.deg2$avg_log2FC <- as.numeric(c1.PAHvD.deg2$avg_log2FC)


c1.PAHvD.deg2$colour <- "#e06606" 



c1.PAHvD.deg2$rank <- abs(c1.PAHvD.deg2$avg_log2FC * -log10(c1.PAHvD.deg2$p_val_adj))

c1.PAHvD.deg2$labels = ifelse(c1.PAHvD.deg2$gene %in% head(c1.PAHvD.deg2[order(c1.PAHvD.deg2$rank, decreasing = T),], n = 100)$gene, ifelse(c1.PAHvD.deg2$colour != "gray",c1.PAHvD.deg2$gene, ""), "")  
c1.PAHvD.deg2$labels = ifelse(c1.PAHvD.deg2$gene %in% head(c1.PAHvD.deg2[order(c1.PAHvD.deg2$rank, decreasing = T),], n = 30)$gene, ifelse(c1.PAHvD.deg2$colour != "gray",c1.PAHvD.deg2$gene, ""), "")  

ggplot(c1.PAHvD.deg, aes(x = c1.PAHvD.deg2$avg_log2FC, y = -log10(c1.PAHvD.deg2$p_val_adj)))+
  geom_text_repel(label = c1.PAHvD.deg2$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_point(color = c1.PAHvD.deg2$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")

#---

c1.PFvD.deg2 <- c1.PFvD.deg

c1.PFvD.deg2 <- rownames_to_column(c1.PFvD.deg2, var="gene")

c1.PFvD.deg2$avg_log2FC <- as.numeric(c1.PFvD.deg2$avg_log2FC)


c1.PFvD.deg2$colour <- "#0491b3" 



c1.PFvD.deg2$rank <- abs(c1.PFvD.deg2$avg_log2FC * -log10(c1.PFvD.deg2$p_val_adj))

c1.PFvD.deg2$labels = ifelse(c1.PFvD.deg2$gene %in% head(c1.PFvD.deg2[order(c1.PFvD.deg2$rank, decreasing = T),], n = 100)$gene, ifelse(c1.PFvD.deg2$colour != "gray",c1.PFvD.deg2$gene, ""), "")  
c1.PFvD.deg2$labels = ifelse(c1.PFvD.deg2$gene %in% head(c1.PFvD.deg2[order(c1.PFvD.deg2$rank, decreasing = T),], n = 30)$gene, ifelse(c1.PFvD.deg2$colour != "gray",c1.PFvD.deg2$gene, ""), "")  

ggplot(c1.PFvD.deg, aes(x = c1.PFvD.deg2$avg_log2FC, y = -log10(c1.PFvD.deg2$p_val_adj)))+
  geom_text_repel(label = c1.PFvD.deg2$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_point(color = c1.PFvD.deg2$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")

#----------

pick <- c1.PFvD.deg2

pick$check <- ifelse(c1.PFvD.deg2$gene %in% c1.PAHvD.deg2$gene, T, F)
pick <- subset(pick, check == T)

reduce.c1.PFvD.deg2 <- subset(c1.PFvD.deg2, gene %in% pick$gene)
reduce.c1.PAHvD.deg2 <- subset(c1.PAHvD.deg2, gene %in% pick$gene)

reduce.c1.PFvD.deg2$colour <- "#e06606"
reduce.c1.PAHvD.deg2$colour <- "#0491b3"

# merge <- rbind(reduce.c1.PFvD.deg2, reduce.c1.PAHvD.deg2)
# merge <- subset(merge, p_val_adj != 1)

# merge1 <- subset(reduce.c1.PFvD.deg2, select = c("gene", "avg_log2FC", "disease"))
# merge2 <- subset(reduce.c1.PAHvD.deg2, select = c("gene", "avg_log2FC", "disease"))
# 
# merge <- rbind(merge1, merge2)

merge1 <- subset(reduce.c1.PFvD.deg2, select = c("gene", "avg_log2FC", "disease","p_val_adj", "colour"))
names(merge1)[names(merge1) == 'avg_log2FC'] <- 'PHPF_L2FC'
#merge1 <- subset(merge1, p_val_adj <= 0.05)

merge2 <- subset(reduce.c1.PAHvD.deg2, select = c("gene", "avg_log2FC", "disease","p_val_adj", "colour"))
names(merge2)[names(merge2) == 'avg_log2FC'] <- 'PAH_L2FC'
#merge2 <- subset(merge2, p_val_adj <= 0.05)

pick <- merge1

pick$check <- ifelse(merge1$gene %in% merge2$gene, T, F)
pick <- subset(pick, check == T)

merge1 <- subset(merge1, gene %in% pick$gene)
merge2 <- subset(merge2, gene %in% pick$gene)

# merge$PFL2FC <- merge1$avg_log2FC
# merge1$PAHL2FC <- merge2$PAHL2FC

names(merge1)[names(merge1) == 'p_val_adj'] <- 'PF_pval_adj'
names(merge2)[names(merge2) == 'p_val_adj'] <- 'PAH_pval_adj'

merge1 <- merge1[order(merge1$gene),]
merge2 <- merge2[order(merge2$gene),]

merge <- subset(merge1, select = c("gene", "PHPF_L2FC","PF_pval_adj"))  
merge$PAH_L2FC <- merge2$PAH_L2FC
merge$PAH_pval_adj  <- merge2$PAH_pval_adj 

merge$rankPF <- abs(merge$PHPF_L2FC * -log10(merge$PF_pval_adj))
merge$rankPAH <- abs(merge$PAH_L2FC * -log10(merge$PAH_pval_adj))

merge$labels = ifelse(merge$gene %in% head(merge[order(merge$rankPF, decreasing = T),], n = 15)$gene, merge$gene, "")  
merge$topPF = ifelse(merge$gene %in% head(merge[order(merge$rankPF, decreasing = T),], n = 15)$gene, T, "") 
merge$topPAH = ifelse(merge$gene %in% head(merge[order(merge$rankPAH, decreasing = T),], n = 15)$gene, T, "")

merge$labels <- ""
merge$labels <- ifelse((merge$topPAH == T | merge$topPF == T), merge$gene, "")

# pick <- c(c1.PFvD.deg2, c1.PAHvD.deg2)
# pick <- as.data.frame(pick)
# pick <- as.data.frame(c1.PFvD.deg2, c1.PAHvD.deg2)

# ggplot(pick, aes(x = reduce.c1.PFvD.deg2$avg_log2FC[reduce.c1.PFvD.deg2$gene %in%reduce.c1.PAHvD.deg2$gene], 
#                  y = reduce.c1.PAHvD.deg2$avg_log2FC[reduce.c1.PFvD.deg2$gene %in%reduce.c1.PAHvD.deg2$gene]))+
#   geom_point(color = merge$colour, alpha = 0.5, show.legend = TRUE)+
#   geom_hline(yintercept = 0.5, linetype = "dashed") +
#   geom_vline(xintercept = 0.5, linetype = "dashed") +
#   geom_hline(yintercept = -0.5, linetype = "dashed") +
#   geom_vline(xintercept = -0.5, linetype = "dashed") +
#   theme_bw()+
#   ylim(-7.5, 7.5)+
#   xlim(-7.5, 7.5)+
#   xlab("l2fc PAH") +
#   ylab("l2fc PHPF")


ggplot(merge, aes(x = merge$PAH_L2FC, 
                  y = merge$PHPF_L2FC))+
  geom_point(color = "red", alpha = 0.5, show.legend = TRUE)+
  geom_text_repel(label = merge$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_bw()+
  ylim(-7.5, 7.5)+
  xlim(-7.5, 7.5)+
  xlab("l2fc PAH") +
  ylab("l2fc PHPF")


mergedPH <- subset(merge,(merge$PF_pval_adj <= 0.05) | (merge$PAH_pval_adj <= 0.05))
#mergedPH <- merge[(merge$PF_pval_adj <= 0.05) | (merge$PAH_pval_adj <= 0.05)]



mergedPH$colour <- ""

#mergedPH$colour[(mergedPH$PF_pval_adj <= 0.05) | (mergedPH$PAH_pval_adj <= 0.05)] <- "black"

mergedPH$colour[(mergedPH$PF_pval_adj < mergedPH$PAH_pval_adj)] <- "#0491b3"
mergedPH$colour[(mergedPH$PAH_pval_adj < mergedPH$PF_pval_adj)] <- "#e06606"
mergedPH$colour[(mergedPH$PF_pval_adj <= 0.05) & (mergedPH$PAH_pval_adj <= 0.05)] <- "#109e09"

mergedPH$topPF = ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPF, decreasing = T),], n = 15)$gene, T, "") 
mergedPH$topPAH = ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPAH, decreasing = T),], n = 15)$gene, T, "")

mergedPH$labels2 <- ""
mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPAH, decreasing = T),], n = 50)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), "")
mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPF, decreasing = T),], n = 50)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), "")
# 
mergedPH$labels2 <- ""
mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPAH, decreasing = T),], n = 60)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPF, decreasing = T),], n = 60)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), ""))


ggplot(mergedPH, aes(x = mergedPH$PAH_L2FC, 
                     y = mergedPH$PHPF_L2FC))+
  geom_point(color = mergedPH$colour, alpha = 0.5, show.legend = TRUE)+
  geom_text_repel(label = mergedPH$labels2, nudge_y = 3, na.rm = FALSE, size=2.5, max.overlaps = 100, color=mergedPH$colour) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_bw()+
  ylim(-7.5, 7.5)+
  xlim(-7.5, 7.5)+
  xlab("l2fc PAH") +
  ylab("l2fc PHPF")



mergedPH2 <- subset(mergedPH, colour %!in% "#109e09")                 
mergedPH2$labels2 <- ifelse(mergedPH2$gene %in% head(mergedPH2[order(mergedPH2$rankPAH, decreasing = T),], n = 30)$gene, ifelse(mergedPH2$colour != "#109e09",mergedPH2$gene, ""), mergedPH2$labels2 <- ifelse(mergedPH2$gene %in% head(mergedPH2[order(mergedPH2$rankPF, decreasing = T),], n = 30)$gene, ifelse(mergedPH2$colour != "#109e09",mergedPH2$gene, ""), ""))

mergedPH3 <- subset(mergedPH2, labels2 %!in% "") 

ggplot(mergedPH3, aes(x = mergedPH3$PAH_L2FC, 
                      y = mergedPH3$PHPF_L2FC))+
  geom_point(color = mergedPH3$colour, alpha = 0.5, show.legend = TRUE)+
  #geom_text_repel(label = mergedPH3$labels2, nudge_y = 0 , na.rm = FALSE, size=2, max.overlaps = 100, color=mergedPH3$colour) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_bw()+
  ylim(-7.5, 7.5)+
  xlim(-7.5, 7.5)+
  xlab("l2fc PAH") +
  ylab("l2fc PHPF")

mergedPH4 <- subset(mergedPH3, select = c("gene","PAH_L2FC","PHPF_L2FC"))

mergedPH5 <- mergedPH4 %>% pivot_longer(!gene, names_to="disease", values_to = "l2fc")

ggplot(mergedPH5, aes(disease, l2fc, fill=disease)) +
  geom_point()+ 
  geom_line(aes(group=gene)) 


mergedPH5$colour <- ""
PAHGOI <- mergedPH3$gene[mergedPH3$colour == "#e06606"]
mergedPH5$colour[(mergedPH5$gene %in% PAHGOI) & (mergedPH5$disease == "PAH_L2FC")] <- "#e06606"
mergedPH5$colour[(mergedPH5$gene %in% PAHGOI) & (mergedPH5$disease != "PAH_L2FC")] <- "grey"

PHPFGOI <- mergedPH3$gene[mergedPH3$colour == "#0491b3"]
mergedPH5$colour[(mergedPH5$gene %in% PHPFGOI) & (mergedPH5$disease == "PHPF_L2FC")] <- "#0491b3"
mergedPH5$colour[(mergedPH5$gene %in% PHPFGOI) & (mergedPH5$disease != "PHPF_L2FC")] <- "grey"  

mergedPH5$label <- ""
mergedPH5$label <- ifelse(mergedPH5$colour != "grey", mergedPH5$gene, "")


ggplot() +
  geom_point(aes(mergedPH5$disease, mergedPH5$l2fc, fill=mergedPH5$disease),
             color =mergedPH5$colour, size=2)+
  geom_point(aes(mergedPH5$disease, mergedPH5$l2fc, fill=mergedPH5$disease),
             color ="black", shape=1, size=2, alpha=0.6)+ 
  geom_line(aes(mergedPH5$disease, mergedPH5$l2fc, group=mergedPH5$gene), 
            color ="black", size=0.4, alpha=0.4, linetype="dashed")+
  geom_text_repel(aes(mergedPH5$disease, mergedPH5$l2fc),
                  label = mergedPH5$label, nudge_x = 0, na.rm = T, size=2.5, max.overlaps = 100, color=mergedPH5$colour) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw()+
  xlab("diseases") +
  ylab("logfoldchange")+
  #ylim(-4,4.5)+
  NoLegend()






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




#cluster 1pick
enrichedplot2 <- subset(enrichedplot2, Term %in% c("Cytoplasmic Translation (GO:0002181)",
                                                   "Peptide Biosynthetic Process (GO:0043043)",
                                                   "Gene Expression (GO:0010467)",
                                                   "Positive Regulation Of Intracellular Signal Transduction (GO:1902533)",
                                                   "Regulation Of Apoptotic Process (GO:0042981)"))


enrichedplot4 <- subset(enrichedplot4, Term %in% c("Membrane Lipid Biosynthetic Process (GO:0046467)",
                                                   "Sphingolipid Biosynthetic Process (GO:0030148)",
                                                   "Regulation Of Nucleotide-Binding Oligomerization Domain Containing 2 Signaling Pathway (GO:0070432)",
                                                   "Regulation Of Endoplasmic Reticulum Stress-Induced Intrinsic Apoptotic Signaling Pathway (GO:1902235)",
                                                   "Organophosphate Ester Transport (GO:0015748)"))

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




#retry only cluster 2
EC.story.2 <- subset(ECstory.subset.fil, ECstory.j_clusters %in% c("2"))
Idents(EC.story.2)
Idents(EC.story.2) <- EC.story.2$samplecondition
#c2.PAHvPF.deg <- FindMarkers(EC.story.2, ident.1 ="IPAH", ident.2 = "IPF", test.use="MAST" )
c2.PAHvD.deg <- FindMarkers(EC.story.2, ident.1 ="IPAH", ident.2 = "Donor", test.use="MAST" )
c2.PFvD.deg <- FindMarkers(EC.story.2, ident.1 ="IPF", ident.2 = "Donor", test.use="MAST" )

c2.PAHvD.deg$disease <- "PAH"
c2.PFvD.deg$disease <- "PFPH"

#---

c2.PAHvD.deg2 <- c2.PAHvD.deg

c2.PAHvD.deg2 <- rownames_to_column(c2.PAHvD.deg2, var="gene")

c2.PAHvD.deg2$avg_log2FC <- as.numeric(c2.PAHvD.deg2$avg_log2FC)


c2.PAHvD.deg2$colour <- "#e06606" 



c2.PAHvD.deg2$rank <- abs(c2.PAHvD.deg2$avg_log2FC * -log10(c2.PAHvD.deg2$p_val_adj))

c2.PAHvD.deg2$labels = ifelse(c2.PAHvD.deg2$gene %in% head(c2.PAHvD.deg2[order(c2.PAHvD.deg2$rank, decreasing = T),], n = 100)$gene, ifelse(c2.PAHvD.deg2$colour != "gray",c2.PAHvD.deg2$gene, ""), "")  
c2.PAHvD.deg2$labels = ifelse(c2.PAHvD.deg2$gene %in% head(c2.PAHvD.deg2[order(c2.PAHvD.deg2$rank, decreasing = T),], n = 30)$gene, ifelse(c2.PAHvD.deg2$colour != "gray",c2.PAHvD.deg2$gene, ""), "")  

ggplot(c2.PAHvD.deg, aes(x = c2.PAHvD.deg2$avg_log2FC, y = -log10(c2.PAHvD.deg2$p_val_adj)))+
  geom_text_repel(label = c2.PAHvD.deg2$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_point(color = c2.PAHvD.deg2$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")

#---

c2.PFvD.deg2 <- c2.PFvD.deg

c2.PFvD.deg2 <- rownames_to_column(c2.PFvD.deg2, var="gene")

c2.PFvD.deg2$avg_log2FC <- as.numeric(c2.PFvD.deg2$avg_log2FC)


c2.PFvD.deg2$colour <- "#0491b3" 



c2.PFvD.deg2$rank <- abs(c2.PFvD.deg2$avg_log2FC * -log10(c2.PFvD.deg2$p_val_adj))

c2.PFvD.deg2$labels = ifelse(c2.PFvD.deg2$gene %in% head(c2.PFvD.deg2[order(c2.PFvD.deg2$rank, decreasing = T),], n = 100)$gene, ifelse(c2.PFvD.deg2$colour != "gray",c2.PFvD.deg2$gene, ""), "")  
c2.PFvD.deg2$labels = ifelse(c2.PFvD.deg2$gene %in% head(c2.PFvD.deg2[order(c2.PFvD.deg2$rank, decreasing = T),], n = 30)$gene, ifelse(c2.PFvD.deg2$colour != "gray",c2.PFvD.deg2$gene, ""), "")  

ggplot(c2.PFvD.deg, aes(x = c2.PFvD.deg2$avg_log2FC, y = -log10(c2.PFvD.deg2$p_val_adj)))+
  geom_text_repel(label = c2.PFvD.deg2$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_point(color = c2.PFvD.deg2$colour, alpha = 0.5, show.legend = TRUE)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw()+
  theme_linedraw() + 
  xlab("Fold change (log2)") +
  ylab("P-Value adjusted (-log10)")

#----------

pick <- c2.PFvD.deg2

pick$check <- ifelse(c2.PFvD.deg2$gene %in% c2.PAHvD.deg2$gene, T, F)
pick <- subset(pick, check == T)

reduce.c2.PFvD.deg2 <- subset(c2.PFvD.deg2, gene %in% pick$gene)
reduce.c2.PAHvD.deg2 <- subset(c2.PAHvD.deg2, gene %in% pick$gene)

reduce.c2.PFvD.deg2$colour <- "#e06606"
reduce.c2.PAHvD.deg2$colour <- "#0491b3"

# merge <- rbind(reduce.c2.PFvD.deg2, reduce.c2.PAHvD.deg2)
# merge <- subset(merge, p_val_adj != 1)

# merge1 <- subset(reduce.c2.PFvD.deg2, select = c("gene", "avg_log2FC", "disease"))
# merge2 <- subset(reduce.c2.PAHvD.deg2, select = c("gene", "avg_log2FC", "disease"))
# 
# merge <- rbind(merge1, merge2)

merge1 <- subset(reduce.c2.PFvD.deg2, select = c("gene", "avg_log2FC", "disease","p_val_adj", "colour"))
names(merge1)[names(merge1) == 'avg_log2FC'] <- 'PHPF_L2FC'
#merge1 <- subset(merge1, p_val_adj <= 0.05)

merge2 <- subset(reduce.c2.PAHvD.deg2, select = c("gene", "avg_log2FC", "disease","p_val_adj", "colour"))
names(merge2)[names(merge2) == 'avg_log2FC'] <- 'PAH_L2FC'
#merge2 <- subset(merge2, p_val_adj <= 0.05)

pick <- merge1

pick$check <- ifelse(merge1$gene %in% merge2$gene, T, F)
pick <- subset(pick, check == T)

merge1 <- subset(merge1, gene %in% pick$gene)
merge2 <- subset(merge2, gene %in% pick$gene)

# merge$PFL2FC <- merge1$avg_log2FC
# merge1$PAHL2FC <- merge2$PAHL2FC

names(merge1)[names(merge1) == 'p_val_adj'] <- 'PF_pval_adj'
names(merge2)[names(merge2) == 'p_val_adj'] <- 'PAH_pval_adj'

merge1 <- merge1[order(merge1$gene),]
merge2 <- merge2[order(merge2$gene),]

merge <- subset(merge1, select = c("gene", "PHPF_L2FC","PF_pval_adj"))  
merge$PAH_L2FC <- merge2$PAH_L2FC
merge$PAH_pval_adj  <- merge2$PAH_pval_adj 

merge$rankPF <- abs(merge$PHPF_L2FC * -log10(merge$PF_pval_adj))
merge$rankPAH <- abs(merge$PAH_L2FC * -log10(merge$PAH_pval_adj))

merge$labels = ifelse(merge$gene %in% head(merge[order(merge$rankPF, decreasing = T),], n = 15)$gene, merge$gene, "")  
merge$topPF = ifelse(merge$gene %in% head(merge[order(merge$rankPF, decreasing = T),], n = 15)$gene, T, "") 
merge$topPAH = ifelse(merge$gene %in% head(merge[order(merge$rankPAH, decreasing = T),], n = 15)$gene, T, "")

merge$labels <- ""
merge$labels <- ifelse((merge$topPAH == T | merge$topPF == T), merge$gene, "")



ggplot(merge, aes(x = merge$PAH_L2FC, 
                  y = merge$PHPF_L2FC))+
  geom_point(color = "red", alpha = 0.5, show.legend = TRUE)+
  geom_text_repel(label = merge$labels, nudge_y = 3, na.rm = FALSE, size=3, max.overlaps = 100) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_bw()+
  ylim(-7.5, 7.5)+
  xlim(-7.5, 7.5)+
  xlab("l2fc PAH") +
  ylab("l2fc PHPF")


mergedPH <- subset(merge,(merge$PF_pval_adj <= 0.05) | (merge$PAH_pval_adj <= 0.05))
#mergedPH <- merge[(merge$PF_pval_adj <= 0.05) | (merge$PAH_pval_adj <= 0.05)]

ggplot(mergedPH, aes(x = mergedPH$PAH_L2FC, 
                     y = mergedPH$PHPF_L2FC))+
  geom_point(color = mergedPH$colour, alpha = 0.5, show.legend = TRUE)+
  #geom_text_repel(label = mergedPH$labels2, nudge_y = 3, na.rm = FALSE, size=2.5, max.overlaps = 100, color=mergedPH$colour) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_bw()+
  ylim(-7.5, 7.5)+
  xlim(-7.5, 7.5)+
  xlab("l2fc PAH") +
  ylab("l2fc PHPF")

mergedPH$colour <- ""


enrichdat <- subset(mergedPH, colour %in% c("#0491b3"))


#mergedPH$colour[(mergedPH$PF_pval_adj <= 0.05) | (mergedPH$PAH_pval_adj <= 0.05)] <- "black"

mergedPH$colour[(mergedPH$PF_pval_adj < mergedPH$PAH_pval_adj)] <- "#0491b3"
mergedPH$colour[(mergedPH$PAH_pval_adj < mergedPH$PF_pval_adj)] <- "#e06606"
mergedPH$colour[(mergedPH$PF_pval_adj <= 0.05) & (mergedPH$PAH_pval_adj <= 0.05)] <- "#109e09"

mergedPH$topPF = ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPF, decreasing = T),], n = 15)$gene, T, "") 
mergedPH$topPAH = ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPAH, decreasing = T),], n = 15)$gene, T, "")

mergedPH$labels <- ""
mergedPH$labels <- ifelse((mergedPH$topPAH == T | mergedPH$topPF == T), mergedPH$gene, "")

mergedPH$labels2 <- ""
mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPAH, decreasing = T),], n = 50)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), "")
mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPF, decreasing = T),], n = 50)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), "")
# 
mergedPH$labels2 <- ""
mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPAH, decreasing = T),], n = 60)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), mergedPH$labels2 <- ifelse(mergedPH$gene %in% head(mergedPH[order(mergedPH$rankPF, decreasing = T),], n = 60)$gene, ifelse(mergedPH$colour != "#109e09",mergedPH$gene, ""), ""))

mergedPH2 <- subset(mergedPH, colour %!in% "#109e09")                 
mergedPH2$labels2 <- ifelse(mergedPH2$gene %in% head(mergedPH2[order(mergedPH2$rankPAH, decreasing = T),], n = 30)$gene, ifelse(mergedPH2$colour != "#109e09",mergedPH2$gene, ""), mergedPH2$labels2 <- ifelse(mergedPH2$gene %in% head(mergedPH2[order(mergedPH2$rankPF, decreasing = T),], n = 30)$gene, ifelse(mergedPH2$colour != "#109e09",mergedPH2$gene, ""), ""))

mergedPH3 <- subset(mergedPH2, labels2 %!in% "") 

ggplot(mergedPH3, aes(x = mergedPH3$PAH_L2FC, 
                      y = mergedPH3$PHPF_L2FC))+
  geom_point(color = mergedPH3$colour, alpha = 0.5, show.legend = TRUE)+
  #geom_text_repel(label = mergedPH3$labels2, nudge_y = 0 , na.rm = FALSE, size=2, max.overlaps = 100, color=mergedPH3$colour) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_bw()+
  ylim(-7.5, 7.5)+
  xlim(-7.5, 7.5)+
  xlab("l2fc PAH") +
  ylab("l2fc PHPF")

mergedPH4 <- subset(mergedPH3, select = c("gene","PAH_L2FC","PHPF_L2FC"))

mergedPH5 <- mergedPH4 %>% pivot_longer(!gene, names_to="disease", values_to = "l2fc")

ggplot(mergedPH5, aes(disease, l2fc, fill=disease)) +
  geom_point()+ 
  geom_line(aes(group=gene)) 


mergedPH5$colour <- ""
PAHGOI <- mergedPH3$gene[mergedPH3$colour == "#e06606"]
mergedPH5$colour[(mergedPH5$gene %in% PAHGOI) & (mergedPH5$disease == "PAH_L2FC")] <- "#e06606"
mergedPH5$colour[(mergedPH5$gene %in% PAHGOI) & (mergedPH5$disease != "PAH_L2FC")] <- "grey"

PHPFGOI <- mergedPH3$gene[mergedPH3$colour == "#0491b3"]
mergedPH5$colour[(mergedPH5$gene %in% PHPFGOI) & (mergedPH5$disease == "PHPF_L2FC")] <- "#0491b3"
mergedPH5$colour[(mergedPH5$gene %in% PHPFGOI) & (mergedPH5$disease != "PHPF_L2FC")] <- "grey"  

mergedPH5$label <- ""
mergedPH5$label <- ifelse(mergedPH5$colour != "grey", mergedPH5$gene, "")

mergedPH5 <- subset(mergedPH5, gene %!in% c("XIST"))#artifact from singular sample

ggplot() +
  geom_point(aes(mergedPH5$disease, mergedPH5$l2fc, fill=mergedPH5$disease),
             color =mergedPH5$colour, size=2)+
  geom_point(aes(mergedPH5$disease, mergedPH5$l2fc, fill=mergedPH5$disease),
             color ="black", shape=1, size=2, alpha=0.6)+ 
  geom_line(aes(mergedPH5$disease, mergedPH5$l2fc, group=mergedPH5$gene), 
            color ="black", size=0.4, alpha=0.4, linetype="dashed")+
  geom_text_repel(aes(mergedPH5$disease, mergedPH5$l2fc),
                  label = mergedPH5$label, nudge_x = 0, na.rm = T, size=2.5, max.overlaps = 100, color=mergedPH5$colour) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw()+
  xlab("diseases") +
  ylab("logfoldchange")+
  #ylim(-4,4.5)+
  NoLegend()







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


FeaturePlot_scCustom(ECstory.subset.fil, features="TGFB1", split.by = "samplecondition")
FeaturePlot_scCustom(ECstory.subset.fil, features="NRG1", split.by = "samplecondition")
FeaturePlot_scCustom(ECstory.subset.fil, features="MDK", split.by = "samplecondition")
FeaturePlot_scCustom(ECstory.subset.fil, features="LAMB2", split.by = "samplecondition")
FeaturePlot_scCustom(ECstory.subset.fil, features="HLA-DQA1", split.by = "samplecondition")
FeaturePlot_scCustom(ECstory.subset.fil, features="GRN", split.by = "samplecondition")
FeaturePlot_scCustom(ECstory.subset.fil, features="GAS6", split.by = "samplecondition")
FeaturePlot_scCustom(ECstory.subset.fil, features="APOE", split.by = "samplecondition")

Ligands_PAH <- c("TGFB1","NRG1","MDK","LAMB2","HLA-DQA1", "GRN", "GAS6", "APOE")
DotPlot_scCustom(ECstory.subset.fil, features=Ligands_PAH, group.by = "samplecondition", scale=FALSE)+coord_flip()

DotPlot_scCustom(
  ECstory.subset.fil,
  features = Ligands_PAH,
  group.by = "samplecondition",
  scale = FALSE
) +
  coord_flip() +
  scale_color_gradient(
    low = "grey80",
    high = "red"
  )



Ligands_PAH <- c("TGFB1","NRG1","MDK","LAMB2","HLA-DQA1", "GRN", "GAS6", "APOE")

DotPlot_scCustom(
  ECstory.subset.fil,
  features = Ligands_PAH,
  group.by = "samplecondition",
  scale = FALSE
) +
  coord_flip() +
  scale_color_gradient(
    low = "grey80",
    high = "orange"
  )

VlnPlot_scCustom(
  ECstory.subset.fil,
  features = Ligands_PAH,
  group.by = "samplecondition",
  pt.size=0.1
) +
  scale_color_gradient(
    low = "grey80",
    high = "orange"
  )


Ligands_Both <- c("SPP1","PTPRM","PPIA","POSTN","PECAM1", "MIF", "JAM2", "IGFBP3",
                  "ICAM2", "HLA-E", "HLA-DRB5", "HLA-DRB1", "HLA-DRA", "HLA-DQB1",
                  "HLA-DPB1", "HLA-DPA1", "HLA-DMA", "GJA1", "FN1", "EFNB2", "COL4A2",
                  "CHOLESTEROL -LIPA", "CD99", "CD55", "CD34", "CCL14", "BMP6", "APP",
                  "ANXA1", "ITGAV+ITGB1")

DotPlot_scCustom(
  ECstory.subset.fil,
  features = Ligands_Both,
  group.by = "samplecondition",
  scale = FALSE
) +
  coord_flip() +
  scale_color_gradient(
    low = "grey80",
    high = "green"
  )

VlnPlot_scCustom(
  ECstory.subset.fil,
  features = Ligands_Both,
  group.by = "samplecondition",
  pt.size=0
) +
  scale_color_gradient(
    low = "grey80",
    high = "green"
  )





Ligands_PHPF <- c("SEMA6A", "PGE2-PTGES3", "NAMPT", "LAMC1", "LAMB1", "ICAM1", 
                  "HLA-F", "HLA-C", "HLA-B", "HLA-A", "CXCL2", "CXCL12", "COL4A1", "CCL2")


DotPlot_scCustom(
  ECstory.subset.fil,
  features = Ligands_PHPF,
  group.by = "samplecondition",
  scale = FALSE
) +
  coord_flip() +
  scale_color_gradient(
    low = "grey80",
    high = "blue"
  )

VlnPlot_scCustom(
  ECstory.subset.fil,
  features = Ligands_PHPF,
  group.by = "samplecondition",
  pt.size=0
) +
  scale_color_gradient(
    low = "grey80",
    high = "blue"
  )


ECstory.j$cluster.annot <- ""
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "0"] <- "Fibroblasts"

ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "24"] <- "SMC2"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "14"] <- "SMC2"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "4"] <- "SMC1"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "9"] <- "SMC1"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "16"] <- "SMC1"

ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "12"] <- "Endo1"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "6"] <- "Endo2"
#ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "27"] <- "Endo MT"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "27"] <- "Endo2"

ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "23"] <- "Lymphatic Endo"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "20"] <- "Granulocytes"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "18"] <- "Dendritic cells"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "19"] <- "Dendritic cells"

ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "13"] <- "Mono/Macs"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "5"] <- "Mono/Macs"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "2"] <- "Mono/Macs"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "7"] <- "Mono/Macs"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "17"] <- "Mono/Macs"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "11"] <- "Mono/Macs"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "25"] <- "Mono/Macs"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "8"] <- "NK"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "1"] <- "T-cell"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "3"] <- "T-cell"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "10"] <- "T-cell"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "21"] <- "B-cell"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "22"] <- "Erythrocyte debris"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "15"] <- "Epithelial"
ECstory.j$cluster.annot[ECstory.j$seurat_clusters == "26"] <- "Epithelial"
DimPlot(ECstory.j, group.by = "cluster.annot", cols=polychrome_pal)
DimPlot(ECstory.j, group.by = "QCcut", cols=polychrome_pal)


cont.calc.pop <- subset(ECstory.j, cluster.annot %in% c("Endo1", "Endo2", "Epithelial"))

DimPlot(cont.calc.pop)


cont.calc.pop$cluster.annot2 <- ""
cont.calc.pop$cluster.annot2[ECstory.subset.fil$clusternames == "Immunomodulatory"] <- "Immunomodulatory"
cont.calc.pop$cluster.annot2[ECstory.subset.fil$clusternames == "Vascular tone-modulatory"] <- "Vascular tone-modulatory"
cont.calc.pop$cluster.annot2[ECstory.subset.fil$clusternames == "Stress response-modulatory"] <- "Stress response-modulatory"
DimPlot(cont.calc.pop, group.by = "cluster.annot2")

ECstory.subset.fil$clusternames

seu.full$celltype[colnames(seu.sub)] <- seu.sub$celltype


cont.calc.pop$cluster.annot2 <- ""
cont.calc.pop$cluster.annot2[colnames(ECstory.j)] <- ECstory.j$cluster.annot
cont.calc.pop$cluster.annot2[colnames(ECstory.subset.fil)] <- ECstory.subset.fil$clusternames

cont.calc.pop.fil <- subset(cont.calc.pop, cluster.annot2 %!in% c("Endo1", "Endo2"))
DimPlot(cont.calc.pop.fil, group.by = "cluster.annot2")


cont.calc.pop.fil <- RunUMAP(cont.calc.pop.fil, dims = 1:30, reduction= "integrated.cca", seed.use = 1)
DimPlot(cont.calc.pop.fil, group.by = "cluster.annot2", reduction = "umap")





#####calculate possible contamination per sample

epi_frac <- ECstory.j@meta.data %>%
  group_by(sample.ident) %>%
  summarise(
    n_total = n(),
    n_epi   = sum(cluster.annot == "Epithelial"),
    epi_pct = n_epi / n_total
  )

endo_counts <- ECstory.j@meta.data %>%
  filter(cluster.annot %in% c("Endo1", "Endo2")) %>%
  count(sample.ident, name = "n_endo")

contam_est <- left_join(endo_counts, epi_frac, by = "sample.ident") %>%
  mutate(
    est_contam_endo = epi_pct * n_endo,
    pct_contam_endo = epi_pct
  )

endo <- subset(ECstory.j, subset = cluster.annot %in% c("Endo1", "Endo2"))

endo$contam_prob <- epi_frac$epi_pct[
  match(endo$sample.ident, epi_frac$sample.ident)
]

endo$contam_bin <- cut(
  endo$contam_prob,
  breaks = c(0, 0.05, 0.15, 1),
  labels = c("Low", "Moderate", "High")
)

contam_summary <- contam_est %>%
  select(
    sample.ident,
    n_endo,
    n_epi,
    n_total,
    epi_pct,
    est_contam_endo,
    pct_contam_endo
  ) %>%
  arrange(desc(epi_pct))

View(contam_summary)

write.csv(contam_summary, "PA_epithelial_contamination_per_sample.csv",
          row.names = FALSE)


ggplot(contam_summary,
       aes(x = reorder(sample.ident, pct_contam_endo),
           y = pct_contam_endo)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    x = "Sample",
    y = "Estimated fraction of potentially contaminated ECs",
    title = "Upper-bound estimate of endothelial contamination per sample",
    subtitle = "Dashed line indicates 10% threshold (>90% ECs unaffected)"
  ) +
  theme_bw()


summary(contam_summary$epi_pct)
hist(contam_summary$epi_pct, breaks = 20,
     main = "Distribution of epithelial contamination",
     xlab = "Fraction epithelial cells")


epi_frac %>%
  mutate(contam_bin = cut(
    epi_pct,
    breaks = c(0, 0.05, 0.15, 1),
    labels = c("Low", "Moderate", "High")
  )) %>%
  count(contam_bin)


ggplot(endo@meta.data,
       aes(x = contam_bin, y = nCount_RNA)) +
  geom_boxplot() +
  labs(title = "Endothelial UMI counts vs epithelial contamination")


kruskal.test(nCount_RNA ~ contam_bin, data = endo@meta.data)

endo@meta.data %>%
  group_by(sample.ident) %>%
  summarise(
    n_endo = n(),
    contam = first(contam_prob)
  ) %>%
  cor.test(~ contam + n_endo, method = "spearman")




#plots

library(ggplot2)

ggplot(contam_summary,
       aes(x = reorder(sample.ident, epi_pct),
           y = epi_pct)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Sample",
    y = "Fraction epithelial cells",
    title = "Sample-level epithelial contamination"
  ) +
  theme_bw()

ggplot(contam_summary, aes(x = epi_pct)) +
  geom_histogram(bins = 20, fill = "grey70", color = "black") +
  labs(
    x = "Fraction epithelial cells",
    y = "Number of samples",
    title = "Distribution of epithelial contamination across samples"
  ) +
  theme_bw()


ggplot(endo@meta.data,
       aes(x = contam_bin, y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = "Epithelial contamination bin",
    y = "Endothelial UMI count",
    title = "Endothelial library complexity vs epithelial contamination"
  ) +
  theme_bw()


ggplot(endo@meta.data,
       aes(x = contam_bin, y = percent.mt)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = "Epithelial contamination bin",
    y = "Percent mitochondrial reads",
    title = "Endothelial mitochondrial content vs contamination"
  ) +
  theme_bw()

FeaturePlot(
  endo,
  features = "contam_prob",
  reduction = "umap.cca"
) +
  ggtitle("Endothelial cells colored by sample-level epithelial contamination")


p<- FeaturePlot(
  endo,
  features = "contam_prob",
  reduction = "umap.cca",
  split.by = "sample.ident",
  pt.size= 0.1,
  combine = F
)

wrap_plots(p, ncol = 4)

endo_yield <- endo@meta.data %>%
  group_by(sample.ident) %>%
  summarise(n_endo = n()) %>%
  left_join(epi_frac, by = "sample.ident")

ggplot(endo_yield,
       aes(x = epi_pct, y = n_endo)) +
  geom_point(size = 3) +
  labs(
    x = "Fraction epithelial cells",
    y = "Number of endothelial cells",
    title = "Endothelial yield vs epithelial contamination"
  ) +
  theme_bw()


write.csv(
  contam_summary,
  "epithelial_contamination_summary_per_sample.csv",
  row.names = FALSE
)


###use subset nomenclature

ECstory.j$cluster.annot2 <- ""
ECstory.j$cluster.annot2[colnames(ECstory.j)] <- ECstory.j$cluster.annot
ECstory.j$cluster.annot2[colnames(ECstory.subset.fil)] <- ECstory.subset.fil$clusternames
DimPlot(ECstory.j, group.by = "cluster.annot2")


epi_frac <- ECstory.j@meta.data %>%
  group_by(sample.ident) %>%
  summarise(
    n_total = n(),
    n_epi   = sum(cluster.annot2 == "Epithelial"),
    epi_pct = n_epi / n_total
  )

endo_counts <- ECstory.j@meta.data %>%
  filter(cluster.annot2 %in% c("Immunomodulatory", "Stress response-modulatory", "Vascular tone-modulatory")) %>%
  count(sample.ident, name = "n_endo")

contam_est <- left_join(endo_counts, epi_frac, by = "sample.ident") %>%
  mutate(
    est_contam_endo = epi_pct * n_endo,
    pct_contam_endo = epi_pct
  )

endo <- subset(ECstory.j, subset = cluster.annot2 %in% c("Immunomodulatory", "Stress response-modulatory", "Vascular tone-modulatory"))

endo$contam_prob <- epi_frac$epi_pct[
  match(endo$sample.ident, epi_frac$sample.ident)
]

endo$contam_bin <- cut(
  endo$contam_prob,
  breaks = c(0, 0.05, 0.15, 1),
  labels = c("Low", "Moderate", "High")
)

contam_summary <- contam_est %>%
  select(
    sample.ident,
    n_endo,
    n_epi,
    n_total,
    epi_pct,
    est_contam_endo,
    pct_contam_endo
  ) %>%
  arrange(desc(epi_pct))

# View(contam_summary)
# 
# write.csv(contam_summary, "PA_epithelial_contamination_per_sample.csv",
#           row.names = FALSE)
# 
# 
# ggplot(contam_summary,
#        aes(x = reorder(sample.ident, pct_contam_endo),
#            y = pct_contam_endo)) +
#   geom_col(fill = "steelblue") +
#   geom_hline(yintercept = 0.10, linetype = "dashed", color = "red") +
#   coord_flip() +
#   labs(
#     x = "Sample",
#     y = "Estimated fraction of potentially contaminated ECs",
#     title = "Upper-bound estimate of endothelial contamination per sample",
#     subtitle = "Dashed line indicates 10% threshold (>90% ECs unaffected)"
#   ) +
#   theme_bw()
# 
# 
# summary(contam_summary$epi_pct)
# hist(contam_summary$epi_pct, breaks = 20,
#      main = "Distribution of epithelial contamination",
#      xlab = "Fraction epithelial cells")
# 
# 
# epi_frac %>%
#   mutate(contam_bin = cut(
#     epi_pct,
#     breaks = c(0, 0.05, 0.15, 1),
#     labels = c("Low", "Moderate", "High")
#   )) %>%
#   count(contam_bin)
# 
# 
# ggplot(endo@meta.data,
#        aes(x = contam_bin, y = nCount_RNA)) +
#   geom_boxplot() +
#   labs(title = "Endothelial UMI counts vs epithelial contamination")
# 
# 
# kruskal.test(nCount_RNA ~ contam_bin, data = endo@meta.data)
# 
# endo@meta.data %>%
#   group_by(sample.ident) %>%
#   summarise(
#     n_endo = n(),
#     contam = first(contam_prob)
#   ) %>%
#   cor.test(~ contam + n_endo, method = "spearman")




#plots

library(ggplot2)

ggplot(contam_summary,
       aes(x = reorder(sample.ident, epi_pct),
           y = epi_pct)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Sample",
    y = "Fraction epithelial cells",
    title = "Sample-level epithelial contamination"
  ) +
  theme_bw()

ggplot(contam_summary, aes(x = epi_pct)) +
  geom_histogram(bins = 20, fill = "grey70", color = "black") +
  labs(
    x = "Fraction epithelial cells",
    y = "Number of samples",
    title = "Distribution of epithelial contamination across samples"
  ) +
  theme_bw()


ggplot(endo@meta.data,
       aes(x = contam_bin, y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = "Epithelial contamination bin",
    y = "Endothelial UMI count",
    title = "Endothelial library complexity vs epithelial contamination"
  ) +
  theme_bw()


ggplot(endo@meta.data,
       aes(x = contam_bin, y = percent.mt)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = "Epithelial contamination bin",
    y = "Percent mitochondrial reads",
    title = "Endothelial mitochondrial content vs contamination"
  ) +
  theme_bw()

FeaturePlot(
  endo,
  features = "contam_prob",
  reduction = "umap.cca"
) +
  ggtitle("Endothelial cells colored by sample-level epithelial contamination")


p<- FeaturePlot(
  endo,
  features = "contam_prob",
  reduction = "umap.cca",
  split.by = "sample.ident",
  pt.size= 0.1,
  combine = F
)

wrap_plots(p, ncol = 4)

endo_yield <- endo@meta.data %>%
  group_by(sample.ident) %>%
  summarise(n_endo = n()) %>%
  left_join(epi_frac, by = "sample.ident")

ggplot(endo_yield,
       aes(x = epi_pct, y = n_endo)) +
  geom_point(size = 3) +
  labs(
    x = "Fraction epithelial cells",
    y = "Number of endothelial cells",
    title = "Endothelial yield vs epithelial contamination"
  ) +
  theme_bw()


write.csv(
  contam_summary,
  "epithelial_contamination_summary_per_sample.csv",
  row.names = FALSE
)





ec_comp <- endo@meta.data %>%
  group_by(sample.ident, cluster.annot2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample.ident) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

ec_comp <- ec_comp %>%
  left_join(epi_frac, by = "sample.ident")

library(dplyr)

cluster_contam_stats <- ec_comp %>%
  group_by(cluster.annot2) %>%
  summarise(
    rho = cor(epi_pct, frac, method = "spearman"),
    pval = cor.test(epi_pct, frac, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  arrange(desc(abs(rho)))


cluster_lm_stats <- ec_comp %>%
  group_by(cluster.annot2) %>%
  summarise(
    r2 = summary(lm(frac ~ epi_pct))$r.squared,
    slope = coef(lm(frac ~ epi_pct))[2],
    pval = summary(lm(frac ~ epi_pct))$coefficients[2, 4],
    .groups = "drop"
  ) %>%
  arrange(desc(r2))

cluster_contam_stats


set.seed(1)

perm_results <- ec_comp %>%
  group_by(cluster.annot2) %>%
  summarise(
    perm_p = {
      df <- cur_data()
      obs <- cor(df$epi_pct, df$frac, method = "spearman")
      perm <- replicate(
        1000,
        cor(sample(df$epi_pct), df$frac, method = "spearman")
      )
      mean(abs(perm) >= abs(obs))
    },
    .groups = "drop"
  )

cluster_stats <- cluster_contam_stats %>%
  left_join(cluster_lm_stats, by = "cluster.annot2") %>%
  left_join(perm_results, by = "cluster.annot2") %>%
  arrange(desc(abs(rho)))



lm(frac ~ epi_pct)

set.seed(1)

perm_test <- function(df, n = 1000) {
  obs <- cor(df$epi_pct, df$frac, method = "spearman")
  perm <- replicate(n, cor(sample(df$epi_pct), df$frac, method = "spearman"))
  mean(abs(perm) >= abs(obs))
}


ggplot(ec_comp, aes(x = epi_pct, y = frac)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ cluster.annot2, scales = "free_y") +
  labs(
    x = "Fraction epithelial cells per sample",
    y = "Fraction of endothelial cells",
    title = "Endothelial subcluster abundance vs epithelial contamination"
  ) +
  theme_classic()

ggplot(cluster_stats, aes(x = cluster.annot2, y = rho)) +
  geom_col(fill = "grey60") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    x = "Endothelial subcluster",
    y = "Spearman correlation with epithelial contamination",
    title = "Coupling between EC subcluster abundance and contamination"
  ) +
  theme_classic()




ggplot(ec_comp, aes(x = epi_pct, y = frac)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
  facet_wrap(~ cluster.annot2, scales = "free_y") +
  labs(
    x = "Fraction epithelial cells per sample",
    y = "Fraction of endothelial cells"
  ) +
  theme_classic()






annot_df <- cluster_stats %>%
  mutate(
    label = paste0(
      "ρ = ", round(rho, 2),
      "\nperm p = ", round(perm_p, 2)
    )
  )

ggplot(ec_comp, aes(x = epi_pct, y = frac)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ cluster.annot2, scales = "free_y") +
  geom_text(
    data = annot_df,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.1,
    inherit.aes = FALSE
  ) +
  theme_classic()



ggplot(ec_comp, aes(x = epi_pct, y = frac)) +
  geom_point(size = 2) +
  facet_wrap(~ cluster.annot2, scales = "free_y") +
  theme_classic()


ggplot(contam_summary,
              aes(x = reorder(sample.ident, pct_contam_endo),
                  y = pct_contam_endo)) +
         geom_col(fill = "steelblue") +
         geom_hline(yintercept = 0.10, linetype = "dashed", color = "red") +
         coord_flip() +
         labs(
           x = "Sample",
           y = "Estimated fraction of potentially contaminated ECs",
           title = "Upper-bound estimate of endothelial contamination per sample",
           subtitle = "Dashed line indicates 10% threshold (>90% ECs unaffected)"
         ) +
         theme_bw()



# EC-only metadata
ec_counts <- endo@meta.data %>%
  group_by(sample.ident, cluster.annot2) %>%
  summarise(n_cells = n(), .groups = "drop")


sample_order <- epi_frac %>%
  arrange(desc(epi_pct)) %>%
  pull(sample.ident)

ec_counts$sample.ident <- factor(
  ec_counts$sample.ident,
  levels = rev(sample_order)
)

library(ggplot2)

ggplot(ec_counts,
       aes(x = sample.ident, y = n_cells, fill = cluster.annot2)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Sample (ordered by epithelial contamination)",
    y = "Number of endothelial cells",
    fill = "EC subcluster",
    title = "Endothelial subcluster cell counts per sample"
  ) +
  theme_classic()


ec_frac_plot <- ec_counts %>%
  group_by(sample.ident) %>%
  mutate(frac = n_cells / sum(n_cells))

ggplot(ec_frac_plot,
       aes(x = sample.ident, y = frac, fill = cluster.annot2)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Sample (ordered by epithelial contamination)",
    y = "Fraction of endothelial cells",
    fill = "EC subcluster",
    title = "Endothelial subcluster composition per sample"
  ) +
  theme_classic()



#no integration for batch correction proof

ECstory.j2 <- ECstory.j

ECstory.j2 <- NormalizeData(ECstory.j2)
ECstory.j2 <- FindVariableFeatures(ECstory.j2)
ECstory.j2 <- ScaleData(ECstory.j2)
ECstory.j2 <- RunPCA(ECstory.j2)


ECstory.j2 <- FindNeighbors(ECstory.j2, reduction = "pca", dims = 1:30)
ECstory.j2 <- FindClusters(ECstory.j2, resolution = 0.5, cluster.name = "pca")

ECstory.j2 <- RunUMAP(ECstory.j2, reduction = "pca", dims = 1:30)
ECstory.j2$cluster.identity[ECstory.j2$cluster.identity == "Endo MT"] <- "Endothelial cells"
DimPlot(ECstory.j2, group.by= "lineage", cols=polychrome_pal, reduction = "umap")
DimPlot(ECstory.j2, group.by= "sample.ident", cols=polychrome_pal, reduction = "umap")

ECstory.j$cluster.identity[ECstory.j$cluster.identity == "Endo MT"] <- "Endothelial cells"
DimPlot(ECstory.j, group.by= "cluster.identity", cols=polychrome_pal)
DimPlot(ECstory.j, group.by= "sample.ident", cols=polychrome_pal, pt.size=0.001, order=F)

DimPlot(ECstory.j, group.by= "cluster.identity", cols=polychrome_pal, pt.size=0.001, order=F)
DimPlot(ECstory.j2, group.by= "cluster.identity", cols=polychrome_pal, pt.size=0.001, order=F, reduction="umap")

DimPlot(ECstory.j, group.by= "samplecondition", pt.size=0.001, order=F)
DimPlot(ECstory.j2, group.by= "samplecondition", pt.size=0.001, order=F, reduction="umap")
