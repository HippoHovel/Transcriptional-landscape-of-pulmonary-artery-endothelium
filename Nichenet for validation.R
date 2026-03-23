library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)

DimPlot(ECstory.fil, group.by = "cellchat.ident")
ECstory.cellchat.fil <- subset(ECstory.fil, cellchat.ident %!in% c("Platelets", "Endo MT", "Erythrocytes", "EC3", "EC4"))
ECstory.cellchat.fil <- JoinLayers(ECstory.cellchat.fil)
DimPlot(ECstory.cellchat.fil, group.by = "cellchat.ident")
#mylevels <- c("gCap", "early fibrotic gCap", "late fibrotic gCap", "aCap", "fibrotic aCap", "Artery", "Vein", "Lymphatic", "B-cell", "T-cell", "Monocytes","alv. Mac", "int. Mac", "Granulocytes", "Dendritic", "NK", "Epithelial", "Mesenchymal")

ECstory.cellchat.fil$nichenet <- ECstory.cellchat.fil$cellchat.ident 
DimPlot(ECstory.cellchat.fil , group.by = "nichenet")

seuratObj <- ECstory.cellchat.fil  
Idents(seuratObj) <- seuratObj$nichenet

DimPlot(seuratObj, group.by = "nichenet")
seuratObj <- JoinLayers(seuratObj)
seuratObj<- subset(seuratObj, nichenet %!in% c("Lymphatic") )
# seuratObj$nichenet2 <- seuratObj$nichenet
# seuratObj$nichenet2 <- ""
#seuratObj$nichenet2<- levels(droplevels(seuratObj$nichenet2))


install.packages("curl") # if not installed
library(curl)


base <- "https://zenodo.org/records/7074291/files/"


curl_download(
  paste0(base, "ligand_target_matrix_nsga2r_final.rds?download=1"),
  "ligand_target_matrix_nsga2r_final.rds"
)


curl_download(
  paste0(base, "weighted_networks_nsga2r_final.rds?download=1"),
  "weighted_networks_nsga2r_final.rds"
)


curl_download(
  paste0(base, "lr_network_human_21122021.rds?download=1"),
  "lr_network_human_21122021.rds"
)

files <- c(
  lt = "ligand_target_matrix_nsga2r_final.rds",
  wn = "weighted_networks_nsga2r_final.rds",
  lr = "lr_network_human_21122021.rds"
)

lr_network <- readRDS(files["lr"])
ligand_target_matrix <- readRDS(files["lt"])
weighted_networks <- readRDS(files["wn"])


lr_network = lr_network %>% distinct(from, to)
head(lr_network)

weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network

###1. Define a “sender/niche” cell population and a “receiver/target”
## receiver
receiver = c("Endothelial")
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10, assay_oi = "RNA")


background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Epithelial cells", "Mono/Macs", "T-cells", "Dendritic cells", "SMC", "Fibroblasts",
                     "B-cells", "NK cells", "Endothelial", "Granulocytes")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10, assay_oi ="RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

###2. Define a gene set of interest:
seurat_obj_receiver= subset(seuratObj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["samplecondition", drop=TRUE]])

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

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")+ theme(legend.text=element_text(angle=45,hjust=1.2, vjust=1.4),legend.title = element_text(vjust = 0.9))

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
DE_table_all = Idents(seuratObj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratObj, condition_colname = "samplecondition", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
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
rotated_dotplot = DotPlot(seuratObj %>% subset(nichenet %in% sender_celltypes), features = order_ligands_adapted) + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
figures_without_legend = cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
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
DE_table <- calculate_de(seuratObj, celltype_colname = "nichenet",
                         condition_colname = "samplecondition", condition_oi = condition_oi,
                         features = union(expressed_ligands, expressed_receptors))

# Average expression information - only for LCMV condition
expression_info <- get_exprs_avg(seuratObj, "nichenet", condition_colname = "samplecondition", condition_oi = condition_oi)

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


rotated_dotplot <- DotPlot(seuratObj %>% subset(nichenet %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") +
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




#adapted for main ligands in PAH and PHPF


run_nichenet <- function(condition_oi, top_n_ligands = 50) {
  
  condition_reference <- "Donor"
  
  # Receiver DE
  seurat_obj_receiver <- subset(seuratObj, idents = receiver)
  seurat_obj_receiver <- SetIdent(seurat_obj_receiver,
                                  value = seurat_obj_receiver[["samplecondition", drop = TRUE]])
  
  DE_table_receiver <- FindMarkers(
    object = seurat_obj_receiver,
    ident.1 = condition_oi,
    ident.2 = condition_reference,
    min.pct = 0.10
  ) %>% rownames_to_column("gene")
  
  geneset_oi <- DE_table_receiver %>%
    dplyr::filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))
  
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ) %>%
    arrange(desc(aupr_corrected))
  
  top_ligands <- ligand_activities %>%
    slice_head(n = top_n_ligands)
  
  return(top_ligands)
}

ipah_results <- run_nichenet("IPAH", 50)
pf_results   <- run_nichenet("IPF", 50)

common <- intersect(ipah_results$test_ligand, pf_results$test_ligand)
common
length(common)

ipah_unique <- setdiff(ipah_results$test_ligand, pf_results$test_ligand)
pf_unique   <- setdiff(pf_results$test_ligand, ipah_results$test_ligand)

comparison <- full_join(
  ipah_results %>% select(test_ligand, aupr_corrected) %>% rename(IPAH = aupr_corrected),
  pf_results   %>% select(test_ligand, aupr_corrected) %>% rename(PF   = aupr_corrected),
  by = "test_ligand"
) %>% arrange(desc(IPAH))

comparison

library(ggvenn)

ggvenn(list(IPAH = ipah_results$test_ligand,
            PF   = pf_results$test_ligand))

#venn too suspicious use other metric
comparison <- full_join(
  ipah_results %>% select(test_ligand, aupr_corrected) %>% rename(IPAH = aupr_corrected),
  pf_results %>% select(test_ligand, aupr_corrected) %>% rename(PF = aupr_corrected),
  by = "test_ligand"
) %>%
  mutate(delta = IPAH - PF)

ggplot(comparison, aes(IPAH, PF)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic()


strong_ligands <- comparison %>%
  filter(IPAH > 0.05 | PF > 0.05)

specific <- comparison %>%
  filter(abs(delta) > 0.02)









#construct grounds for comparison

ipah_results
pf_results


ipah_rank <- ipah_results %>%
  mutate(rank_IPAH = rank(-aupr_corrected)) %>%
  select(test_ligand, rank_IPAH)


pf_rank <- pf_results %>%
  mutate(rank_PF = rank(-aupr_corrected)) %>%
  select(test_ligand, rank_PF)


comparison <- full_join(ipah_rank, pf_rank, by = "test_ligand")

comparison <- comparison %>%
  mutate(rank_diff = rank_IPAH - rank_PF)

ipah_specific <- comparison %>%
  filter(rank_IPAH <= 20 & rank_PF > 40)


pf_specific <- comparison %>%
  filter(rank_PF <= 20 & rank_IPAH > 40)


shared <- comparison %>%
  filter(rank_IPAH <= 30 & rank_PF <= 30)


ipah_specific_lr <- ipah_specific %>%
  rename(ligand = test_ligand) %>%
  left_join(lr_network_renamed, by = "ligand")


#cross validation
#Do CellChat ligands score high in NicheNet?
run_nichenet <- function(condition_oi) {
  
  
  condition_reference <- "Donor"
  
  
  seurat_obj_receiver <- subset(seuratObj, idents = receiver)
  seurat_obj_receiver <- SetIdent(seurat_obj_receiver,
                                  value = seurat_obj_receiver[["samplecondition", drop = TRUE]])
  
  
  DE_table_receiver <- FindMarkers(
    object = seurat_obj_receiver,
    ident.1 = condition_oi,
    ident.2 = condition_reference,
    min.pct = 0.10
  ) %>% rownames_to_column("gene")
  
  
  geneset_oi <- DE_table_receiver %>%
    dplyr::filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))
  
  
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ) %>%
    arrange(desc(aupr_corrected))
  
  
  return(ligand_activities) # ← full table
}

expressed_genes_sender <- get_expressed_genes(
  sender_celltypes,
  seuratObj,
  pct = 0.10
)

expressed_genes_receiver <- get_expressed_genes(
  receiver,
  seuratObj,
  pct = 0.10
)

potential_ligands <- lr_network %>%
  dplyr::filter(from %in% expressed_genes_sender &
                  to %in% expressed_genes_receiver) %>%
  dplyr::pull(from) %>%
  unique()



ipah_results <- run_nichenet("IPAH")
pf_results <- run_nichenet("IPF")


validate_set <- function(nichenet_df, ligand_set, label) {
  
  
  df <- nichenet_df %>%
    mutate(in_set = test_ligand %in% ligand_set)
  
  
  test <- wilcox.test(
    aupr_corrected ~ in_set,
    data = df
  )
  
  
  summary <- df %>%
    group_by(in_set) %>%
    summarise(
      mean_score = mean(aupr_corrected),
      median_score = median(aupr_corrected),
      n = n()
    )
  
  
  list(label = label,
       stats = summary,
       p_value = test$p.value)
}


val_ipah <- validate_set(ipah_results, PAH.exc.ligands, "IPAH exclusive")
val_ipah

val_pf <- validate_set(pf_results, PF.exc.ligands, "PF exclusive")
val_pf

val_shared_ipah <- validate_set(ipah_results, overlap.ligands, "Shared in IPAH")
val_shared_ipah 
val_shared_pf   <- validate_set(pf_results, overlap.ligands, "Shared in PF")
val_shared_pf




mean_score (in_set = TRUE)  >  mean_score (FALSE)
p_value = 0.003

#results suggests that  



#pah vs pf nichenet

run_nichenet <- function(condition_oi, reference) {
  
  # Receiver subset
  seurat_obj_receiver <- subset(seuratObj, idents = receiver)
  seurat_obj_receiver <- SetIdent(
    seurat_obj_receiver,
    value = seurat_obj_receiver[["samplecondition", drop = TRUE]]
  )
  
  # Differential expression: disease vs disease
  DE_table_receiver <- FindMarkers(
    object = seurat_obj_receiver,
    ident.1 = condition_oi,
    ident.2 = reference,
    min.pct = 0.10
  ) %>%
    rownames_to_column("gene")
  
  # Gene set of interest
  geneset_oi <- DE_table_receiver %>%
    dplyr::filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))
  
  # NicheNet ligand activity
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ) %>%
    arrange(desc(aupr_corrected))
  
  return(ligand_activities)
}



pf_vs_ipah <- run_nichenet("IPF", "IPAH")
ipah_vs_pf <- run_nichenet("IPAH", "IPF")

pf_vs_ipah <- pf_vs_ipah %>%
  mutate(direction = "IPF > IPAH")

ipah_vs_pf <- ipah_vs_pf %>%
  mutate(direction = "IPAH > IPF")


pf_specific <- pf_vs_ipah %>% slice_head(n = 50)
ipah_specific <- ipah_vs_pf %>% slice_head(n = 50)


head(pf_vs_ipah)


val_pf_specific <- validate_set(
  pf_vs_ipah,
  PF.exc.ligands,
  "PF exclusive vs IPAH"
)

val_pf_specific



val_ipah_specific <- validate_set(
  ipah_vs_pf,
  PAH.exc.ligands,
  "IPAH exclusive vs PF"
)

val_ipah_specific





val_shared_pf <- validate_set(
  pf_vs_ipah,
  overlap.ligands,
  "Shared ligands (PF)"
)

val_shared_ipah <- validate_set(
  ipah_vs_pf,
  overlap.ligands,
  "Shared ligands (IPAH)"
)


plot_validation <- function(nichenet_df, ligand_set, title) {
  
  
  nichenet_df %>%
    mutate(group = ifelse(test_ligand %in% ligand_set,
                          "CellChat ligands", "Background")) %>%
    ggplot(aes(group, aupr_corrected, fill = group)) +
    geom_boxplot() +
    theme_classic() +
    ggtitle(title)
}


plot_validation(pf_vs_ipah, PF.exc.ligands, "PF vs IPAH validation")
plot_validation(ipah_vs_pf, PAH.exc.ligands, "IPAH vs PF validation")

#graph1
pf_vs_ipah %>%
  mutate(group =
           case_when(
             test_ligand %in% PF.exc.ligands ~ "PF-exclusive",
             test_ligand %in% PAH.exc.ligands ~ "PAH-exclusive",
             test_ligand %in% overlap.ligands ~ "Shared",
             TRUE ~ "Background"
           ),
         rank = rank(-aupr_corrected)) %>%
  ggplot(aes(rank, aupr_corrected, color = group)) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  labs(x = "Ligand rank (NicheNet)",
       y = "Ligand activity score")


pf_vs_ipah %>%
  mutate(group =
           case_when(
             test_ligand %in% PF.exc.ligands ~ "PF-exclusive",
             test_ligand %in% PAH.exc.ligands ~ "PAH-exclusive",
             test_ligand %in% overlap.ligands ~ "Shared",
             TRUE ~ "Background"
           ),
         group = factor(group,
                        levels = c("PF-exclusive",
                                   "PAH-exclusive",
                                   "Shared",
                                   "Background")),
         rank = rank(-aupr_corrected)) %>%
  ggplot(aes(rank, aupr_corrected, color = group)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "PF-exclusive" = "#377eb8",
    "PAH-exclusive" = "#e67e22",
    "Shared" = "#7fd34e",
    "Background" = "grey80"
  )) +
  theme_classic() +
  labs(x = "Ligand rank (NicheNet)",
       y = "Ligand activity score",
       color = "Ligand class")


df <- pf_vs_ipah %>%
  mutate(group =
           case_when(
             test_ligand %in% PF.exc.ligands ~ "PF-exclusive",
             test_ligand %in% PAH.exc.ligands ~ "PAH-exclusive",
             test_ligand %in% overlap.ligands ~ "Shared",
             TRUE ~ "Background"
           ),
         rank = rank(-aupr_corrected))

ggplot() +
  
  # background first
  geom_point(data = subset(df, group == "Background"),
             aes(rank, aupr_corrected),
             color = "grey85",
             alpha = 0.4,
             size = 1.5) +
  
  # highlighted ligands on top
  geom_point(data = subset(df, group != "Background"),
             aes(rank, aupr_corrected, color = group),
             size = 2.5,
             alpha = 0.95) +
  
  scale_color_manual(values = c(
    "PF-exclusive" = "#377eb8",
    "PAH-exclusive" = "#e67e22",
    "Shared" = "#7fd34e"
  )) +
  
  theme_classic() +
  labs(x = "Ligand rank (NicheNet)",
       y = "Ligand activity score",
       color = "Ligand class")



df <- pf_vs_ipah %>%
  mutate(group =
           case_when(
             test_ligand %in% PF.exc.ligands ~ "PF-exclusive",
             test_ligand %in% PAH.exc.ligands ~ "PAH-exclusive",
             test_ligand %in% overlap.ligands ~ "Shared",
             TRUE ~ "Background"
           ),
         group = factor(group,
                        levels = c("Background",
                                   "PF-exclusive",
                                   "PAH-exclusive",
                                   "Shared")),
         rank = rank(-aupr_corrected))

ggplot() +
  
  # background layer (mapped to group!)
  geom_point(data = subset(df, group == "Background"),
             aes(rank, aupr_corrected, color = group),
             alpha = 0.3,
             size = 1.5) +
  
  # highlighted ligands
  geom_point(data = subset(df, group != "Background"),
             aes(rank, aupr_corrected, color = group),
             size = 2.5,
             alpha = 0.95) +
  
  scale_color_manual(values = c(
    "Background" = "grey80",
    "PF-exclusive" = "#377eb8",
    "PAH-exclusive" = "#e67e22",
    "Shared" = "#7fd34e"
  )) +
  
  theme_classic() +
  labs(x = "Ligand rank (NicheNet)",
       y = "Ligand activity score",
       color = "Ligand class")



#graph2
pf_vs_ipah %>%
  mutate(group =
           case_when(
             test_ligand %in% PF.exc.ligands ~ "PF-exclusive",
             test_ligand %in% PAH.exc.ligands ~ "PAH-exclusive",
             test_ligand %in% overlap.ligands ~ "Shared",
             TRUE ~ "Background"
           )) %>%
  ggplot(aes(group, aupr_corrected, fill = group)) +
  geom_boxplot() +
  theme_classic() +
  labs(y = "Ligand activity score")

#graph3
library(tidygraph)
library(ggraph)
library(ggplot2)


nodes <- data.frame(
  name = c(PF.exc.ligands, PAH.exc.ligands, "Shared endothelial program"),
  group = c(rep("PF", length(PF.exc.ligands)),
            rep("PAH", length(PAH.exc.ligands)),
            "Endpoint")
)


edges <- data.frame(
  from = c(PF.exc.ligands, PAH.exc.ligands),
  to = "Shared endothelial program"
)


graph <- tbl_graph(nodes = nodes, edges = edges)


ggraph(graph, layout = "linear") +
  geom_edge_link(alpha = 0.6) +
  geom_node_point(aes(color = group), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(
    PF = "#377eb8",
    PAH = "#e41a1c",
    Endpoint = "black"
  )) +
  theme_void()


library(ggalluvial)
library(dplyr)


flow_data <- data.frame(
  disease = c(rep("PF", length(PF.exc.ligands)),
              rep("PAH", length(PAH.exc.ligands))),
  ligand = c(PF.exc.ligands, PAH.exc.ligands),
  endpoint = "Common pathological vascular endpoint"
)


ggplot(flow_data,
       aes(axis1 = disease,
           axis2 = ligand,
           axis3 = endpoint)) +
  geom_alluvium(aes(fill = disease), width = 0.3) +
  geom_stratum(width = 0.3, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Disease", "Ligand", "Endpoint"),
                   expand = c(.1, .1)) +
  theme_void()







#Cell chat

ECstory.fil$cellchat.ident <- ""
ECstory.fil$ec.cluster.sub <- ECstory.subset$seurat_clusters
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "27")] <- "Fibroblasts"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "0" )] <- "Fibroblasts"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "24")] <- "Fibroblasts"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "4" )] <- "SMC"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "14")] <- "SMC"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "16")] <- "SMC" #pericytes seem to be a subset of these
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "9" )] <- "SMC"    
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "21")] <- "B-cells"  
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "23")] <- "Platelets" 
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "20")] <- "Granulocytes"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "22")] <- "Erythrocytes"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "13")] <- "Mono/Macs"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "7" )] <- "Mono/Macs"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "2" )] <- "Mono/Macs"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "25")] <- "Dendritic cells"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "19")] <- "Dendritic cells"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "18")] <- "Dendritic cells"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "11")] <- "Mono/Macs"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "5")] <- "Mono/Macs"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "17")] <- "Mono/Macs"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "15")] <- "Epithelial cells"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "26")] <- "Epithelial cells"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "8")] <- "NK cells"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "1")] <- "T-cells"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "3")] <- "T-cells"
ECstory.fil$cellchat.ident[(ECstory.fil$seurat_clusters == "10")] <- "T-cells"
ECstory.fil$cellchat.ident[(ECstory.fil$ec.cluster.sub == "0")] <- "Endothelial"
ECstory.fil$cellchat.ident[(ECstory.fil$ec.cluster.sub == "1")] <- "Endothelial"
ECstory.fil$cellchat.ident[(ECstory.fil$ec.cluster.sub == "2")] <- "Endothelial"
ECstory.fil$cellchat.ident[(ECstory.fil$ec.cluster.sub == "3")] <- "EC3"
ECstory.fil$cellchat.ident[(ECstory.fil$ec.cluster.sub == "4")] <- "EC4"

DimPlot(ECstory.fil, group.by = "cellchat.ident")

ECstory.cellchat.fil2 <- subset(ECstory.fil, cellchat.ident %!in% c("Platelets", "Endo MT", "Erythrocytes", "EC3", "EC4"))


DimPlot(ECstory.cellchat.fil2, group.by = "cellchat.ident")





Donor <- subset(ECstory.cellchat.fil2,  samplecondition %in% c("Donor"))
Donor <- JoinLayers(Donor)
PF <- subset(ECstory.cellchat.fil2,  samplecondition %in% c("IPF"))
PF <- JoinLayers(PF)
PAH <- subset(ECstory.cellchat.fil2,  samplecondition %in% c("IPAH"))
PAH <- JoinLayers(PAH)


DimPlot(Donor, group.by = "cellchat.ident")

data.input <- GetAssayData(PF, assay = "RNA", slot = "data") # normalized data matrix
cellchat <- createCellChat(object = PF, group.by = "cellchat.ident", assay = "RNA")
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data #PPI.mouse for protein prtoein interactions
showDatabaseCategory(CellChatDB)

cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, nboot = 20)


cellchat <- filterCommunication(cellchat, min.cells = 5)

df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

pal3 <- c("#7ED7D1FF", "#683B79FF", "#F8A19FFF", "#16FF32FF", "#FA0087FF", "#325A9BFF", "#2ED9FFFF","#f39d10","#038637", "#F6222EFF")


names(pal3) <- unique(Donor$cellchat.ident)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count,top=0.3,color.use =pal3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=0.3,color.use =pal3,vertex.weight = groupSize, weight.scale = T, title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)

netVisual_circle(cellchat@net$count,top=0.5,color.use =pal3,sources.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=0.5,color.use =pal3,sources.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)

netVisual_circle(cellchat@net$count,top=1,color.use =pal3,targets.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Number of interactions", vertex.label.cex=0.01,label.edge=T)
netVisual_circle(cellchat@net$weight,top=1,color.use =pal3,targets.use = 3 ,vertex.weight = groupSize, weight.scale = T,title.name = "Interaction weights/strength", vertex.label.cex=0.01,label.edge=T)


#netVisual_bubble(cellchat, sources.use = 3, targets.use = ECstory.cellchat.fil$cellchat.ident, remove.isolate = TRUE, sort.by.source = T)


#cellchat_Donor_allcells <- cellchat
#cellchat_PAH_allcells <- cellchat
#cellchat_PF_allcells <- cellchat


object.list <- list(Donor = cellchat_Donor_allcells, PF = cellchat_PF_allcells)

object.list <- list(Donor = cellchat_Donor_allcells, PAH = cellchat_PAH_allcells)

cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

cellchat.merge@meta$cellchat.ident <- factor(cellchat.merge@meta$cellchat.ident, 
                                              levels = c("Endothelial",
                                                         "Fibroblasts", "SMC",
                                                         "Epithelial cells",
                                                         "Mono/Macs", "NK cells", "T-cells", "B-cells", "Granulocytes", "Dendritic cells"))

levels(cellchat.merge@meta$cellchat.ident)
unique(cellchat.merge@meta$cellchat.ident)
names(pal3) <- unique(cellchat.merge@meta$cellchat.ident)

par(mfrow=c(1,1))
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("count"), color.use = pal3)
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("weight"), color.use = pal3)
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("count"))
netVisual_heatmap(cellchat.merge, signaling = NULL,  measure = c("weight"))


netVisual_bubble(cellchat.merge, sources.use = 3, targets.use = NULL ,  comparison = c(1, 2), title.name = "Increased signaling in Disease", 
                 angle.x = 45, remove.isolate = F,color.text=c("black","#0491b3"), font.size = 8,color.heatmap = c("Spectral"), max.dataset = 2 )

df <- subsetCommunication(cellchat.merge)
df <- do.call(rbind, df)
source.names <- c("Endothelial")


df <- df[df$source %in% source.names, ]

comparison = c(1, 2)
max.dataset = 2

df <- df[df$prob > 0, ]

genes <- unique(unlist(strsplit(c(df$ligand, df$receptor), "\\|")))
genes <- genes[genes != ""]



#ligands only
df <- subsetCommunication(cellchat.merge)
df <- do.call(rbind, df)


df$source <- as.character(df$source)


df <- df[df$source %in% c("Endothelial"), ]
df <- df[df$prob > 0, ]

ligands <- unique(unlist(strsplit(df$ligand, "\\|")))
ligands <- ligands[ligands != ""]
ligands <- unique(toupper(ligands))
ligands







#spec list 

ligand.list <- c(
  "TGFB1",
  "NRG1",
  "MDK",
  "LAMB2",
  "HLA-DQA1",
  "GRN",
  "GAS6",
  "APOE"
)


ligand_list_both <- c(
  "SPP1",
  "PTPRM",
  "PPIA",
  "POSTN",
  "PECAM1",
  "MIF",
  "JAM2",
  "IGFBP3",
  "ICAM2",
  "HLA-E",
  "HLA-DRB5",
  "HLA-DRB1",
  "HLA-DRA",
  "GJA1",
  "FN1",
  "EFNB2",
  "COL4A2",
  "LIPA", # from "CHOLESTEROL -LIPA"
  "CD99",
  "CD55",
  "CD34",
  "CCL14",
  "BMP6",
  "APP",
  "ANXA1",
  "ITGAV",
  "HLA-DQB1",
  "HLA-DPB1",
  "HLA-DPA1",
  "HLA-DMA"
)

ligand_list_PF <- c(
  "SEMA6A",
  "PTGES3", # from "PGE2-PTGES3"
  "NAMPT",
  "LAMC1",
  "LAMB1",
  "ICAM1",
  "HLA-F",
  "HLA-C",
  "HLA-B",
  "HLA-A",
  "CXCL2",
  "CXCL12",
  "COL4A1",
  "CCL2"
)


df <- subsetCommunication(cellchat.merge)
df <- do.call(rbind, df)

df$source <- as.character(df$source)

df <- df[df$source %in% c("Endothelial"), ]

ligand_list_PF <- toupper(ligand.list)
df <- df[
  sapply(df$ligand, function(x)
    any(strsplit(x, "\\|")[[1]] %in% ligand.list)
  ),
]

unique(df$ligand)

pairLR.use <- unique(df$interaction_name)
pairLR.use <- data.frame(
  interaction_name = unique(df$interaction_name),
  stringsAsFactors = FALSE
)


netVisual_bubble(
  cellchat.merge,
  sources.use = c("Endothelial"),
  targets.use = NULL,
  comparison = c(1, 2),
  max.dataset = 2,
  pairLR.use = pairLR.use,
  title.name = "PF-increased interactions (selected ligands)",
  angle.x = 45,
  remove.isolate = TRUE,
  font.size = 8
)



#
overlap <- intersect(PF.ligands, PAH.ligands)
overlap.ligands <- overlap

PF.exc.ligands <- setdiff(PF.ligands, PAH.ligands)

PAH.exc.ligands <- setdiff(PAH.ligands, PF.ligands)


#receptors
df <- subsetCommunication(cellchat.merge)
df <- do.call(rbind, df)


df$target <- as.character(df$target)


df <- df[df$target %in% c("Endothelial"), ]
df <- df[df$prob > 0, ]


receptors <- unique(unlist(strsplit(df$receptor, "\\|")))
receptors <- receptors[receptors != ""]
receptors <- unique(toupper(receptors))


receptors

#PF.receptors <- receptors
#PAH.receptors <- receptors

overlap <- intersect(PF.receptors, PAH.receptors)
overlap.receptors <- overlap

PF.exc.receptors <- setdiff(PF.receptors, PAH.receptors)

PAH.exc.receptors <- setdiff(PAH.receptors, PF.receptors)


#graph1
netVisual_diffInteraction(
  cellchat.merge,
  weight.scale = TRUE,
  measure = "weight",
  label.edge = FALSE
)



#graph2
library(ggvenn)
library(ggVennDiagram)


ggvenn(
  list(
    PF = PF.ligands,
    PAH = PAH.ligands
  ),
  fill_color = c("#377eb8","darkorange"),
  stroke_size = 0.5
)

ggVennDiagram(
  list(PF = PF.ligands, PAH = PAH.ligands),
  label_alpha = 1, # box opacity (1 = solid)
  label_color = "white", # brighter text
  label_size = 5 # optional: bigger text
) +
  scale_fill_gradientn(
    colours = c("darkorange","#377eb8",  "lightgreen"),
    values = scales::rescale(c(5, 10, 34)),
    limits = c(5, 34)
  )




#graph
library(tidygraph)
library(ggraph)
library(dplyr)


edges <- data.frame(
  from = c(PF.exc.ligands, PAH.exc.ligands),
  to = "Endothelial",
  disease = c(rep("PF", length(PF.exc.ligands)),
              rep("PAH", length(PAH.exc.ligands)))
)


nodes <- data.frame(
  name = unique(c(edges$from, edges$to))
) %>%
  mutate(type = name == "Endothelial") # bipartite partition


graph <- tbl_graph(nodes = nodes, edges = edges)


ggraph(graph, layout = "bipartite") +
  geom_edge_link(aes(color = disease), alpha = 0.8) +
  geom_node_point(size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(PF = "#377eb8", PAH = "#e41a1c")) +
  theme_void()


#graph
library(tidygraph)
library(ggraph)
library(dplyr)

edges <- data.frame(
  from = c(PF.exc.ligands,
           overlap.ligands,
           overlap.ligands,
           PAH.exc.ligands),
  to = "Endothelial",
  group = c(rep("PF", length(PF.exc.ligands)),
            rep("Shared", length(overlap.ligands)),
            rep("Shared", length(overlap.ligands)),
            rep("PAH", length(PAH.exc.ligands)))
)

nodes <- data.frame(
  name = unique(c(edges$from, edges$to))
)

graph <- tbl_graph(nodes = nodes, edges = edges)

ggraph(graph, layout = "linear") +
  geom_edge_link(aes(color = group), alpha = 0.7) +
  geom_node_point(size = 3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(
    PF = "#377eb8",
    Shared = "lightgreen",
    PAH = "darkorange"
  )) +
  theme_void()

library(ggalluvial)
library(dplyr)


flow_data <- bind_rows(
  data.frame(disease = "PF", class = "PF-exclusive", ligand = PF.exc.ligands),
  data.frame(disease = "PF", class = "Shared", ligand = overlap.ligands),
  data.frame(disease = "PAH", class = "Shared", ligand = overlap.ligands),
  data.frame(disease = "PAH", class = "PAH-exclusive", ligand = PAH.exc.ligands)
)


flow_data$endpoint <- "" #Convergent pulmonary vascular pathology


ggplot(flow_data,
       aes(axis1 = disease,
           axis2 = class,
           axis3 = endpoint)) +
  geom_alluvium(aes(fill = class), width = 0.3) +
  geom_stratum(width = 0.3, fill = "grey85", color = "black") +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Disease", "Ligand class", "Endpoint"),
                   expand = c(.1, .1)) +
  scale_fill_manual(values = c(
    "PF-exclusive" = "#377eb8", #"PF-exclusive"
    "Shared" = "lightgreen", #"Shared"
    "PAH-exclusive" = "darkorange" #"PAH-exclusive"
  )) +
  theme_void()

#clusterprofiler
library(clusterProfiler)
library(pheatmap)
library(dplyr)
library(org.Hs.eg.db)


pf_go <- enrichGO(PF.exc.ligands, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
pah_go <- enrichGO(PAH.exc.ligands, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
shared_go <- enrichGO(overlap.ligands, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")


merge_go <- bind_rows(
  pf_go@result %>% mutate(group = "PF"),
  pah_go@result %>% mutate(group = "PAH"),
  shared_go@result %>% mutate(group = "Shared")
)


mat <- merge_go %>%
  dplyr::select(Description, group, p.adjust) %>%
  mutate(score = -log10(p.adjust)) %>%
  tidyr::pivot_wider(names_from = group, values_from = score)

mat[is.na(mat)] <- 0
mat[is.infinite(mat)] <- max(mat[is.finite(mat)], na.rm = TRUE)


rownames(mat) <- mat$Description
mat <- as.matrix(mat[,-1])


pheatmap(mat, labels_row = rownames(mat))


#ligand fingerprint

ligands <- unique(c(PAH.exc.ligands,
                    overlap.ligands,
                    PF.exc.ligands))

ligand_df <- data.frame(
  ligand = ligands,
  PAH = ligands %in% PAH.exc.ligands,
  Shared = ligands %in% overlap.ligands,
  PHPF = ligands %in% PF.exc.ligands
)

library(tidyr)
ligand_long <- pivot_longer(ligand_df,
                            -ligand,
                            names_to = "group",
                            values_to = "present")

ligand_long$ligand <- factor(
  ligand_long$ligand,
  levels = c(overlap.ligands,
             PF.exc.ligands,
             PAH.exc.ligands)
)

library(ggplot2)

ggplot(ligand_long,
       aes(group, ligand, fill = present)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("FALSE" = "grey95",
                               "TRUE" = "#e41a1c")) +
  theme_classic() +
  labs(x = "", y = "Ligands") +
  theme(axis.text.y = element_text(size = 7))





#-#--#-#
mat_df <- merge_go %>%
  dplyr::select(Description, group, p.adjust) %>%
  dplyr::mutate(score = -log10(p.adjust)) %>%
  tidyr::pivot_wider(names_from = group, values_from = score)


mat_df$Description <- make.unique(mat_df$Description)

mat_num <- mat_df[, -1]
mat_num <- as.data.frame(lapply(mat_num, as.numeric))

mat_num[is.na(mat_num)] <- 0
mat_num[is.infinite(mat_num)] <- max(mat_num[is.finite(mat_num)], na.rm = TRUE)

mat <- as.matrix(mat_num)


pheatmap(mat,
         color = colorRampPalette(c("white", "darkred"))(50),
         border_color = NA,
         fontsize_row = 8,
         fontsize_col = 10)


top_terms <- merge_go %>%
  group_by(group) %>%
  slice_min(p.adjust, n = 8) %>%   # top 8 per group
  pull(Description) %>%
  unique()

merge_go2 <- merge_go %>% filter(Description %in% top_terms)

shorten <- function(x, n = 40) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n), "..."), x)
}

merge_go2$Description <- shorten(merge_go2$Description)


mat_df <- merge_go2 %>%
  dplyr::select(Description, group, p.adjust) %>%
  dplyr::mutate(score = -log10(p.adjust)) %>%
  tidyr::pivot_wider(names_from = group, values_from = score)

mat_df$Description <- make.unique(mat_df$Description)



mat_num <- as.matrix(mat_df[,-1])
mat_num[is.na(mat_num)] <- 0

mat_df$id <- paste0("GO_", seq_len(nrow(mat_df)))
pretty_labels <- wrap(mat_df$id)

rownames(mat_num) <- mat_df$Description

pheatmap(mat_num,
         labels_row = pretty_labels,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#b2182b"))(50),
         fontsize_row = 9,
         fontsize_col = 11,
         border_color = NA)

pheatmap(mat_num,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#b2182b"))(50),
         fontsize_row = 9,
         fontsize_col = 11,
         border_color = NA,
         legend_breaks = c(0, 5, 10),
         legend_labels = c("0", "5", "10"),
         main = "GO enrichment (-log10 adj. p)")




wrap <- function(x, width = 35) {
  vapply(strwrap(x, width = width, simplify = T),
         paste, collapse = "\n", FUN.VALUE = character(1))
}


rownames(mat_num) <- wrap(rownames(mat_num))

pheatmap(mat_num,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#e41a1c"))(100),
         fontsize_row = 7,
         cellheight = 12,
         cellwidth = 25,
         main = "GO enrichment (-log10 adj. p)")


library(clusterProfiler)


pf_go <- simplify(pf_go, cutoff = 0.7, by = "p.adjust",select_fun = min)
pah_go <- simplify(pah_go, cutoff = 0.7, by = "p.adjust",select_fun = min)
shared_go <- simplify(shared_go, cutoff = 0.7, by = "p.adjust",select_fun = min)


#-------
pf_top <- pf_go@result %>% arrange(p.adjust) %>% head(10)
pah_top <- pah_go@result %>% arrange(p.adjust) %>% head(10)
shared_top <- shared_go@result %>% arrange(p.adjust) %>% head(10)

merge_go <- bind_rows(
  pf_top %>% mutate(group = "PF"),
  pah_top %>% mutate(group = "PAH"),
  shared_top %>% mutate(group = "Shared")
)

mat_df <- merge_go %>%
  dplyr::select(Description, group, p.adjust) %>%
  dplyr::mutate(score = -log10(p.adjust)) %>%
  tidyr::pivot_wider(names_from = group, values_from = score)


mat_num <- as.matrix(mat_df[, -1])
rownames(mat_num) <- wrap(mat_df$Description)
mat_num[is.na(mat_num)] <- 0


#------------------------


mat_df$id <- paste0("GO_", seq_len(nrow(mat_df)))
rownames(mat_num) <- mat_df$id

pretty_labels <- wrap(mat_df$Description)
pretty_labels <- wrap(mat_df$id)

pheatmap(mat_num,
         labels_row = pretty_labels,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#e41a1c"))(100),
         fontsize_row = 6,
         fontsize_col = 10,
         cellheight = 14,
         cellwidth = 30,
         border_color = NA,
         main = "GO enrichment (-log10 adj. p)")
#-##--#-
library(dplyr)
library(tidyr)

extract_go <- function(go, label) {
  go@result %>%
    dplyr::select(Description, p.adjust) %>%
    mutate(group = label,
           score = -log10(p.adjust))
}

pf_tbl <- extract_go(pf_go, "PF")
pah_tbl <- extract_go(pah_go, "PAH")
shared_tbl <- extract_go(shared_go, "Shared")

merged <- bind_rows(pf_tbl, pah_tbl, shared_tbl)



top_terms <- merged %>%
  group_by(Description) %>%
  summarise(max_score = max(score, na.rm = TRUE)) %>%
  arrange(desc(max_score)) %>%
  slice_head(n = 15) %>%
  pull(Description)

merged <- merged %>% filter(Description %in% top_terms)


mat <- merged %>%
  dplyr::select(Description, group, score) %>%
  pivot_wider(names_from = group, values_from = score) %>%
  replace(is.na(.), 0)

rownames(mat) <- mat$Description
mat <- as.matrix(mat[,-1])


library(pheatmap)

pheatmap(
  mat,
  color = colorRampPalette(c("white", "#1f78b4"))(50),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Ligand functional enrichment"
)

pheatmap(
  mat,
  color = colorRampPalette(c("white", "darkred"))(50),
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = 0
)



mat_df <- merge_go2 %>%
  dplyr::select(Description, group, p.adjust) %>%
  dplyr::mutate(score = -log10(p.adjust)) %>%
  dplyr::select(-p.adjust) %>%     # ← remove raw column
  tidyr::pivot_wider(names_from = group, values_from = score)

mat_num <- as.matrix(mat_df[,-1])
mat_num[is.na(mat_num)] <- 0


pheatmap(mat_num,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#e41a1c"))(100),
         fontsize_row = 6,
         fontsize_col = 10,
         cellheight = 14,
         cellwidth = 30,
         border_color = NA,
         main = "GO enrichment (-log10 adj. p)")


mat_num <- mat_num[rowSums(mat_num) > 1, ]

pheatmap(mat_num,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#e41a1c"))(100),
         fontsize_row = 6,
         fontsize_col = 10,
         cellheight = 14,
         cellwidth = 60, # ← important
         border_color = NA,
         main = "GO enrichment (-log10 adj. p)")







