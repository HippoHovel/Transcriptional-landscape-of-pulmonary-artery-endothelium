#cell cycle

ECstory.subset.fil$namedclusters <- ""
ECstory.subset.fil$namedclusters[ECstory.subset.fil$seurat_clusters == "0"] <- "Immunomodulatory"
ECstory.subset.fil$namedclusters[ECstory.subset.fil$seurat_clusters == "1"] <- "Vascular tone-modulatory"
ECstory.subset.fil$namedclusters[ECstory.subset.fil$seurat_clusters == "2"] <- "Vascular plasticity modulatory"



ECstory.subset.fil$clusternames <- ""
ECstory.subset.fil$clusternames[ECstory.subset.fil$seurat_clusters == "0"] <- "Immunomodulatory"
ECstory.subset.fil$clusternames[ECstory.subset.fil$seurat_clusters == "1"] <- "Vascular tone-modulatory"
ECstory.subset.fil$clusternames[ECstory.subset.fil$seurat_clusters == "2"] <- "Vascular plasticity modulatory"
#ECstory.subset.fil <- subset(ECstory.subset.fil, sample.ident %!in% c("PAH4"))
DimPlot(ECstory.subset.fil, group.by = "namedclusters")
DimPlot(ECstory.subset.fil, group.by = "clusternames")

s.genes <- cc.genes$genenames[(cc.genes$phase == "S")]
g2m.genes <- cc.genes$genenames[(cc.genes$phase == "G2/M")]

ECstory.subset.fil <- FindVariableFeatures(ECstory.subset.fil, selection.method = "vst")
ECstory.subset.fil <- CellCycleScoring(ECstory.subset.fil, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=TRUE)

DimPlot(ECstory.subset.fil, label=F,group.by = "Phase", cols=myphasecol)


DimPlot(ECstory.subset.fil, label=F,split.by = "Phase", group.by = "seurat_clusters", cols=c("#b20505","#0312ad","#03ad46"))


DimPlot(ECstory.subset.fil, split.by = "sample.ident", group.by = "namedclusters")

ECstory.subset.fil.Donor <- subset(ECstory.subset.fil, samplecondition %in% c("Donor"))
ECstory.subset.fil.PAH <- subset(ECstory.subset.fil, samplecondition %in% c("IPAH"))
ECstory.subset.fil.PF <- subset(ECstory.subset.fil, samplecondition %in% c("IPF"))






#Donor
meta_df_Donor <- ECstory.subset.fil.Donor@meta.data

#Count cells per clusternames × samplecondition × sample.ident
cell_counts_Donor <- meta_df_Donor %>%
  count(clusternames, samplecondition, sample.ident, name = "cell_count")

#Calculate percentage per sample.ident
cell_counts_Donor <- cell_counts_Donor %>%
  group_by(sample.ident) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

#View results
print(cell_counts_Donor, n = 35)


summary_stats_Donor <- cell_counts_Donor %>%
  group_by(clusternames) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )

print(summary_stats_Donor)


summary_stats_Donor <- cell_counts_Donor %>%
  group_by(clusternames) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )


ggplot(summary_stats_Donor, aes(x="", y=median_percent, fill=clusternames)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=c("#b20505","#0312ad","#03ad46"))



ggplot(cell_counts_Donor, aes(x = clusternames, y = percent)) +
  geom_boxplot(fill = c("#b20505","#0312ad","#03ad46"), color = "black", outlier.shape = 21, outlier.fill = "red") +
  labs(
    title = "Distribution of Cluster Percentages Across Samples",
    x = "Cluster Name",
    y = "Percentage of Cells per Sample"
  ) +
  geom_jitter(
    aes(color = clusternames),
    width = 0.15,
    size = 2.5,
    alpha = 0.9
  ) +
  scale_fill_manual(values = c("#CC4B4B","#454FD1","#57D487")) +
  scale_color_manual(values = c("#CC4B4B","#454FD1","#57D487")) +
  labs(
    title = "Cluster Proportions Across Donor Samples",
    x = "Cluster",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
    )



stats_donor <- cell_counts_Donor %>%
  group_by(clusternames) %>%
  summarise(
    median = median(percent),
    q1 = quantile(percent, 0.25),
    q3 = quantile(percent, 0.75),
    ymax = max(percent),
    .groups = "drop"
  ) %>%
  mutate(label_y = 55)   # fixed offset

ggplot(cell_counts_Donor, aes(x = clusternames, y = percent)) +
  geom_boxplot(fill = c("#b20505","#03ad46","#0312ad"), color = "black", outlier.shape = 21, outlier.fill = "red") +
  geom_jitter(
    aes(color = clusternames),
    width = 0.15,
    size = 2.5,
    alpha = 0.9
  ) +
  # geom_text(
  #   data = stats_donor,
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
    title = "Cluster Proportions Across Donor Samples",
    x = "Cluster",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )







#PAH
meta_df_PAH <- ECstory.subset.fil.PAH@meta.data

# Step 1: Count cells per clusternames × samplecondition × sample.ident
cell_counts_PAH <- meta_df_PAH %>%
  count(clusternames, samplecondition, sample.ident, name = "cell_count")

# Step 2: Calculate percentage per sample.ident
cell_counts_PAH <- cell_counts_PAH %>%
  group_by(sample.ident) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

# View results
print(cell_counts_PAH, n = 35)


summary_stats_PAH <- cell_counts_PAH %>%
  group_by(clusternames) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )

print(summary_stats_PAH)


ggplot(summary_stats_PAH, aes(x="", y=median_percent, fill=clusternames)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=c("#b20505","#0312ad","#03ad46"))



ggplot(cell_counts_PAH, aes(x = clusternames, y = percent)) +
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


stats_PAH <- cell_counts_PAH %>%
  group_by(clusternames) %>%
  summarise(
    median = median(percent),
    q1 = quantile(percent, 0.25),
    q3 = quantile(percent, 0.75),
    ymax = max(percent),
    .groups = "drop"
  ) %>%
  mutate(label_y =50 )   # adaptive offset

ggplot(cell_counts_PAH, aes(x = clusternames, y = percent)) +
  geom_boxplot(fill = c("#b20505","#03ad46","#0312ad"), color = "black", outlier.shape = 21, outlier.fill = "red") +
  geom_jitter(
    aes(color = clusternames),
    width = 0.15,
    size = 2.5,
    alpha = 0.9
  ) +
  # geom_text(
  #   data = stats_PAH,
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
    title = "Cluster Proportions Across PAH Samples",
    x = "Cluster",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  expand_limits(y = max(stats_PAH$label_y) + 2)

#PHPF
meta_df_PHPF <- ECstory.subset.fil.PF@meta.data

# Step 1: Count cells per clusternames × samplecondition × sample.ident
cell_counts_PHPF <- meta_df_PHPF %>%
  count(clusternames, samplecondition, sample.ident, name = "cell_count")

# Step 2: Calculate percentage per sample.ident
cell_counts_PHPF <- cell_counts_PHPF %>%
  group_by(sample.ident) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

# View results
print(cell_counts_PHPF, n = 35)


summary_stats_PHPF <- cell_counts_PHPF %>%
  group_by(clusternames) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )

print(summary_stats_PHPF)


ggplot(summary_stats_PHPF, aes(x="", y=median_percent, fill=clusternames)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=c("#b20505","#0312ad","#03ad46"))



ggplot(cell_counts_PHPF, aes(x = clusternames, y = percent)) +
  geom_boxplot(fill = c("#b20505","#0312ad","#03ad46"), color = "black", outlier.shape = 21, outlier.fill = "red") +
  scale_x_discrete(limits = rev) +
  labs(
    title = "Distribution of Cluster Percentages Across Samples",
    x = "Cluster Name",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

stats_PHPF <- cell_counts_PHPF %>%
  group_by(clusternames) %>%
  summarise(
    median = median(percent),
    q1 = quantile(percent, 0.25),
    q3 = quantile(percent, 0.75),
    ymax = max(percent),
    .groups = "drop"
  ) %>%
  mutate(label_y = 50)   # adaptive offset

ggplot(cell_counts_PHPF, aes(x = clusternames, y = percent)) +
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




summary_Donor <- cell_counts_Donor %>%
  group_by(clusternames) %>%
  summarise(
    mean = mean(percent),
    median = median(percent),
    sd = sd(percent),
    IQR = IQR(percent),
    n = n()
  )

summary_PAH <- cell_counts_PAH %>%
  group_by(clusternames) %>%
  summarise(
    mean = mean(percent),
    median = median(percent),
    sd = sd(percent),
    IQR = IQR(percent),
    n = n()
  )

summary_PHPF <- cell_counts_PHPF %>%
  group_by(clusternames) %>%
  summarise(
    mean = mean(percent),
    median = median(percent),
    sd = sd(percent),
    IQR = IQR(percent),
    n = n()
  )






write_xlsx(list(
  Donor_raw = cell_counts_Donor,
  PAH_raw = cell_counts_PAH,
  PHPF_raw = cell_counts_PHPF,
  Donor_summary = summary_Donor,
  PAH_summary = summary_PAH,
  PHPF_summary = summary_PHPF
), "EC_cluster_proportions_ALL_VALUES.xlsx")



Donor_raw <- cell_counts_Donor %>% mutate(condition = "Donor")
PAH_raw <- cell_counts_PAH %>% mutate(condition = "PAH")
PHPF_raw <- cell_counts_PHPF %>% mutate(condition = "PHPF")

combined_raw <- bind_rows(Donor_raw, PAH_raw, PHPF_raw)

write_xlsx(combined_raw, "EC_cluster_proportions_ALL_RAW.xlsx")


combined_summary <- combined_raw %>%
  group_by(condition, clusternames) %>%
  summarise(
    mean = mean(percent),
    median = median(percent),
    sd = sd(percent),
    IQR = IQR(percent),
    .groups = "drop"
  )

write_xlsx(list(
  Raw_values = combined_raw,
  Summary_stats = combined_summary
), "EC_cluster_proportions_combined.xlsx")


combined_full <- left_join(combined_raw, combined_summary,
                           by = c("condition","clusternames"))

write_xlsx(combined_full, "EC_cluster_proportions_SINGLE_TABLE.xlsx")



#Phase
#Donor
meta_df_Donor <- ECstory.subset.fil.Donor@meta.data

# Step 1: Count cells per clusternames × samplecondition × sample.ident
cell_counts_Donor <- meta_df_Donor %>%
  count(Phase, samplecondition, sample.ident, name = "cell_count")

# Step 2: Calculate percentage per sample.ident
cell_counts_Donor <- cell_counts_Donor %>%
  group_by(sample.ident) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

# View results
print(cell_counts_Donor, n = 35)


summary_stats_Donor <- cell_counts_Donor %>%
  group_by(Phase) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )

print(summary_stats_Donor)


ggplot(summary_stats_Donor, aes(x="", y=median_percent, fill=Phase)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=myphasecol)



ggplot(cell_counts_Donor, aes(x = Phase, y = percent)) +
  geom_boxplot(fill = myphasecol, color = "black", outlier.shape = 21, outlier.fill = "red") +
  labs(
    title = "Distribution of Cluster Percentages Across Samples",
    x = "Cluster Name",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


stats_Donor <- cell_counts_Donor %>%
  group_by(Phase) %>%
  summarise(
    median = median(percent),
    q1 = quantile(percent, 0.25),
    q3 = quantile(percent, 0.75),
    ymax = max(percent),
    .groups = "drop"
  ) %>%
  mutate(label_y = 50)   # adaptive offset

ggplot(cell_counts_Donor, aes(x = Phase, y = percent)) +
  geom_boxplot(
    fill = c("#804daa", "#c7c881", "#525250"),
    color = "black",
    outlier.shape = 21,
    outlier.fill = "red"
  ) +
  geom_jitter(
    aes(color = Phase),
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
  scale_fill_manual(values = c("#C99DED","#E0E0BF","#7D7D7D")) +
  scale_color_manual(values = c("#C99DED","#E0E0BF","#7D7D7D")) +
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





#PAH
meta_df_PAH <- ECstory.subset.fil.PAH@meta.data

# Step 1: Count cells per clusternames × samplecondition × sample.ident
cell_counts_PAH <- meta_df_PAH %>%
  count(Phase, samplecondition, sample.ident, name = "cell_count")

# Step 2: Calculate percentage per sample.ident
cell_counts_PAH <- cell_counts_PAH %>%
  group_by(sample.ident) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

# View results
print(cell_counts_PAH, n = 35)


summary_stats_PAH <- cell_counts_PAH %>%
  group_by(Phase) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )

print(summary_stats_PAH)


ggplot(summary_stats_PAH, aes(x="", y=median_percent, fill=Phase)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=myphasecol)



ggplot(cell_counts_PAH, aes(x = Phase, y = percent)) +
  geom_boxplot(fill = myphasecol, color = "black", outlier.shape = 21, outlier.fill = "red") +
  labs(
    title = "Distribution of Cluster Percentages Across Samples",
    x = "Cluster Name",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

stats_PAH <- cell_counts_PAH %>%
  group_by(Phase) %>%
  summarise(
    median = median(percent),
    q1 = quantile(percent, 0.25),
    q3 = quantile(percent, 0.75),
    ymax = max(percent),
    .groups = "drop"
  ) %>%
  mutate(label_y = 50)   # adaptive offset

ggplot(cell_counts_PAH, aes(x = Phase, y = percent)) +
  geom_boxplot(
    fill = c("#804daa", "#c7c881", "#525250"),
    color = "black",
    outlier.shape = 21,
    outlier.fill = "red"
  ) +
  geom_jitter(
    aes(color = Phase),
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
  scale_fill_manual(values = c("#C99DED","#E0E0BF","#7D7D7D")) +
  scale_color_manual(values = c("#C99DED","#E0E0BF","#7D7D7D")) +
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


#PHPF
meta_df_PHPF <- ECstory.subset.fil.PF@meta.data

# Step 1: Count cells per clusternames × samplecondition × sample.ident
cell_counts_PHPF <- meta_df_PHPF %>%
  count(Phase, samplecondition, sample.ident, name = "cell_count")

# Step 2: Calculate percentage per sample.ident
cell_counts_PHPF <- cell_counts_PHPF %>%
  group_by(sample.ident) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  ungroup()

# View results
print(cell_counts_PHPF, n = 35)


summary_stats_PHPF <- cell_counts_PHPF %>%
  group_by(Phase) %>%
  summarise(
    mean_percent = mean(percent),
    median_percent = median(percent),
    sd_percent = sd(percent),
    variance_percent = var(percent),
    IQR_percent = IQR(percent),
    .groups = "drop"
  )

print(summary_stats_PHPF)


ggplot(summary_stats_PHPF, aes(x="", y=median_percent, fill=Phase)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  scale_fill_manual(values=myphasecol)



ggplot(cell_counts_PHPF, aes(x = Phase, y = percent)) +
  geom_boxplot(fill = myphasecol, color = "black", outlier.shape = 21, outlier.fill = "red") +
  scale_x_discrete(limits = rev) +
  labs(
    title = "Distribution of Cluster Percentages Across Samples",
    x = "Cluster Name",
    y = "Percentage of Cells per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


stats_PHPF <- cell_counts_PHPF %>%
  group_by(Phase) %>%
  summarise(
    median = median(percent),
    q1 = quantile(percent, 0.25),
    q3 = quantile(percent, 0.75),
    ymax = max(percent),
    .groups = "drop"
  ) %>%
  mutate(label_y = 50)   # adaptive offset

ggplot(cell_counts_PHPF, aes(x = Phase, y = percent)) +
  geom_boxplot(
    fill = c("#804daa", "#c7c881", "#525250"),
    color = "black",
    outlier.shape = 21,
    outlier.fill = "red"
  ) +
  geom_jitter(
    aes(color = Phase),
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
  scale_fill_manual(values = c("#C99DED","#E0E0BF","#7D7D7D")) +
  scale_color_manual(values = c("#C99DED","#E0E0BF","#7D7D7D")) +
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










#bubblechart
Donorcells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.Donor$Phase == "S")) ,
          sum((ECstory.subset.fil.Donor$Phase == "G1")),
          sum((ECstory.subset.fil.Donor$Phase == "G2M")))
)

ECstory.subset.fil.Donor.c0 <- subset(ECstory.subset.fil.Donor, seurat_clusters %in% c("0"))  
ECstory.subset.fil.PAH.c0 <- subset(ECstory.subset.fil.PAH, seurat_clusters %in% c("0")) 
ECstory.subset.fil.PF.c0 <- subset(ECstory.subset.fil.PF, seurat_clusters %in% c("0"))


ECstory.subset.fil.Donor.c1 <- subset(ECstory.subset.fil.Donor, seurat_clusters %in% c("1"))  
ECstory.subset.fil.PAH.c1 <- subset(ECstory.subset.fil.PAH, seurat_clusters %in% c("1")) 
ECstory.subset.fil.PF.c1 <- subset(ECstory.subset.fil.PF, seurat_clusters %in% c("1"))

ECstory.subset.fil.Donor.c2 <- subset(ECstory.subset.fil.Donor, seurat_clusters %in% c("2"))  
ECstory.subset.fil.PAH.c2 <- subset(ECstory.subset.fil.PAH, seurat_clusters %in% c("2")) 
ECstory.subset.fil.PF.c2 <- subset(ECstory.subset.fil.PF, seurat_clusters %in% c("2"))


Donor.c0.cells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.Donor.c0$Phase == "S")) ,
          sum((ECstory.subset.fil.Donor.c0$Phase == "G1")),
          sum((ECstory.subset.fil.Donor.c0$Phase == "G2M")))
)

PAH.c0.cells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.PAH.c0$Phase == "S")) ,
          sum((ECstory.subset.fil.PAH.c0$Phase == "G1")),
          sum((ECstory.subset.fil.PAH.c0$Phase == "G2M")))
)

PF.c0.cells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.PF.c0$Phase == "S")) ,
          sum((ECstory.subset.fil.PF.c0$Phase == "G1")),
          sum((ECstory.subset.fil.PF.c0$Phase == "G2M")))
)


Donor.c0.cells <- Donor.c0.cells %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(Donor.c0.cells$value) *100)

PAH.c0.cells <- PAH.c0.cells %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(PAH.c0.cells$value) *100)

PF.c0.cells <- PF.c0.cells %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(PF.c0.cells$value) *100)

Donor.c0.cells$condition <- "Donor"
PAH.c0.cells$condition <- "PAH"
PF.c0.cells$condition <- "PF"

Donor.c0.cells$cluster <- "0"
PAH.c0.cells$cluster <- "0"
PF.c0.cells$cluster <- "0"




Donor.c1.cells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.Donor.c1$Phase == "S")) ,
          sum((ECstory.subset.fil.Donor.c1$Phase == "G1")),
          sum((ECstory.subset.fil.Donor.c1$Phase == "G2M")))
)

PAH.c1.cells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.PAH.c1$Phase == "S")) ,
          sum((ECstory.subset.fil.PAH.c1$Phase == "G1")),
          sum((ECstory.subset.fil.PAH.c1$Phase == "G2M")))
)

PF.c1.cells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.PF.c1$Phase == "S")) ,
          sum((ECstory.subset.fil.PF.c1$Phase == "G1")),
          sum((ECstory.subset.fil.PF.c1$Phase == "G2M")))
)

Donor.c1.cells <- Donor.c1.cells %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(Donor.c1.cells$value) *100)

PAH.c1.cells <- PAH.c1.cells %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(PAH.c1.cells$value) *100)

PF.c1.cells <- PF.c1.cells %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(PF.c1.cells$value) *100)

Donor.c1.cells$condition <- "Donor"
PAH.c1.cells$condition <- "PAH"
PF.c1.cells$condition <- "PF"

Donor.c1.cells$cluster <- "1"
PAH.c1.cells$cluster <- "1"
PF.c1.cells$cluster <- "1"


Donor.c2.cells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.Donor.c2$Phase == "S")) ,
          sum((ECstory.subset.fil.Donor.c2$Phase == "G1")),
          sum((ECstory.subset.fil.Donor.c2$Phase == "G2M")))
)

PAH.c2.cells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.PAH.c2$Phase == "S")) ,
          sum((ECstory.subset.fil.PAH.c2$Phase == "G1")),
          sum((ECstory.subset.fil.PAH.c2$Phase == "G2M")))
)

PF.c2.cells <- data.frame(
  group=c("S", "G1", "G2M"),
  value=c(sum((ECstory.subset.fil.PF.c2$Phase == "S")) ,
          sum((ECstory.subset.fil.PF.c2$Phase == "G1")),
          sum((ECstory.subset.fil.PF.c2$Phase == "G2M")))
)


Donor.c2.cells <- Donor.c2.cells %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(Donor.c2.cells$value) *100)

PAH.c2.cells <- PAH.c2.cells %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(PAH.c2.cells$value) *100)

PF.c2.cells <- PF.c2.cells %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(PF.c2.cells$value) *100)

Donor.c2.cells$condition <- "Donor"
PAH.c2.cells$condition <- "PAH"
PF.c2.cells$condition <- "PF"

Donor.c2.cells$cluster <- "2"
PAH.c2.cells$cluster <- "2"
PF.c2.cells$cluster <- "2"

Bubblechart <- rbind(Donor.c0.cells,PAH.c0.cells,PF.c0.cells,
                     Donor.c1.cells,PAH.c1.cells,PF.c1.cells,
                     Donor.c2.cells,PAH.c2.cells,PF.c2.cells)

Bubblechart.D <- rbind(Donor.c0.cells,Donor.c1.cells,Donor.c2.cells)
Bubblechart.D$color <- ""
Bubblechart.D$color[Bubblechart.D$cluster == 0] <- "#b20505"
Bubblechart.D$color[Bubblechart.D$cluster == 1] <- "#0312ad"
Bubblechart.D$color[Bubblechart.D$cluster == 2] <- "#03ad46"

# ggplot(Bubblechart.D, aes(x=cluster, y=group, size=prop, color=color))+
#   geom_point(aes(alpha=0.9))+
#   geom_point(shape = 1,colour = "black")+
#   scale_size(range=c(1,10),breaks=c(0,10,20,30,40,50,60,70,80,90,100),labels=c("0%","10%","20%","30%","40%","50%", "60%","70%","80%","90%","100%"),guide="legend")+
#   theme_light()

ggplot(Bubblechart.D, aes(x=cluster, y=group, size=prop, color=color))+
  geom_point(aes(alpha=0.9))+
  geom_point(shape = 1,colour = "black")+
  scale_size(range=c(1,10),breaks=c(20,40,60),labels=c("20%","40%","60%"),guide="legend", limits=c(20, 60))+
  theme_light()






Bubblechart.PAH <- rbind(PAH.c0.cells,PAH.c1.cells,PAH.c2.cells)
Bubblechart.PAH$color <- ""
Bubblechart.PAH$color[Bubblechart.PAH$cluster == 0] <- "#b20505"
Bubblechart.PAH$color[Bubblechart.PAH$cluster == 1] <- "#0312ad"
Bubblechart.PAH$color[Bubblechart.PAH$cluster == 2] <- "#03ad46"

ggplot(Bubblechart.PAH, aes(x=cluster, y=group, size=prop, color=color))+
  geom_point(aes(alpha=0.9))+
  geom_point(shape = 1,colour = "black")+
  scale_size(range=c(1,10),breaks=c(20,40,60),labels=c("20%","40%","60%"),guide="legend", limits=c(20, 60))+
  theme_light()



Bubblechart.PF <- rbind(PF.c0.cells,PF.c1.cells,PF.c2.cells)
Bubblechart.PF$color <- ""
Bubblechart.PF$color[Bubblechart.PF$cluster == 0] <- "#b20505"
Bubblechart.PF$color[Bubblechart.PF$cluster == 1] <- "#0312ad"
Bubblechart.PF$color[Bubblechart.PF$cluster == 2] <- "#03ad46"

ggplot(Bubblechart.PF, aes(x=cluster, y=group, size=prop, color=color))+
  geom_point(aes(alpha=0.9))+
  geom_point(shape = 1,colour = "black")+
  scale_size(range=c(1,10),breaks=c(20,40,60),labels=c("20%","40%","60%"),guide="legend", limits=c(20, 60))+
  theme_light()



ggplot(Bubblechart.D,
       aes(x = cluster, y = group, size = prop, color = factor(cluster))) +
  geom_point(alpha = 0.9) +
  geom_point(shape = 1, colour = "black") +
  scale_color_manual(
    values = c(
      "0" = "#b20505",
      "1" = "#0312ad",
      "2" = "#03ad46"
    )
  ) +
  scale_size(
    range = c(1,10),
    breaks = c(20,40,60),
    labels = c("20%","40%","60%"),
    limits = c(20, 60)
  ) +
  theme_light()


ggplot(Bubblechart.PAH,
       aes(x = cluster, y = group, size = prop, color = factor(cluster))) +
  geom_point(alpha = 0.9) +
  geom_point(shape = 1, colour = "black") +
  scale_color_manual(
    values = c(
      "0" = "#b20505",
      "1" = "#0312ad",
      "2" = "#03ad46"
    )
  ) +
  scale_size(
    range = c(1,10),
    breaks = c(20,40,60),
    labels = c("20%","40%","60%"),
    limits = c(15, 60)
  ) +
  theme_light()


ggplot(Bubblechart.PF,
       aes(x = cluster, y = group, size = prop, color = factor(cluster))) +
  geom_point(alpha = 0.9) +
  geom_point(shape = 1, colour = "black") +
  scale_color_manual(
    values = c(
      "0" = "#b20505",
      "1" = "#0312ad",
      "2" = "#03ad46"
    )
  ) +
  scale_size(
    range = c(1,10),
    breaks = c(20,40,60),
    labels = c("20%","40%","60%"),
    limits = c(20, 60)
  ) +
  theme_light()


#03022026
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










plot_cluster_boxpoints <- function(df, title = NULL) {
  
  comparisons <- list(
    c("Immunomodulatory", "Vascular plasticity modulatory"),
    c("Immunomodulatory", "Vascular tone-modulatory"),
    c("Vascular plasticity modulatory", "Vascular tone-modulatory")
  )
  
  ggplot(df, aes(x = clusternames, y = percent)) +
    
    geom_boxplot(
      aes(fill = clusternames),
      width = 0.6,
      alpha = 0.9,
      color = "black",
      outlier.shape = NA
    ) +
    
    geom_jitter(
      aes(color = clusternames),
      width = 0.15,
      size = 1.5,
      alpha = 0.95
    ) +
    
    # global 3-group test
    stat_compare_means(
      method = "kruskal.test",
      label.y = 110
    ) +
    
    # pairwise post-hoc tests (adjusted)
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      method.args = list(exact = FALSE),
      label = "p.format",
      tip.length = 0.01,
      size = 4,
      label.y = c(85, 92, 99)
    ) +
    
    scale_fill_manual(values = c(
      "Immunomodulatory" = "#b20505",
      "Vascular plasticity modulatory" = "#03ad46",
      "Vascular tone-modulatory" = "#0312ad"
    )) +
    
    scale_color_manual(values = c(
      "Immunomodulatory" = "#CC4B4B",
      "Vascular plasticity modulatory" = "#57D487",
      "Vascular tone-modulatory" = "#454FD1"
    )) +
    
    coord_cartesian(ylim = c(0, 150)) +
    
    labs(
      title = title,
      x = "Cluster",
      y = "Percentage of Cells per Sample"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


plot_cluster_boxpoints(
  cell_counts_Donor,
  title = "Donor – Endothelial compartment composition"
)


plot_cluster_boxpoints(
  cell_counts_PHPF,
  title = "PHPF – Endothelial compartment composition"
)











plot_cluster_boxpoints <- function(df, title = NULL) {
  
  comparisons <- list(
    c("Immunomodulatory", "Vascular plasticity modulatory"),
    c("Immunomodulatory", "Vascular tone-modulatory"),
    c("Vascular plasticity modulatory", "Vascular tone-modulatory")
  )
  
  ggplot(df, aes(x = clusternames, y = percent)) +
    
    geom_boxplot(
      aes(fill = clusternames),
      width = 0.6,
      alpha = 0.9,
      color = "black",
      outlier.shape = NA
    ) +
    
    geom_jitter(
      aes(color = clusternames),
      width = 0.15,
      size = 1.5,
      alpha = 0.95
    ) +
    
    
    scale_fill_manual(values = c(
      "Immunomodulatory" = "#b20505",
      "Vascular plasticity modulatory" = "#03ad46",
      "Vascular tone-modulatory" = "#0312ad"
    )) +
    
    scale_color_manual(values = c(
      "Immunomodulatory" = "#CC4B4B",
      "Vascular plasticity modulatory" = "#57D487",
      "Vascular tone-modulatory" = "#454FD1"
    )) +
    
    coord_cartesian(ylim = c(0, 100)) +
    
    labs(
      title = title,
      x = "Cluster",
      y = "Percentage of Cells per Sample"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#PHASE plots

plot_cluster_boxpoints <- function(df, title = NULL) {
  
  comparisons <- list(
    c("G1", "G2M"),
    c("G1", "S"),
    c("S", "G2M")
  )
  
  ggplot(df, aes(x = Phase, y = percent)) +
    
    geom_boxplot(
      aes(fill = Phase),
      width = 0.6,
      alpha = 0.9,
      color = "black",
      outlier.shape = NA
    ) +
    
    geom_jitter(
      aes(color = Phase),
      width = 0.15,
      size = 1.5,
      alpha = 0.95
    ) +
    
    # global 3-group test
    stat_compare_means(
      method = "kruskal.test",
      label.y = 110
    ) +
    
    # pairwise post-hoc tests (adjusted)
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      method.args = list(exact = FALSE),
      label = "p.format",
      tip.length = 0.01,
      size = 4,
      label.y = c(85, 92, 99)
    ) +
    
    scale_fill_manual(values = c(
      "G1" = "#804daa",
      "G2M" = "#c7c881",
      "S" = "#323227"
    )) +
   
    scale_color_manual(values = c(
      "G1" = "#C99DED",
      "G2M" = "#E0E0BF",
      "S" = "#7D7D7D"
    )) +
    
    coord_cartesian(ylim = c(0, 150)) +
    
    labs(
      title = title,
      x = "Cluster",
      y = "Percentage of Cells per Sample"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


plot_cluster_boxpoints(
  cell_counts_PHPF,
  title = "Donor – Endothelial compartment composition"
)














plot_cluster_boxpoints <- function(df, title = NULL) {
  
  comparisons <- list(
    c("G1", "G2M"),
    c("G1", "S"),
    c("S", "G2M")
  )
  
  ggplot(df, aes(x = Phase, y = percent)) +
    
    geom_boxplot(
      aes(fill = Phase),
      width = 0.6,
      alpha = 0.9,
      color = "black",
      outlier.shape = NA
    ) +
    
    geom_jitter(
      aes(color = Phase),
      width = 0.15,
      size = 1.5,
      alpha = 0.95
    ) +
    
    
    scale_fill_manual(values = c(
      "G1" = "#804daa",
      "G2M" = "#c7c881",
      "S" = "#323227"
    )) +
    
    scale_color_manual(values = c(
      "G1" = "#C99DED",
      "G2M" = "#E0E0BF",
      "S" = "#7D7D7D"
    )) +
    
    coord_cartesian(ylim = c(0, 100)) +
    
    labs(
      title = title,
      x = "Cluster",
      y = "Percentage of Cells per Sample"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

### box


meta_df <- ECstory.subset.fil@meta.data

cell_counts_phase <- meta_df %>%
  count(sample.ident, samplecondition, Phase, name = "cell_count") %>%
  group_by(sample.ident) %>%
  mutate(percent = cell_count / sum(cell_count) * 100) %>%
  ungroup()

library(ggpubr)

ggplot(cell_counts_phase,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  
  geom_jitter(aes(color = Phase),
              width = 0.15,
              size = 2) +
  
  stat_compare_means(
    aes(group = Phase),
    method = "kruskal.test"
  ) +
  
  facet_wrap(~Phase) +
  
  scale_fill_manual(values = c(
    "G1" = "#804daa",
    "S" = "#323227",
    "G2M" = "#c7c881"
  )) +
  
  scale_color_manual(values = c(
    "G1" = "#C99DED",
    "S" = "#7D7D7D",
    "G2M" = "#E0E0BF"
  )) +
  
  labs(
    x = "Condition",
    y = "Percentage of cells",
    title = "Cell cycle distribution across conditions"
  ) +
  
  theme_classic(base_size = 14)






comparisons <- list(
  c("Donor","IPAH"),
  c("Donor","IPF"),
  c("IPAH","IPF")
)

ggplot(cell_counts_phase,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test") +
  
  facet_wrap(~Phase) +
  
  theme_classic()




comparisons <- list(
  c("Donor","IPAH"),
  c("Donor","IPF"),
  c("IPAH","IPF")
)

ggplot(cell_counts_phase,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Phase), width = 0.15, size = 2) +
  
  facet_wrap(~Phase) +
  
  # GLOBAL test per facet
  stat_compare_means(method = "kruskal.test", label.y = 70) +
  
  # Pairwise tests
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format"
  ) +
  
  scale_fill_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  scale_color_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  theme_classic()





stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",
  label = "p.signif",
  hide.ns = TRUE
)


ggplot(cell_counts_phase,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Phase), width = 0.15, size = 2) +
  
  facet_wrap(~Phase) +
  
  stat_compare_means(method = "kruskal.test") +
  
  stat_compare_means(
    aes(group = samplecondition),
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  
  scale_fill_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  scale_color_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  theme_classic(base_size = 14)



comparisons <- list(
  c("Donor","IPAH"),
  c("Donor","IPF"),
  c("IPAH","IPF")
)

ggplot(cell_counts_phase,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Phase), width = 0.15, size = 2) +
  
  facet_wrap(~Phase) +
  
  # GLOBAL test per phase
  stat_compare_means(
    aes(x = samplecondition),
    method = "kruskal.test"
  ) +
  
  # Pairwise tests within each phase
  stat_compare_means(
    aes(x = samplecondition),
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    step.increase = 0.08
  ) +
  
  scale_fill_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  scale_color_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  theme_classic(base_size = 14)










library(rstatix)

# Global Kruskal per Phase
kruskal_res <- cell_counts_phase %>%
  group_by(Phase) %>%
  kruskal_test(percent ~ samplecondition)

# Pairwise Wilcoxon per Phase
pairwise_res <- cell_counts_phase %>%
  group_by(Phase) %>%
  pairwise_wilcox_test(percent ~ samplecondition, p.adjust.method = "BH") %>%
  add_y_position()

ggplot(cell_counts_phase,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Phase), width = 0.15, size = 2) +
  
  facet_wrap(~Phase) +
  
  stat_pvalue_manual(pairwise_res, label = "p.adj.signif", hide.ns = TRUE) +
  
  scale_fill_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  scale_color_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  theme_classic(base_size = 14)




kruskal_res <- kruskal_res %>% add_y_position()

stat_pvalue_manual(kruskal_res, label = "p")


ggplot(cell_counts_phase,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Phase), width = 0.15, size = 2) +
  
  facet_wrap(~Phase) +
  
  stat_pvalue_manual(kruskal_res, label = "p") +
  stat_pvalue_manual(pairwise_res, label = "p.adj") +
  
  scale_fill_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  scale_color_manual(values = c(
    "G1" = "#F8766D",
    "G2M" = "#00BA38",
    "S" = "#619CFF"
  )) +
  
  theme_classic(base_size = 14)


pairwise_res <- cell_counts_phase %>%
  group_by(Phase) %>%
  pairwise_wilcox_test(percent ~ samplecondition, p.adjust.method = "BH") %>%
  add_xy_position(x = "samplecondition")


ggplot(cell_counts_phase,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Phase), width = 0.15, size = 2) +
  
  facet_wrap(~Phase) +
    stat_pvalue_manual(
    pairwise_res,
    label = "p.adj",
    hide.ns = FALSE
  ) +

  
  scale_fill_manual(values = c(
    "G1" = "#804daa",
    "S" = "#7D7D7D",
    "G2M" = "#c7c881"
  )) +
  
  scale_color_manual(values = c(
    "G1" = "#C99DED",
    "G2M" = "#E0E0BF",
    "S" = "#ABA9A9"
  )) +
  
  theme_classic(base_size = 14)

kruskal_res %>%
  mutate(
    p = signif(p, 3),
    statistic = round(statistic, 2)
  )


stack_df <- cell_counts_phase %>%
  group_by(samplecondition, Phase) %>%
  summarise(mean_percent = mean(percent), .groups = "drop")

ggplot(stack_df,
       aes(x = samplecondition, y = mean_percent, fill = Phase)) +
  
  geom_bar(stat = "identity", color = "black") +
  
  scale_fill_manual(values = c(
    "G1" = "#804daa",
    "S" = "#7D7D7D",
    "G2M" = "#c7c881"
  )) +
  
  labs(
    x = "Condition",
    y = "Mean % of cells",
    title = "Cell cycle composition across conditions"
  ) +
  
  theme_classic(base_size = 14)


#alt bubble

phase_cluster_df <- ECstory.subset.fil@meta.data %>%
  count(samplecondition, clusternames, Phase) %>%
  group_by(samplecondition, clusternames) %>%
  mutate(percent = n / sum(n) * 100)

ggplot(phase_cluster_df,
       aes(x = Phase, y = clusternames, fill = percent)) +
  
  geom_tile(color = "white") +
  
  facet_wrap(~samplecondition) +
  
  scale_fill_gradient(low = "white", high = "#804daa") +
  
  labs(
    x = "Cell cycle phase",
    y = "EC subpopulation",
    fill = "% cells"
  ) +
  
  theme_classic(base_size = 14)


ggplot(phase_cluster_df,
       aes(x = clusternames, y = percent, fill = Phase)) +
  
  geom_bar(stat = "identity", color ="black") +
  
  facet_wrap(~samplecondition) +
  
  scale_fill_manual(values = c(
    "G1" = "#804daa",
    "S" = "#7D7D7D",
    "G2M" = "#c7c881"
  )) +
  
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))











phase_cluster_sample <- ECstory.subset.fil@meta.data %>%
  count(sample.ident, samplecondition, clusternames, Phase) %>%
  group_by(sample.ident, clusternames) %>%
  mutate(percent = n / sum(n) * 100)

ggplot(phase_cluster_sample,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15) +
  
  facet_grid(Phase ~ clusternames) +
  
  theme_classic(base_size = 14)



ggplot(phase_cluster_sample,
       aes(x = samplecondition, y = percent, fill = Phase)) +
  
  geom_boxplot(outlier.shape = NA) +
  
  geom_jitter(aes(color = Phase), width = 0.15,size = 1, color="black", shape=21) +
  
  facet_grid(Phase ~ clusternames) +
  
  scale_fill_manual(values = c(
    "G1" = "#804daa",
    "G2M" = "#c7c881",
    "S" = "#A3A3A3"
  )) +
  
  scale_color_manual(values = c(
    "G1" = "#C99DED",
    "G2M" = "#E0E0BF",
    "S" = "#DEDCDC"
  )) +
  
  theme_classic(base_size = 14)



library(dplyr)

phase_cluster_sample %>%
  group_by(samplecondition, clusternames, Phase) %>%
  summarise(
    median = median(percent),
    mean = mean(percent),
    sd = sd(percent),
    IQR = IQR(percent),
    n = n()
  )


stats_table <- phase_cluster_sample %>%
  group_by(samplecondition, clusternames, Phase) %>%
  summarise(
    median = median(percent),
    IQR = IQR(percent)
  )

write.csv(stats_table, "CellCycle_cluster_statistics.csv")




library(dplyr)
library(rstatix)
library(tidyr)
library(purrr)
library(writexl)

# Use replicate-level dataset
df <- phase_cluster_sample

df$samplecondition <- as.factor(df$samplecondition)
df$clusternames <- as.factor(df$clusternames)
df$Phase <- as.factor(df$Phase)

# ============================================
# Function: stats per Phase × Cluster
# ============================================

analyze_group <- function(data) {
  
  data$samplecondition <- droplevels(data$samplecondition)
  
  phase_val <- unique(data$Phase)
  cluster_val <- unique(data$clusternames)
  
  # -------------------------
  # Descriptive statistics
  # -------------------------
  desc <- data %>%
    group_by(samplecondition) %>%
    summarise(
      n = n(),
      Median = median(percent),
      Q1 = quantile(percent, 0.25),
      Q3 = quantile(percent, 0.75),
      IQR = IQR(percent),
      .groups = "drop"
    )
  
  # Join descriptive stats back to RAW rows
  data_joined <- data %>%
    left_join(desc, by = "samplecondition")
  
  # If fewer than 2 groups → skip inferential stats
  if (nlevels(data$samplecondition) < 2) {
    
    data_joined <- data_joined %>%
      mutate(
        KW_H = NA,
        KW_df = NA,
        KW_p = NA,
        KW_epsilon2 = NA
      )
    
    return(data_joined)
  }
  
  # -------------------------
  # Kruskal–Wallis
  # -------------------------
  kw <- kruskal_test(data, percent ~ samplecondition)
  
  H <- kw$statistic
  df_kw <- kw$df
  p_kw <- kw$p
  
  k <- nlevels(data$samplecondition)
  n_total <- nrow(data)
  epsilon_sq <- (H - k + 1) / (n_total - k)
  
  # -------------------------
  # Pairwise Wilcoxon
  # -------------------------
  pw <- pairwise_wilcox_test(
    data,
    percent ~ samplecondition,
    p.adjust.method = "holm"
  )
  
  eff <- wilcox_effsize(
    data,
    percent ~ samplecondition,
    comparisons = combn(levels(data$samplecondition), 2, simplify = FALSE)
  )
  
  pw_full <- pw %>%
    left_join(eff, by = c("group1", "group2")) %>%
    mutate(
      comparison = paste(group1, "vs", group2)
    ) %>%
    select(comparison, statistic, p, p.adj, effsize) %>%
    pivot_wider(
      names_from = comparison,
      values_from = c(statistic, p, p.adj, effsize),
      names_sep = "_"
    )
  
  # Add inferential stats to every replicate row
  data_final <- data_joined %>%
    mutate(
      Phase = phase_val,
      Cluster = cluster_val,
      KW_H = H,
      KW_df = df_kw,
      KW_p = p_kw,
      KW_epsilon2 = epsilon_sq
    ) %>%
    bind_cols(pw_full)
  
  return(data_final)
}

# ============================================
# Run across Phase × Cluster
# ============================================

final_table <- df %>%
  group_by(Phase, clusternames) %>%
  group_split() %>%
  map_df(analyze_group)

# ============================================
# Export ONE single sheet
# ============================================

write_xlsx(
  final_table,
  "CellCycle_FINAL_Nature_Compliant_Table.xlsx"
)




##--#-#-#-##-#-#-##--#-#-#-#-#--#-



library(dplyr)
library(rstatix)
library(tidyr)
library(purrr)
library(writexl)

# Use your replicate-level dataset
df <- PA_cluster_sample

df$samplecondition <- as.factor(df$samplecondition)
df$clusternames <- as.factor(df$clusternames)

# ============================================
# Function: stats per condition (raw + stats)
# ============================================

analyze_condition <- function(data) {
  
  data$clusternames <- droplevels(data$clusternames)
  cond_val <- unique(data$samplecondition)
  
  # -------------------------
  # Descriptive stats
  # -------------------------
  desc <- data %>%
    group_by(clusternames) %>%
    summarise(
      n = n(),
      Median = median(percent),
      Q1 = quantile(percent, 0.25),
      Q3 = quantile(percent, 0.75),
      IQR = IQR(percent),
      .groups = "drop"
    )
  
  # Join descriptive stats back to RAW rows
  data_joined <- data %>%
    left_join(desc, by = "clusternames")
  
  # If only 1 cluster → skip inferential stats
  if (nlevels(data$clusternames) < 2) {
    
    data_joined <- data_joined %>%
      mutate(
        KW_H = NA,
        KW_df = NA,
        KW_p = NA,
        KW_epsilon2 = NA
      )
    
    return(data_joined)
  }
  
  # -------------------------
  # Kruskal–Wallis
  # -------------------------
  kw <- kruskal_test(data, percent ~ clusternames)
  
  H <- kw$statistic
  df_kw <- kw$df
  p_kw <- kw$p
  
  k <- nlevels(data$clusternames)
  n_total <- nrow(data)
  epsilon_sq <- (H - k + 1) / (n_total - k)
  
  # -------------------------
  # Pairwise Wilcoxon
  # -------------------------
  pw <- pairwise_wilcox_test(
    data,
    percent ~ clusternames,
    p.adjust.method = "holm"
  )
  
  eff <- wilcox_effsize(
    data,
    percent ~ clusternames,
    comparisons = combn(levels(data$clusternames), 2, simplify = FALSE)
  )
  
  pw_full <- pw %>%
    left_join(eff, by = c("group1", "group2")) %>%
    mutate(
      comparison = paste(group1, "vs", group2)
    ) %>%
    select(comparison, statistic, p, p.adj, effsize) %>%
    pivot_wider(
      names_from = comparison,
      values_from = c(statistic, p, p.adj, effsize),
      names_sep = "_"
    )
  
  # Add inferential stats to every replicate row
  data_final <- data_joined %>%
    mutate(
      KW_H = H,
      KW_df = df_kw,
      KW_p = p_kw,
      KW_epsilon2 = epsilon_sq
    ) %>%
    bind_cols(pw_full)
  
  return(data_final)
}

# ============================================
# Run per condition
# ============================================

final_table <- df %>%
  group_by(samplecondition) %>%
  group_split() %>%
  map_df(analyze_condition)

# ============================================
# Export ONE single sheet
# ============================================

write_xlsx(
  final_table,
  "PA_cluster_distribution_FINAL_Nature_Compliant_Table.xlsx"
)



#graph fix

ECstory.subset.fil$clusternames <- recode(
  ECstory.subset.fil$seurat_clusters,
  "0" = "Vascular tone-modulatory",
  "1" = "Immunomodulatory",
  "2" = "Vascular plasticity modulatory"
)

ECstory.subset.fil$clusternames <- factor(
  ECstory.subset.fil$clusternames,
  levels = c(
    "Immunomodulatory",
    "Vascular plasticity modulatory",
    "Vascular tone-modulatory"
  )
)


Donor <- subset(ECstory.subset.fil, samplecondition == "Donor")
PAH   <- subset(ECstory.subset.fil, samplecondition == "IPAH")
PF    <- subset(ECstory.subset.fil, samplecondition == "IPF")


library(dplyr)

calc_cluster_percent <- function(obj){
  
  meta <- obj@meta.data
  
  meta %>%
    count(sample.ident, clusternames) %>%
    group_by(sample.ident) %>%
    mutate(percent = n/sum(n)*100) %>%
    ungroup()
  
}

cell_counts_Donor <- calc_cluster_percent(Donor)
cell_counts_PAH   <- calc_cluster_percent(PAH)
cell_counts_PF    <- calc_cluster_percent(PF)


plot_cluster_box <- function(df, title){
  
  ggplot(df, aes(x=clusternames, y=percent, fill=clusternames)) +
    
    geom_boxplot(width=0.6, outlier.shape=NA) +
    
    geom_jitter(
      aes(color=clusternames),
      width=0.15,
      size=2
    ) +
    
    scale_fill_manual(values=c(
      "#b20505",
      "#03ad46",
      "#0312ad"
    )) +
    
    scale_color_manual(values=c(
      "#CC4B4B",
      "#57D487",
      "#454FD1"
    )) +
    
    coord_cartesian(ylim=c(0,100)) +
    
    labs(
      title=title,
      x="Cluster",
      y="Percentage of cells"
    ) +
    
    theme_classic(base_size=14) +
    theme(
      legend.position="none",
      axis.text.x = element_text(angle=45,hjust=1)
    )
  
}


plot_cluster_box(cell_counts_Donor,
                 "Donor – EC composition")

plot_cluster_box(cell_counts_PAH,
                 "PAH – EC composition")

plot_cluster_box(cell_counts_PF,
                 "PF – EC composition")






#-#-#-#-#-#-#-#-

library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------------------------------
# Fix cluster naming
# ---------------------------------------
DimPlot(ECstory.subset.fil, split.by = "sample.ident")
ECstory.subset.fil$clusternames <- recode(
  ECstory.subset.fil$seurat_clusters,
  "0" = "Immunomodulatory",
  "1" = "Vascular tone-modulatory",
  "2" = "Vascular plasticity modulatory"
)


# ---------------------------------------
# Split conditions
# ---------------------------------------

Donor <- subset(ECstory.subset.fil, samplecondition == "Donor")
PAH   <- subset(ECstory.subset.fil, samplecondition == "IPAH")
PF    <- subset(ECstory.subset.fil, samplecondition == "IPF")

# ---------------------------------------
# Function to calculate cluster %
# ---------------------------------------

calc_cluster_percent <- function(obj){
  
  obj@meta.data %>%
    
    count(sample.ident, clusternames) %>%
    
    # ensure every cluster exists for every sample
    complete(
      sample.ident,
      clusternames,
      fill = list(n = 0)
    ) %>%
    
    group_by(sample.ident) %>%
    
    mutate(percent = n / sum(n) * 100) %>%
    
    ungroup()
}

cell_counts_Donor <- calc_cluster_percent(Donor)
cell_counts_PAH   <- calc_cluster_percent(PAH)
cell_counts_PF    <- calc_cluster_percent(PF)

# ---------------------------------------
# Plot function
# ---------------------------------------

plot_cluster_box <- function(df, title){
  
  ggplot(df, aes(x = clusternames, y = percent, fill = clusternames)) +
    
    geom_boxplot(
      width = 0.6,
      color = "black",
      outlier.shape = NA
    ) +
    
    geom_jitter(
      aes(color = clusternames),
      width = 0.15,
      size = 2,
      alpha = 0.9
    ) +
    
    scale_fill_manual(values = c(
      "Immunomodulatory" = "#b20505",
      "Vascular plasticity modulatory" = "#03ad46",
      "Vascular tone-modulatory" = "#0312ad"
    )) +
    
    scale_color_manual(values = c(
      "Immunomodulatory" = "#CC4B4B",
      "Vascular plasticity modulatory" = "#57D487",
      "Vascular tone-modulatory" = "#454FD1"
    )) +
    
    coord_cartesian(ylim = c(0,100)) +
    
    labs(
      title = title,
      x = "Cluster",
      y = "Percentage of cells"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# ---------------------------------------
# Generate plots
# ---------------------------------------

plot_cluster_box(cell_counts_Donor,
                 "Donor – EC composition")

plot_cluster_box(cell_counts_PAH,
                 "PAH – EC composition")

plot_cluster_box(cell_counts_PF,
                 "PF – EC composition")



#recreate statistics sheet

library(dplyr)
library(tidyr)
library(rstatix)
library(purrr)
library(writexl)

# ----------------------------------------------------
# 1. Extract metadata
# ----------------------------------------------------

meta_df <- ECstory.subset.fil@meta.data

# ----------------------------------------------------
# 2. Calculate cluster proportions per sample
# ----------------------------------------------------

cell_counts <- meta_df %>%
  count(sample.ident, samplecondition, clusternames, name = "n") %>%
  
  group_by(sample.ident, samplecondition) %>%
  
  tidyr::complete(
    clusternames,
    fill = list(n = 0)
  ) %>%
  
  mutate(percent = n / sum(n) * 100) %>%
  
  ungroup()

# ----------------------------------------------------
# 3. Function performing stats per condition
# ----------------------------------------------------

analyze_condition <- function(data) {
  
  data$clusternames <- droplevels(data$clusternames)
  
  # Descriptive stats
  desc <- data %>%
    group_by(clusternames) %>%
    summarise(
      Median = median(percent),
      Q1 = quantile(percent,0.25),
      Q3 = quantile(percent,0.75),
      IQR = IQR(percent),
      .groups="drop"
    )
  
  data_joined <- left_join(data, desc, by="clusternames")
  
  # Kruskal–Wallis
  kw <- kruskal_test(data, percent ~ clusternames)
  
  H <- kw$statistic
  df_kw <- kw$df
  p_kw <- kw$p
  
  k <- nlevels(data$clusternames)
  n_total <- nrow(data)
  epsilon_sq <- (H - k + 1) / (n_total - k)
  
  # Pairwise Wilcoxon
  pw <- pairwise_wilcox_test(
    data,
    percent ~ clusternames,
    p.adjust.method="holm"
  )
  
  eff <- wilcox_effsize(
    data,
    percent ~ clusternames,
    comparisons=combn(levels(data$clusternames),2,simplify=FALSE)
  )
  
  pw_full <- pw %>%
    left_join(eff, by=c("group1","group2")) %>%
    mutate(comparison=paste(group1,"vs",group2)) %>%
    select(comparison,statistic,p,p.adj,effsize) %>%
    pivot_wider(
      names_from=comparison,
      values_from=c(statistic,p,p.adj,effsize),
      names_sep="_"
    )
  
  data_final <- data_joined %>%
    mutate(
      KW_H = H,
      KW_df = df_kw,
      KW_p = p_kw,
      KW_epsilon2 = epsilon_sq
    ) %>%
    bind_cols(pw_full)
  
  return(data_final)
}

# ----------------------------------------------------
# 4. Run analysis per condition
# ----------------------------------------------------

final_table <- cell_counts %>%
  group_by(samplecondition) %>%
  group_split() %>%
  map_df(analyze_condition)

# ----------------------------------------------------
# 5. Export table
# ----------------------------------------------------

write_xlsx(
  final_table,
  "ST19_Condition_subpop_stats.xlsx"
)


#inter condition comparison

library(dplyr)
library(rstatix)
library(purrr)

# ---------------------------------------
# Ensure condition order
# ---------------------------------------

cell_counts$samplecondition <- factor(
  cell_counts$samplecondition,
  levels = c("Donor", "IPAH", "IPF")
)

# ---------------------------------------
# Function: stats per cluster
# ---------------------------------------

analyze_cluster <- function(df) {
  
  df$samplecondition <- droplevels(df$samplecondition)
  
  # Kruskal-Wallis across conditions
  kw <- kruskal_test(df, percent ~ samplecondition)
  
  H <- kw$statistic
  df_kw <- kw$df
  p_kw <- kw$p
  
  k <- nlevels(df$samplecondition)
  n_total <- nrow(df)
  epsilon_sq <- (H - k + 1) / (n_total - k)
  
  # Pairwise Wilcoxon
  pw <- pairwise_wilcox_test(
    df,
    percent ~ samplecondition,
    p.adjust.method = "holm"
  )
  
  # Effect sizes
  eff <- wilcox_effsize(
    df,
    percent ~ samplecondition,
    comparisons = combn(levels(df$samplecondition), 2, simplify = FALSE)
  )
  
  pw_full <- pw %>%
    left_join(eff, by = c("group1", "group2")) %>%
    mutate(comparison = paste(group1, "vs", group2)) %>%
    select(comparison, statistic, p, p.adj, effsize) %>%
    pivot_wider(
      names_from = comparison,
      values_from = c(statistic, p, p.adj, effsize),
      names_sep = "_"
    )
  
  # Descriptive stats
  desc <- df %>%
    group_by(samplecondition) %>%
    summarise(
      Median = median(percent),
      Q1 = quantile(percent, 0.25),
      Q3 = quantile(percent, 0.75),
      IQR = IQR(percent),
      .groups = "drop"
    )
  
  df_final <- df %>%
    left_join(desc, by = "samplecondition") %>%
    mutate(
      KW_H = H,
      KW_df = df_kw,
      KW_p = p_kw,
      KW_epsilon2 = epsilon_sq
    ) %>%
    bind_cols(pw_full)
  
  return(df_final)
}

# ---------------------------------------
# Run per cluster
# ---------------------------------------

cluster_stats <- cell_counts %>%
  group_by(clusternames) %>%
  group_split() %>%
  map_df(analyze_cluster)

# ---------------------------------------
# Optional: save
# ---------------------------------------

write.csv(cluster_stats,
          "Between_condition_cluster_stats.csv",
          row.names = FALSE)





#inter condition comparison


library(dplyr)
library(rstatix)
library(purrr)

# ---------------------------------------
# Ensure correct order
# ---------------------------------------

cell_counts$samplecondition <- factor(
  cell_counts$samplecondition,
  levels = c("Donor", "IPAH", "IPF")
)

# ---------------------------------------
# Function: stats per cluster
# ---------------------------------------

analyze_cluster <- function(df) {
  
  df$samplecondition <- droplevels(df$samplecondition)
  
  # -------------------------------
  # Global Kruskal-Wallis
  # -------------------------------
  
  kw <- kruskal_test(df, percent ~ samplecondition)
  
  H <- kw$statistic
  df_kw <- kw$df
  p_kw <- kw$p
  
  k <- nlevels(df$samplecondition)
  n_total <- nrow(df)
  epsilon_sq <- (H - k + 1) / (n_total - k)
  
  # -------------------------------
  # Pairwise: ONLY vs Donor
  # -------------------------------
  
  pw <- pairwise_wilcox_test(
    df,
    percent ~ samplecondition,
    ref.group = "Donor",   # 👈 key line
    p.adjust.method = "holm"
  )
  
  # Effect sizes ONLY vs Donor
  eff <- wilcox_effsize(
    df,
    percent ~ samplecondition,
    comparisons = list(
      c("Donor", "IPAH"),
      c("Donor", "IPF")
    )
  )
  
  pw_full <- pw %>%
    left_join(eff, by = c("group1", "group2")) %>%
    mutate(comparison = paste(group1, "vs", group2)) %>%
    select(comparison, statistic, p, p.adj, effsize) %>%
    tidyr::pivot_wider(
      names_from = comparison,
      values_from = c(statistic, p, p.adj, effsize),
      names_sep = "_"
    )
  
  # -------------------------------
  # Descriptive stats
  # -------------------------------
  
  desc <- df %>%
    group_by(samplecondition) %>%
    summarise(
      Median = median(percent),
      Q1 = quantile(percent, 0.25),
      Q3 = quantile(percent, 0.75),
      IQR = IQR(percent),
      .groups = "drop"
    )
  
  # -------------------------------
  # Final table
  # -------------------------------
  
  df_final <- df %>%
    left_join(desc, by = "samplecondition") %>%
    mutate(
      KW_H = H,
      KW_df = df_kw,
      KW_p = p_kw,
      KW_epsilon2 = epsilon_sq
    ) %>%
    bind_cols(pw_full)
  
  return(df_final)
}

# ---------------------------------------
# Run per cluster
# ---------------------------------------

cluster_stats_vs_donor <- cell_counts %>%
  group_by(clusternames) %>%
  group_split() %>%
  map_df(analyze_cluster)

# ---------------------------------------
# Save
# ---------------------------------------

write.csv(cluster_stats_vs_donor,
          "Cluster_stats_vs_Donor.csv",
          row.names = FALSE)

