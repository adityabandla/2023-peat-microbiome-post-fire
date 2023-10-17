## name: 5_btp_functional_profiles.R
## purpose: to analyse and visualise predicted and sequenced
## metagenome profiles
## author: Dr Aditya Bandla
## email: adityabandla@u.nus.edu

## load packages
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(emmeans)
  library(ggh4x)
  library(ggpattern)
  library(ggplotify)
  library(ggsignif)
  library(Hmisc)
  library(janitor)
  library(patchwork)
  library(phyloseq)
  library(scales)
  library(tidyverse)
})

## paths to directories
repo <- file.path("/Users/abandla/Desktop/2_research/1_manuscripts/2_2020_brunei_peat_fire/")
data <- file.path(repo, "1_data")
figures <- file.path(repo, "3_figures")

## set global theme options for plots
btp_theme <- theme(
  axis.text = element_text(size = 16, color = "black"),
  axis.text.y = element_text(margin = margin(0, 10, 0, 10)),
  axis.text.x = element_text(margin = margin(10, 0, 10, 0)),
  axis.title = element_text(size = 18),
  axis.ticks.length = unit(.25, "cm"),
  panel.border = element_rect(linewidth = 0.5, fill = NA),
  panel.background = element_rect(fill = NA),
  panel.grid = element_blank(),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 18),
  legend.key = element_rect(fill = NA),
  legend.background = element_rect(fill = NA)
)

## import phyloseq object
btp_fire_ps <- readRDS(file.path(data, "3_phyloseq", "2020_btp_fire_ps.rds"))
btp_fire_ps

## import nsti data
btp_nsti <- read.table(
  file.path(data, "5_picrust2", "2020_btp_fire_16S_count_nsti_pred.tsv"),
  sep = "\t",
  header = TRUE
) %>%
  filter(NSTI <= 2)

## import PICRUST2 results
## ASVs with methanogenesis modules ≥0.75 completeness
btp_ko_methanogenesis_modules <- read.table(
    file.path(data, "5_picrust2", "2020_btp_ko_methanogenesis_modules.tsv"),
    sep = "\t",
    header = TRUE
) %>%
  clean_names %>%
  filter(genome_name %in% pull(btp_nsti, ASV))

## abundance of methanogenic archaea
## 72 ASVs from class Bathyarchaeia
## 8/72 ASVs ~11% predicted to be methanogens
methanogen_abundance <- btp_fire_ps %>%
  subset_taxa(., Kingdom == "Archaea") %>%
  transform_sample_counts(., function(x) x / sum(x)) %>%
  psmelt %>%
  filter(OTU %in% pull(distinct(btp_ko_methanogenesis_modules, genome_name))) %>%
  group_by(Sample, psf_type, depth) %>%
  reframe(Abundance = round((sum(Abundance) * 100), 2))

## ASV-pathway map
asv_pathway_map <- btp_ko_methanogenesis_modules %>%
  mutate(pathway = case_when(
    str_detect(module_name, "CO2") ~ "Hydrogenotrophic",
    str_detect(module_name, "acetate") ~ "Acetoclastic",
    TRUE ~ "Methylotrophic"
  )) %>%
  mutate(id = case_when(
    pathway == "Hydrogenotrophic" ~ 1,
    pathway == "Acetoclastic" ~ 2,
    TRUE ~ 3
  )) %>%
  distinct(genome_name, pathway, id) %>%
  pivot_wider(names_from = "pathway", values_from = "id", values_fill = 0)

## abundance table for heatmap
methanogen_abundance_tab <- btp_fire_ps %>%
  subset_taxa(., Kingdom == "Archaea") %>%
  transform_sample_counts(., function(x) x / sum(x)) %>%
  psmelt %>%
  filter(OTU %in% pull(distinct(btp_ko_methanogenesis_modules, genome_name))) %>%
  group_by(psf_type, depth, OTU, Class) %>%
  reframe(Abundance = sqrt(mean(Abundance))) %>%
  arrange(psf_type, depth, Class) %>%
  group_by(psf_type, depth) %>%
  mutate(group = cur_group_id()) %>%
  ungroup %>%
  select(OTU, group, Abundance) %>% 
  mutate(OTU = parse_number(OTU)) %>%
  pivot_wider(names_from = OTU, values_from = Abundance) %>%
  column_to_rownames("group") %>%
  as.matrix

## extract class labels for heatmap
methanogen_class <- enframe(colnames(methanogen_abundance_tab), value = "ASV", name = NULL) %>%
  mutate(ASV = paste0("ASV_", ASV)) %>%
  inner_join(
    tax_table(btp_fire_ps) %>%
      data.frame() %>%
      rownames_to_column("ASV")
  ) %>%
  mutate(ASV = parse_number(ASV)) %>%
  mutate(Class = factor(Class)) %>%
  pull(Class)

## plot methanogen distribution
ht_opt$ROW_ANNO_PADDING <- unit(0.1, "cm")
options(repr.plot.width = 12, repr.plot.height = 3)
methanogen_distribution_figure <- draw(
  Heatmap(methanogen_abundance_tab,
    col = c("#ffffff", "#000000"),
    row_gap = unit(4, "mm"),
    column_gap = unit(4, "mm"),
    rect_gp = gpar(col = "black", lwd = 1.5),
    show_heatmap_legend = FALSE,
    show_row_dend = T,
    row_names_side = "right",
    column_dend_gp = gpar(lwd = 1.5),
    heatmap_legend_param = list(
      title = "Z-score",
      legend_direction = "horizontal",
      legend_width = unit(40, "mm"),
      legend_height = unit(55, "mm"),
      labels_gp = gpar(fontsize = 12),
      title_gp = gpar(fontsize = 13, fontface = "bold")
    ),
    border_gp = gpar(col = "black"),
    row_split = factor(rep(letters[1:2], each = 3), levels = letters[1:2]),
    column_split = methanogen_class,
    column_title = NULL,
    cluster_rows = FALSE,
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    clustering_distance_columns = "spearman",
    row_title = NULL,
    show_column_names = TRUE,
    top_annotation = HeatmapAnnotation(
      pathway = anno_block(
        gp = gpar(lwd = 1.5, fill = c("#ffffb3", "#8dd3c7", rep("#fb8072", 5)))),
      empty = anno_empty(border = FALSE, height = unit(1.6, "mm"))
    ),
    left_annotation = rowAnnotation(
      psf_type = anno_block(
        gp = gpar(col = 0),
        labels = c("Burnt", "Intact"),
        labels_gp = gpar(col = "Black", fontsize = 18)
      )
    )
  ) +
    rowAnnotation(
      labels = anno_text(
        c("Surface", "Mid", "Deep", "Surface", "Mid", "Deep"),
        which = "row",
        gp = gpar(fontsize = 16)
      )
    ),
  heatmap_legend_side = "bottom"
)

## export plot
pdf(file.path(figures, "1_main", "3b_methanogen_heatmap.pdf"), width = 12, height = 3)
methanogen_distribution_figure
dev.off()

## import PICRUST2 results
## ASVs with methanotrophy modules ≥0.75 completeness
btp_ko_methanotrophy_modules <- read.table(
    file.path(data, "5_picrust2", "2020_btp_ko_methanotrophy_modules.tsv"),
    sep = "\t",
    header = TRUE
) %>%
  clean_names %>%
  filter(genome_name %in% pull(btp_nsti, ASV))

## number of distinct methanogens = 37
## from seven archaeal classes
btp_ko_methanogenesis_modules %>%
  select(genome_name) %>%
  distinct %>%
  summarise(n = n())

## ANOVA
methanogen_abundance_anova <- aov(Abundance ~ psf_type * depth, data = methanogen_abundance)
summary(methanogen_abundance_anova)

## pairwise comparisons
emmeans(methanogen_abundance_anova, pairwise ~ psf_type | depth)$contrasts

## methanogen abundance figure
methanogen_abundance_figure <- methanogen_abundance %>%
  ggplot(., aes(x = depth, y = Abundance, color = psf_type)) +
  geom_point(position = position_dodge(width = 0.8), shape = 1, size = 4, stroke = 1.4) +
  geom_signif(
    xmin = c(0.8), 
    xmax = c(1.2), 
    y_position = c(10), 
    annotations = c("ns"), 
    size = 0.5, 
    textsize = 6, 
    vjust = -0.3,
    family = "Helvetica",
    color = "black"
  ) +
  geom_signif(
    xmin = c(1.8, 2.8), 
    xmax = c(2.2, 3.2), 
    y_position = c(16, 19),
    annotations = c("*", "***"), 
    size = 0.5, 
    textsize = 8, 
    vjust = 0.3,
    family = "Helvetica",
    color = "black"
  ) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_x_discrete(labels = c("Surface", "Mid", "Deep")) +
  xlab("Depth") +
  ylab("Relative Abundance (%)") +
  guides(shape = guide_legend(title = "Depth")) + 
  scale_color_manual(
    "PSF", 
    values = c("#e66101", "#1a9641"), 
    labels = c("Burnt", "Intact")
  )

## abundance of methanotrophs
methanotroph_abundance <- btp_fire_ps %>%
  subset_taxa(., Kingdom == "Bacteria") %>%
  transform_sample_counts(., function(x) x / sum(x)) %>%
  psmelt %>%
  filter(OTU %in% pull(distinct(btp_ko_methanotrophy_modules, genome_name))) %>%
  group_by(Sample, psf_type, depth) %>%
  reframe(Abundance = round((sum(Abundance) * 100), 2))

## ANOVA
methanotroph_abundance_anova <- aov(Abundance ~ psf_type * depth, data = methanotroph_abundance)
summary(methanotroph_abundance_anova)

## pairwise comparisons
emmeans(methanotroph_abundance_anova, pairwise ~ psf_type | depth)$contrasts

## abundance table for heatmap
methanotroph_abundance_tab <- btp_fire_ps %>%
  subset_taxa(., Kingdom == "Bacteria") %>%
  transform_sample_counts(., function(x) x / sum(x)) %>%
  psmelt %>%
  filter(OTU %in% pull(distinct(btp_ko_methanotrophy_modules, genome_name))) %>%
  group_by(psf_type, depth, OTU, Family) %>%
  reframe(Abundance = sqrt(mean(Abundance))) %>%
  arrange(psf_type, depth, Family) %>%
  group_by(psf_type, depth) %>%
  mutate(group = cur_group_id()) %>%
  ungroup %>%
  select(OTU, group, Abundance) %>% 
  mutate(OTU = parse_number(OTU)) %>%
  pivot_wider(names_from = OTU, values_from = Abundance) %>%
  column_to_rownames("group") %>%
  as.matrix

## extract family labels for heatmap
methanotroph_family <- enframe(colnames(methanotroph_abundance_tab), value = "ASV", name = NULL) %>%
  mutate(ASV = paste0("ASV_", ASV)) %>%
  inner_join(
    tax_table(btp_fire_ps) %>%
      data.frame() %>%
      rownames_to_column("ASV")
  ) %>%
  mutate(ASV = parse_number(ASV)) %>%
  mutate(Family = factor(Family)) %>%
  pull(Family)

## plot methanotroph distribution
ht_opt$ROW_ANNO_PADDING <- unit(0.1, "cm")
options(repr.plot.width = 12, repr.plot.height = 2.8)
methanotroph_distribution_figure <- draw(
  Heatmap(methanotroph_abundance_tab,
    col = c("#ffffff", "#000000"),
    row_gap = unit(4, "mm"),
    column_gap = unit(4, "mm"),
    rect_gp = gpar(col = "black", lwd = 1.5),
    show_heatmap_legend = FALSE,
    show_row_dend = T,
    row_names_side = "right",
    column_dend_gp = gpar(lwd = 1.5),
    heatmap_legend_param = list(
      title = "Z-score",
      legend_direction = "horizontal",
      legend_width = unit(40, "mm"),
      legend_height = unit(55, "mm"),
      labels_gp = gpar(fontsize = 12),
      title_gp = gpar(fontsize = 13, fontface = "bold")
    ),
    border_gp = gpar(col = "black"),
    row_split = factor(rep(letters[1:2], each = 3), levels = letters[1:2]),
    column_split = methanotroph_family,
    column_title = NULL,
    cluster_rows = FALSE,
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    clustering_distance_columns = "spearman",
    row_title = NULL,
    show_column_names = TRUE,
    top_annotation = HeatmapAnnotation(
      pathway = anno_block(
        gp = gpar(lwd = 1.5, fill = c("#ffffb3", "#8dd3c7", rep("#fb8072", 5)))),
      empty = anno_empty(border = FALSE, height = unit(1.6, "mm"))
    ),
    left_annotation = rowAnnotation(
      psf_type = anno_block(
        gp = gpar(col = 0),
        labels = c("Burnt", "Intact"),
        labels_gp = gpar(col = "Black", fontsize = 18)
      )
    )
  ) +
    rowAnnotation(
      labels = anno_text(
        c("Surface", "Mid", "Deep", "Surface", "Mid", "Deep"),
        which = "row",
        gp = gpar(fontsize = 16)
      )
    ),
  heatmap_legend_side = "bottom"
)

## export plot
pdf(file.path(figures, "1_main", "3b_methanotroph_heatmap.pdf"), width = 12, height = 2.8)
methanotroph_distribution_figure
dev.off()

## methanogen abundance figure
methanotroph_abundance_figure <- methanotroph_abundance %>%
  ggplot(., aes(x = depth, y = Abundance, color = psf_type)) +
  geom_point(position = position_dodge(width = 0.8), shape = 1, size = 4, stroke = 1.4) +
  geom_signif(
    xmin = c(0.8), 
    xmax = c(1.2), 
    y_position = c(7), 
    annotations = c("ns"), 
    size = 0.5, 
    textsize = 6, 
    vjust = -0.3,
    family = "Helvetica",
    color = "black"
  ) +
  geom_signif(
    xmin = c(1.8, 2.8), 
    xmax = c(2.2, 3.2), 
    y_position = c(13, 8),
    annotations = c("*", "*"), 
    size = 0.5, 
    textsize = 8, 
    vjust = 0.3,
    family = "Helvetica",
    color = "black"
  ) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_x_discrete(labels = c("Surface", "Mid", "Deep")) +
  xlab("Depth") +
  ylab("Relative Abundance (%)") +
  guides(shape = guide_legend(title = "Depth")) + 
  scale_color_manual(
    "PSF", 
    values = c("#e66101", "#1a9641"), 
    labels = c("Burnt", "Intact")
  )

## methanotroph:methanogen ratio
methanotroph_methanogen_ratio <- btp_fire_ps %>%
  transform_sample_counts(., function(x) x / sum(x)) %>%
  psmelt() %>%
  filter(OTU %in% c(
    pull(distinct(btp_ko_methanotrophy_modules, genome_name)),
    pull(distinct(btp_ko_methanogenesis_modules, genome_name))
  )) %>%
  group_by(Sample, Kingdom, psf_type, depth) %>%
  reframe(Abundance = sum(Abundance)) %>%
  pivot_wider(names_from = Kingdom, values_from = Abundance) %>%
  filter(Archaea > 0) %>%
  mutate(ratio = round((Bacteria / Archaea), 2))

## ANOVA
methanotroph_methanogen_ratio_anova <- aov(ratio ~ psf_type * depth, data = methanotroph_methanogen_ratio)
summary(methanotroph_methanogen_ratio_anova)

## pairwise comparisons
emmeans(methanotroph_methanogen_ratio_anova, pairwise ~ psf_type | depth)$contrasts

## methanogen abundance figure
ratio_figure <- methanotroph_methanogen_ratio %>%
  ggplot(., aes(x = depth, y = ratio, color = psf_type)) +
  geom_point(position = position_dodge(width = 0.8), shape = 1, size = 4, stroke = 1.4) +
  geom_signif(
    xmin = c(0.8), 
    xmax = c(1.2), 
    y_position = c(16), 
    annotations = c("ns"), 
    size = 0.5, 
    textsize = 6, 
    vjust = -0.3,
    family = "Helvetica",
    color = "black"
  ) +
  geom_signif(
    xmin = c(1.8, 2.8), 
    xmax = c(2.2, 3.2), 
    y_position = c(35, 18),
    annotations = c("**", "*"), 
    size = 0.5, 
    textsize = 8, 
    vjust = 0.3,
    family = "Helvetica",
    color = "black"
  ) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.title.y = element_text(size = 15.5),
  ) +
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(labels = c("Surface", "Mid", "Deep")) +
  xlab("Depth") +
  ylab("Methanotrophs:Methanogens") +
  guides(shape = guide_legend(title = "Depth")) + 
  scale_color_manual(
    "PSF", 
    values = c("#e66101", "#1a9641"), 
    labels = c("Burnt", "Intact")
  )

## compose plot
options(repr.plot.width = 12, repr.plot.height = 4)
methanogen_abundance_figure + methanotroph_abundance_figure + ratio_figure

## export plot
pdf(file.path(figures, "1_main", "3b_methanogens_methanotrophs.pdf"), width = 12, height = 4)
methanogen_abundance_figure + methanotroph_abundance_figure + ratio_figure
dev.off()

## import shotgun metagenome marker profiles
mcra_profile <- read.csv(file.path(data, "7_metagenomes", "2020_btp_mcrA_graftM.csv"))
rpsj_profile <- read.csv(file.path(data, "7_metagenomes", "2020_btp_rpsJ_graftM.csv"))
pmoa_profile <- read.csv(file.path(data, "7_metagenomes", "2020_btp_pmoA_graftM.csv"))

## mcrA aggregated counts
mcra_counts <- mcra_profile %>%
  pivot_longer(-c("MSV", "consensus_lineage"), names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  reframe(mcra_counts = sum(counts))

## rpsJ aggregated counts for archaea
rpsj_counts_archaea <- rpsj_profile %>%
  pivot_longer(-c("MSV", "consensus_lineage"), names_to = "sample", values_to = "counts") %>%
  filter(str_detect(consensus_lineage, "Archaea")) %>%
  group_by(sample) %>%
  reframe(rpsj_counts = sum(counts))

## rpsJ aggregated counts for bacteria
rpsj_counts_bacteria <- rpsj_profile %>%
  pivot_longer(-c("MSV", "consensus_lineage"), names_to = "sample", values_to = "counts") %>%
  filter(str_detect(consensus_lineage, "Bacteria")) %>%
  group_by(sample) %>%
  reframe(rpsj_counts = sum(counts))

## pmoA aggregated counts
## retain only pmoA hits
pmoa_counts <- pmoa_profile %>%
  pivot_longer(-c("MSV", "consensus_lineage"), names_to = "sample", values_to = "counts") %>%
  filter(str_detect(consensus_lineage, "pmoA")) %>%
  group_by(sample) %>%
  reframe(pmoa_counts = sum(counts))

## plot correlation between mcrA and predicted methanogen counts
methanogens_predicted_sequenced_correlation_figure <- inner_join(mcra_counts, rpsj_counts_archaea) %>% 
  inner_join(methanogen_abundance, by = c("sample" = "Sample")) %>%
  mutate(norm_mcra_counts = (mcra_counts / rpsj_counts) * 100) %>%
  ggplot(., aes(x = Abundance, y = norm_mcra_counts)) +
  geom_point(size = 5, aes(color = psf_type, shape = depth), stroke = 1.3) +
  scale_shape_manual("Depth", values = c(0, 1, 3), labels = c("Surface", "Mid", "Deep")) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(aspect.ratio = 1, legend.position = "none") +
  geom_smooth(method = "lm", alpha = 0.15) +
  scale_color_manual("Fire", values = c("#e66101", "#1a9641"), labels = c("Burnt", "Intact")) +
  xlab("Predicted Relative Abundance (%)") +
  ylab("Sequenced Relative Abundance (%)") +
  annotate("text", x = 2.5, y = 53, label = "rho = 0.88", size = 7) +
  annotate("text", x = 3.2, y = 48, label = 'italic("p") < 0.001', size = 7, parse = TRUE)

## plot correlation between pmoA and predicted methanotroph counts
methanotrophs_predicted_sequenced_correlation_figure <- inner_join(pmoa_counts, rpsj_counts_bacteria) %>% 
  inner_join(methanotroph_abundance, by = c("sample" = "Sample")) %>%
  mutate(norm_pmoa_counts = (pmoa_counts / rpsj_counts) * 100) %>%
  ggplot(., aes(x = Abundance, y = norm_pmoa_counts)) +
  geom_point(size = 5, aes(color = psf_type, shape = depth), stroke = 1.3) +
  scale_shape_manual("Depth", values = c(0, 1, 3), labels = c("Surface", "Mid", "Deep")) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(legend.spacing = unit(1, "cm")) +
  geom_smooth(method = "lm", alpha = 0.15) +
  scale_color_manual("Fire", values = c("#e66101", "#1a9641"), labels = c("Burnt", "Intact")) +
  xlab("Predicted Relative Abundance (%)") +
  ylab("Sequenced Relative Abundance (%)") +
  annotate("text", x = 1.7, y = 6.5, label = "rho = 0.63", size = 7) +
  annotate("text", x = 1.1, y = 6, label = "italic(p)", size = 7, parse = TRUE) +
  annotate("text", x = 2.2, y = 6.05, label = "= 0.02", size = 7)

## compose plot
options(repr.plot.width = 12, repr.plot.height = 6)
methanogens_predicted_sequenced_correlation_figure + methanotrophs_predicted_sequenced_correlation_figure

## export plot
pdf(file.path(figures, "2_supplementary", "2_predicted_vs_sequenced.pdf"), width = 12, height = 6)
methanogens_predicted_sequenced_correlation_figure + methanotrophs_predicted_sequenced_correlation_figure
dev.off()