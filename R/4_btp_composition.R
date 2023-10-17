## name: 4_btp_composition
## purpose: to analyse and visualise microbiome composition and
## to identify differentially abundant taxa
## author: Dr Aditya Bandla
## email: adityabandla@u.nus.edu

## load packages
suppressPackageStartupMessages({
  library(DESeq2)
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

## identify top five classes in 
## bacteria and archaea
top_five_microbial_classes <- btp_fire_ps %>%
  psmelt %>%
  group_by(Sample, Kingdom) %>%
  mutate(relative_abundance = (Abundance / sum(Abundance)) * 100) %>%
  group_by(Sample, Kingdom, Class, depth, psf_type) %>%
  reframe(class_relative_abundance = sum(relative_abundance)) %>%
  group_by(Class) %>% 
  mutate(
    mean_abundance = mean(class_relative_abundance),
    mean_prevalence = sum(class_relative_abundance > 0) / 24,
    F1 = 2 * ((mean_abundance * mean_prevalence) / (mean_abundance + mean_prevalence))
  ) %>%
  ungroup %>%
  group_by(Sample, Kingdom) %>%
  top_n(6, wt = F1) %>%
  ungroup

## top five bacterial classes
top_five_bacterial_classes <- top_five_microbial_classes %>%
  filter(Kingdom == "Bacteria" & !is.na(Class)) %>%
  group_by(psf_type, depth, Class, Kingdom) %>%
  reframe(mean_abundance = mean(class_relative_abundance)) %>%
  pivot_wider(names_from = "Class", values_from = mean_abundance) %>%
  mutate(Others = 100 - rowSums(select_if(., is.numeric))) %>%
  pivot_longer(-c("psf_type", "depth", "Kingdom"), names_to = "Class", values_to = "Abundance") %>%
    mutate(Class = factor(
    Class,
    levels = c(
      "Acidobacteriae",
      "Actinobacteria",
      "Alphaproteobacteria",
      "Gammaproteobacteria",
      "Planctomycetes",
      "Others"
    )
  ))

## top five archaeal classes
top_five_archaeal_classes <- top_five_microbial_classes %>%
  filter(Kingdom != "Bacteria" & !is.na(Class)) %>%
  group_by(psf_type, depth, Class, Kingdom) %>%
  reframe(mean_abundance = mean(class_relative_abundance)) %>%
  pivot_wider(names_from = "Class", values_from = mean_abundance) %>%
  mutate(Others = 100 - rowSums(select_if(., is.numeric))) %>%
  pivot_longer(-c("psf_type", "depth", "Kingdom"), names_to = "Class", values_to = "Abundance") %>%
  mutate(Class = factor(
    Class,
    levels = c(
      "Bathyarchaeia",
      "Methanomicrobia",
      "Methanosarcinia",
      "Nitrososphaeria",
      "Thermoplasmata",
      "Others"
    )
  ))

## plot top five archaeal classes
top_five_archaeal_classes_figure <- top_five_archaeal_classes %>%
  ggplot(., aes(x = depth, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  facet_wrap(~psf_type) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.spacing.x = unit(0.8, "lines"),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(size = 16),
    legend.key.size = unit(0.8, 'cm'),
    legend.spacing.y = unit(0.24, 'cm'),
    legend.spacing.x = unit(0.15, 'cm'),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  scale_x_discrete(expand = c(0, 0), labels = c("Surface", "Mid", "Deep")) +
  scale_y_continuous(expand = c(0, 0)) +
   scale_fill_manual("Class",
    values = c(
      "Bathyarchaeia" = "#fb8072",
      "Methanomicrobia" = "#80b1d3",
      "Methanosarcinia" = "#b3de69",
      "Nitrososphaeria" = "#fccde5",
      "Thermoplasmata" = "#ccebc5",
      "Others" = "#d9d9d9"
    )
  ) +
guides(fill = guide_legend(byrow = TRUE)) +
ylab("Relative Abundance (%)") +
xlab("Depth")

## plot top five bacterial classes
top_five_bacterial_classes_figure <- top_five_bacterial_classes %>%
  ggplot(., aes(x = depth, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  facet_wrap(~psf_type) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    axis.text.y = element_text(margin = margin(0, 10, 0, 0)),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.spacing.x = unit(0.8, "lines"),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(size = 16),
    legend.key.size = unit(0.8, 'cm'),
    legend.spacing.y = unit(0.24, 'cm'),
    legend.spacing.x = unit(0.15, 'cm'),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  scale_x_discrete(expand = c(0, 0), labels = c("Surface", "Mid", "Deep")) +
  scale_y_continuous(expand = c(0, 0)) +
   scale_fill_manual("Class",
    values = c(
      "Acidobacteriae" = "#8dd3c7",
      "Actinobacteria" = "#ffffb3",
      "Alphaproteobacteria" = "#bebada",
      "Gammaproteobacteria" = "#fdb462",
      "Planctomycetes" = "#bc80bd",
      "Others" = "#d9d9d9"
    )
  ) +
guides(fill = guide_legend(byrow = TRUE)) +
ylab("Relative Abundance (%)") +
xlab("Depth")

## view plot
options(repr.plot.width = 12, repr.plot.height = 4.7)
top_five_archaeal_classes_figure + top_five_bacterial_classes_figure

## save plot
pdf(file.path(figures, "1_main", "2_composition.pdf"), width = 12, height = 4.7)
top_five_archaeal_classes_figure + top_five_bacterial_classes_figure
dev.off()

## fire responders at the
## class level
btp_fire_class_dds <- btp_fire_ps %>%
  tax_glom(., "Class") %>%
  phyloseq_to_deseq2(., ~psf_type + depth + psf_type:depth) %>%
  DESeq(., "Wald")

## check reference levels
## no differentials at depths 0-5 & 35-40
resultsNames(btp_fire_class_dds)

## relevel
btp_fire_class_dds$depth <- relevel(btp_fire_class_dds$depth, "95-100")
btp_fire_class_dds <- nbinomWaldTest(btp_fire_class_dds)
resultsNames(btp_fire_class_dds)

## shrink effect sizes
btp_fire_class_depth_95_100 <- lfcShrink(
  dds = btp_fire_class_dds, 
  res = results(btp_fire_class_dds, contrast=c("psf_type","Intact","Burnt")),
  coef = "psf_type_Intact_vs_Burnt"
) %>%
  data.frame %>%
  filter(padj <= 0.05 & abs(log2FoldChange) >= 1) %>%
  rownames_to_column("ASV") %>%
  inner_join(
    tax_table(btp_fire_ps) %>%
    data.frame %>%
    rownames_to_column("ASV") %>%
    select(ASV, Kingdom, Phylum, Class)
  ) %>%
  arrange(-baseMean)

## save
write.csv(
  btp_fire_class_depth_95_100,
  file.path(data, "4_differentials", "2020_btp_fire_class_depth_95_100.csv"),
  row.names = FALSE,
  quote = FALSE
)
## fire responders at the
## ASV level
btp_fire_asv_dds <- btp_fire_ps %>%
  phyloseq_to_deseq2(., ~psf_type + depth + psf_type:depth) %>%
  DESeq(., "Wald")

## check reference levels
## no differentials at depth 35-40
resultsNames(btp_fire_asv_dds)

## relevel reference levels
btp_fire_asv_dds$depth <- relevel(btp_fire_asv_dds$depth, "95-100")
btp_fire_asv_dds <- nbinomWaldTest(btp_fire_asv_dds)
resultsNames(btp_fire_asv_dds)

## shrink effect sizes
btp_fire_asv_depth_0_5 <- lfcShrink(
  dds = btp_fire_asv_dds, 
  res = results(btp_fire_asv_dds, contrast=c("psf_type","Intact","Burnt")),
  coef = "psf_type_Intact_vs_Burnt"
) %>%
  data.frame %>%
  filter(padj <= 0.05 & abs(log2FoldChange) >= 1) %>%
  rownames_to_column("ASV") %>%
  inner_join(
    tax_table(btp_fire_ps) %>%
    data.frame %>%
    rownames_to_column("ASV") %>%
    select(ASV, Kingdom, Phylum, Class, Order, Family, Genus)
  ) %>%
  arrange(-baseMean)

## shrink effect sizes
btp_fire_asv_depth_95_100 <- lfcShrink(
  dds = btp_fire_asv_dds, 
  res = results(btp_fire_asv_dds, contrast=c("psf_type","Intact","Burnt")),
  coef = "psf_type_Intact_vs_Burnt"
) %>%
  data.frame %>%
  filter(padj <= 0.05 & abs(log2FoldChange) >= 1) %>%
  rownames_to_column("ASV") %>%
  inner_join(
    tax_table(btp_fire_ps) %>%
    data.frame %>%
    rownames_to_column("ASV") %>%
    select(ASV, Kingdom, Phylum, Class, Order, Family, Genus)
  ) %>%
  arrange(-baseMean)

## save
write.csv(
  btp_fire_asv_depth_95_100,
  file.path(data, "4_differentials", "2020_btp_fire_asv_depth_95_100.csv"),
  row.names = FALSE,
  quote = FALSE
)
