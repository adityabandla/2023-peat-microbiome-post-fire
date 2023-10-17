## name: 3_btp_diversity
## purpose: to compute, analyse and visualise
## biodiversity metrics for archael and bacterial communities
## author: Dr Aditya Bandla
## email: adityabandla@u.nus.edu

## load packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(emmeans)
  library(ggsignif)
  library(Hmisc)
  library(janitor)
  library(patchwork)
  library(phyloseq)
  library(scales)
  library(tidyverse)
  library(vegan)
})

## paths to directories
repo <- file.path("/Users/abandla/Desktop/2_research/1_manuscripts/2_2020_brunei_tropical_peat")
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

## export hellinger-transformed sequence counts table
## for primer-E analysis
btp_fire_ps %>%
  subset_taxa(., Kingdom == "Bacteria") %>%
  transform_sample_counts(., function(x) sqrt(x / sum(x))) %>%
  otu_table(.) %>%
  as.matrix() %>%
  t() %>%
  data.frame() %>%
  mutate(across(where(is.numeric), as.character)) %>%
  bind_rows(
    sample_data(btp_fire_ps) %>%
      data.frame() %>%
      select(psf_type, depth, plot, distance_from_canal) %>%
      t() %>%
      data.frame
  ) %>%
  write.csv(
    file.path(constructs, "2020_btp_bacteria_seqtab.csv"), 
    row.names = TRUE
  )

## export hellinger-transformed sequence counts table
## for primer-E analysis
btp_fire_ps %>%
  subset_taxa(., Kingdom != "Bacteria") %>%
  transform_sample_counts(., function(x) sqrt(x / sum(x))) %>%
  otu_table(.) %>%
  as.matrix() %>%
  t() %>%
  data.frame() %>%
  mutate(across(where(is.numeric), as.character)) %>%
  bind_rows(
    sample_data(btp_fire_ps) %>%
      data.frame() %>%
      select(psf_type, depth, plot, distance_from_canal) %>%
      t() %>%
      data.frame
  ) %>%
  write.csv(
    file.path(constructs, "2020_btp_archaea_seqtab.csv"), 
    row.names = TRUE
  )

## compute archea:bacteria ratio
archaea_to_bacteria_ratio <- btp_fire_ps %>%
  tax_glom("Kingdom") %>%
  transform_sample_counts(., function(x) x / sum(x)) %>%
  psmelt() %>%
  select(Sample, Abundance, psf_type, depth, Kingdom) %>%
  pivot_wider(
    id_cols = c("Sample", "psf_type", "depth"),
    names_from = "Kingdom",
    values_from = "Abundance"
  ) %>%
  group_by(Sample) %>%
  mutate(ratio = Archaea / Bacteria) %>%
  ungroup

## ANOVA
archaea_to_bacteria_ratio_anova <- aov(ratio ~ psf_type * depth, data = archaea_to_bacteria_ratio)
summary(archaea_to_bacteria_ratio_anova)

## pairwise comparisons
archaea_to_bacteria_ratio_anova_pairs <- emmeans(archaea_to_bacteria_ratio_anova, pairwise ~ psf_type | depth)$contrasts
archaea_to_bacteria_ratio_anova_pairs

## plot archaea:bacteria ratio
archaea_to_bacteria_ratio_figure <- archaea_to_bacteria_ratio %>%
  group_by(psf_type, depth) %>%
  ggplot(., aes(x = depth, y = ratio, color = psf_type)) +
  geom_point(position = position_dodge(width = 0.8), shape = 1, size = 4, stroke = 1.25) +
  geom_signif(
    xmin = c(0.8, 1.8), 
    xmax = c(1.2, 2.2), 
    y_position = c(0.3, 0.5), 
    annotations = c("ns", "ns"), 
    size = 0.5, 
    textsize = 6, 
    vjust = -0.3,
    family = "Helvetica",
    color = "black"
  ) +
  geom_signif(
    xmin = 2.8, 
    xmax = 3.2, 
    y_position = 1.05, 
    annotations = c("***"), 
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
  scale_y_continuous(limits = c(0, 1.15)) +
  scale_x_discrete(labels = c("Surface", "Mid", "Deep")) +
  xlab("Depth") +
  ylab("Archaea:Bacteria ratio") +
  guides(shape = guide_legend(title = "Depth (cm)")) + 
  scale_color_manual(
    "PSF", 
    values = c("#e66101", "#1a9641"), 
    labels = c("Burnt", "Intact")
  )

## compute shannon diversity for archaea
shannon_diversity_archaea <- btp_fire_ps %>%
  subset_taxa(., Kingdom == "Archaea") %>%
  estimate_richness(., measures = "Shannon") %>%
  rownames_to_column("sample_name") %>%
  inner_join(
    .,
    sample_data(btp_fire_ps) %>%
      data.frame() %>%
      rownames_to_column("sample_name")
  ) %>%
  select(sample_name, Shannon, psf_type, depth)

## ANOVA
shannon_diversity_archaea_anova <- aov(Shannon ~ psf_type * depth, data = shannon_diversity_archaea)
summary(shannon_diversity_archaea_anova)

## plot archaeal diversity
shannon_diversity_archaea_figure <- shannon_diversity_archaea %>%
  group_by(depth, psf_type) %>%
  ggplot(., aes(x = depth, y = Shannon, color = psf_type)) +
  geom_point(position = position_dodge(width = 0.8), shape = 1, size = 4, stroke = 1.25) +
  geom_signif(
    xmin = c(0.8, 1.8, 2.8), 
    xmax = c(1.2, 2.2, 3.2), 
    y_position = c(5, 5, 5),  
    annotations = c("ns", "ns", "ns"), 
    size = 0.5, 
    vjust = -0.3,
    textsize = 6, 
    family = "Helvetica",
    color = "black",
    tip_length = 0.1
  ) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  ) +
  xlab("Depth") +
  ylab("Archaeal diversity") +
  scale_y_continuous(limits = c(0, 8)) +
  scale_x_discrete(labels = c("Surface", "Mid", "Deep")) +
  guides(
    shape = guide_legend(title = "Depth")
  ) + 
  scale_color_manual(
    "PSF", 
    values = c("#e66101", "#1a9641"), 
    labels = c("Burnt", "Intact")
  )

## compute shannon diversity for bacteria
shannon_diversity_bacteria <- btp_fire_ps %>%
  subset_taxa(., Kingdom != "Archaea") %>%
  estimate_richness(., measures = "Shannon") %>%
  rownames_to_column("sample_name") %>%
  inner_join(
    .,
    sample_data(btp_fire_ps) %>%
      data.frame() %>%
      rownames_to_column("sample_name")
  ) %>%
  select(sample_name, Shannon, psf_type, depth)

## ANOVA
shannon_diversity_bacteria_anova <- aov(Shannon ~ psf_type * depth, data = shannon_diversity_bacteria)
summary(shannon_diversity_bacteria_anova)

## pairwise comparisons
shannon_diversity_bacteria_anova_pairs <- emmeans(shannon_diversity_bacteria_anova, pairwise ~ psf_type | depth)$contrasts
shannon_diversity_bacteria_anova_pairs

## plot bacterial diversity
shannon_diversity_bacteria_figure <- shannon_diversity_bacteria %>%
  group_by(depth, psf_type) %>%
  ggplot(., aes(x = depth, y = Shannon, color = psf_type)) +
  geom_point(position = position_dodge(width = 0.8), shape = 1, size = 4, stroke = 1.25) +
  geom_signif(
    xmin = c(0.8, 1.8), 
    xmax = c(1.2, 2.2), 
    y_position = c(7.3, 7.3),  
    annotations = c("ns", "ns"), 
    size = 0.5,
    vjust = -0.3,
    textsize = 6, 
    family = "Helvetica",
    color = "black",
    tip_length = 0.07
  ) +
  geom_signif(
    xmin = 2.8, 
    xmax = 3.2, 
    y_position = 7.3,  
    annotations = c("***"), 
    size = 0.5,
    vjust = 0.3,
    textsize = 8, 
    family = "Helvetica",
    color = "black",
    tip_length = 0.07
  ) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  ) +
  xlab("Depth") +
  ylab("Bacterial diversity") +
  scale_y_continuous(limits = c(0, 8)) +
  scale_x_discrete(labels = c("Surface", "Mid", "Deep")) +
  guides(
    shape = guide_legend(title = "Depth")
  ) + 
  scale_color_manual(
    "PSF", 
    values = c("#e66101", "#1a9641"), 
    labels = c("Burnt", "Intact")
  )

## function to compute ordination ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

## archaea constrained analysis of principle coordinates
archaea_cap <- ordinate(
  btp_fire_ps %>%
    subset_taxa(., Kingdom == "Archaea") %>%
    transform_sample_counts(., function(x) sqrt(x / sum(x))),
  "CAP", "bray", ~psf_type
)

## constrained analysis of principle coordinates ANOVA
anova.cca(archaea_cap)

## merge with metadata
archaea_cap_df <- archaea_cap %>%
  scores(.) %>%
  .$sites %>%
  data.frame() %>%
  cbind(sample_data(btp_fire_ps)) %>%
  group_by(psf_type, depth) %>%
  mutate(index = factor(cur_group_id()))

## estimate centroids
archaea_cap_centroids <- archaea_cap_df %>%
  summarise(CAP1 = mean(CAP1), MDS1 = mean(MDS1))

## extract data to plot
archaea_cap_ellipse <- data.frame()
for (g in levels(archaea_cap_df$index)) {
  archaea_cap_ellipse <- rbind(archaea_cap_ellipse, cbind(as.data.frame(with(
    archaea_cap_df[archaea_cap_df$index == g, ],
    veganCovEllipse(
      cov.wt(cbind(CAP1, MDS1),
        wt = rep(1 / length(CAP1), length(CAP1)))$cov,
      center = c(mean(CAP1), mean(MDS1))
    )
  )),
    index = g
  ))
}

## plot constrained analysis of principle coordinates
archaea_cap_figure <- ggplot(data = archaea_cap_ellipse) +
  geom_polygon(
    aes(x = CAP1, y = MDS1, group = index),
    fill = "#f0f0f0",
    alpha = 0.9
  ) +
  geom_point(
    data = archaea_cap_centroids,
    aes(x = CAP1, y = MDS1, color = psf_type, shape = depth),
    size = 5,
    stroke = 1.25
  ) +
  coord_fixed(clip = "off") +
  btp_theme +
  theme(legend.position = "none") +
  annotate("text", label = "Archaea", x = 1.12, y = 0.95, size = 6) +
  scale_color_manual("Fire", values = c("#e66101", "#1a9641"), labels = c("Burnt", "Intact")) +
  scale_shape_manual("Depth", values = c(0, 1, 3)) +
  scale_y_continuous(limits = c(-1.5, 1), breaks = c(-1, 0, 1)) +
  scale_x_continuous(limits = c(-2.2, 1.5)) +
  guides(shape = guide_legend(title = "Depth")) +
  xlab("CAP1 [10.7%]") +
  ylab("MDS1 [38.9%]")

## bacteria constrained analysis of principle coordinates
bacteria_cap <- ordinate(
  btp_fire_ps %>%
    subset_taxa(., Kingdom == "Bacteria") %>%
    transform_sample_counts(., function(x) sqrt(x / sum(x))),
  "CAP", "bray", ~psf_type
)

## constrained analysis of principle coordinates ANOVA
anova.cca(bacteria_cap)

## merge with metadata
bacteria_cap_df <- bacteria_cap %>%
  scores(.) %>%
  .$sites %>%
  data.frame() %>%
  cbind(sample_data(btp_fire_ps)) %>%
  group_by(psf_type, depth) %>%
  mutate(index = factor(cur_group_id()))

## estimate centroids
bacteria_cap_centroids <- bacteria_cap_df %>%
  summarise(CAP1 = mean(CAP1), MDS1 = mean(MDS1))

## extract data to plot
bacteria_cap_ellipse <- data.frame()
for (g in levels(bacteria_cap_df$index)) {
  bacteria_cap_ellipse <- rbind(bacteria_cap_ellipse, cbind(as.data.frame(with(
    bacteria_cap_df[bacteria_cap_df$index == g, ],
    veganCovEllipse(
      cov.wt(cbind(CAP1, MDS1),
        wt = rep(1 / length(CAP1), length(CAP1)))$cov,
      center = c(mean(CAP1), mean(MDS1))
    )
  )),
    index = g
  ))
}

## plot constrained analysis of principle coordinates
bacteria_cap_figure <- ggplot(data = bacteria_cap_ellipse) +
  geom_polygon(
    aes(x = CAP1, y = MDS1, group = index),
    fill = "#f0f0f0",
    alpha = 0.9
  ) +
  geom_point(
    data = bacteria_cap_centroids,
    aes(x = CAP1, y = MDS1, color = psf_type, shape = depth),
    size = 5,
    stroke = 1.25
  ) +
  coord_fixed(clip = "off") +
  btp_theme +
  theme(legend.position = "none") +
  annotate("text", label = "Bacteria", x = 1.12, y = 0.95, size = 6) +
  scale_color_manual("Fire", values = c("#e66101", "#1a9641"), labels = c("Burnt", "Intact")) +
  scale_shape_manual("Depth", values = c(0, 1, 3)) +
  scale_y_continuous(limits = c(-1.5, 1), breaks = c(-1, 0, 1)) +
  scale_x_continuous(limits = c(-2.4, 1.5)) +
  guides(shape = guide_legend(title = "Depth")) +
  xlab("CAP1 [9.9%]") +
  ylab("MDS1 [27%]")

## compose plot
options(repr.plot.width = 10, repr.plot.height = 7.5)
(archaea_to_bacteria_ratio_figure + shannon_diversity_archaea_figure + shannon_diversity_bacteria_figure) /
(archaea_cap_figure + bacteria_cap_figure) + plot_layout(heights = c(3, 4))

## export plot
pdf(file.path(figures, "1_main", "1a_diversity.pdf"), width = 10, height = 7.5)
(archaea_to_bacteria_ratio_figure + shannon_diversity_archaea_figure + shannon_diversity_bacteria_figure) /
(archaea_cap_figure + bacteria_cap_figure) + plot_layout(heights = c(3, 4))
dev.off()
