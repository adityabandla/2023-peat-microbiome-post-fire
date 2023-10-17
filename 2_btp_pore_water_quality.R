## name: 2_btp_pore_water_quality
## purpose: to analyse peat pore water quality data 
## author: Dr Aditya Bandla
## email: adityabandla@u.nus.edu

## load packages
suppressPackageStartupMessages({
  library(emmeans)
  library(ggh4x)
  library(Hmisc)
  library(janitor)
  library(phyloseq)
  library(scales)
  library(tidyverse)
})

## paths to directories
repo <- file.path("/Users/abandla/Desktop/2_research/1_manuscripts/2_2020_brunei_tropical_peat")
data <- file.path(repo, "1_data")
figures <- file.path(repo, "3_figures")

## define global theme options for plotting
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

## import environmental data
btp_env_depth <- read.csv(file.path(data, "1_metadata", "2020_btp_depth_env_data.csv"))

## summarise
btp_env_depth_figure <- btp_env_depth %>%
  pivot_longer(-c(sample, psf_type, depth), names_to = "variable", values_to = "value") %>%
  group_by(psf_type, variable, depth) %>%
  ggplot(., aes(x = depth, y = value, color = psf_type)) +
  geom_point(position = position_dodge(width = 0.8), shape = 1, size = 4, stroke = 1.3) +
  facet_wrap(
    ~variable,
    scales = "free",
    strip.position = "left",
    labeller = as_labeller(c(
      "DO" = "DO (mg/L)",
      "EC" = "EC (µS/cm)",
      "TDS" = "TDS (mg/L)",
      "pH" = "pH",
      "salinity" = "Salinity (‱)",
      "water_temperature" = "Temperature (°C)"
    ))
  ) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    strip.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    strip.background = element_rect(fill = NA),
    strip.placement = "outside",
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_line(color = "#d3d3d3"),
    panel.grid.minor = element_line(color = "#d3d3d3"),
  ) +
  xlab("Depth") +
  ylab("") +
  guides(
    shape = guide_legend(title = "Depth")
  ) +
  scale_color_manual(
    "PSF",
    values = c("#e66101", "#1a9641"),
    labels = c("Burnt", "Intact")
  ) +
  scale_x_discrete(labels = c("Surface", "Mid", "Deep")) +
  facetted_pos_scales(y = list(
    variable == "DO (mg/" ~ scale_y_continuous(limits = c(0, 5), n.breaks = 3),
    variable == "EC" ~ scale_y_continuous(limits = c(0, 200), n.breaks = 3),
    variable == "TDS" ~ scale_y_continuous(limits = c(0, 120), n.breaks = 3),
    variable == "pH" ~ scale_y_continuous(limits = c(0, 6), n.breaks = 3),
    variable == "salinity" ~ scale_y_continuous(limits = c(0, 0.1), n.breaks = 3),
    variable == "water_temperature" ~ scale_y_continuous(limits = c(24, 28), n.breaks = 3)
  ))

## view plot
options(repr.plot.width = 12, repr.plot.height = 8)
btp_env_depth_figure

## export figure
pdf(file.path(figures, "2_supplementary", "1_porewater_quality.pdf"), width = 12, height = 8)
btp_env_depth_figure
dev.off()

## ANOVA water temperature
water_temperature_anova <- aov(water_temperature ~ psf_type + depth, data = btp_env_depth)
summary(water_temperature_anova)

## pairwise comparisons
## reference: https://cran.r-project.org/web/packages/emmeans/vignettes/AQuickStart.html#additive
emmeans(water_temperature_anova, pairwise ~ depth)$contrasts

## ANOVA pH
ph_anova <- aov(pH ~ psf_type + depth, data = btp_env_depth)
summary(ph_anova)

## ANOVA EC
ec_anova <- aov(EC ~ psf_type + depth, data = btp_env_depth)
summary(ec_anova)

## ANOVA DO
do_anova <- aov(DO ~ psf_type + depth, data = btp_env_depth)
summary(do_anova)
