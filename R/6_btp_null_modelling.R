## name: 6_btp_null_modelling.R
## purpose: to estimate relative importance of assembly processes governing
## archaeal & bacterial communities
## author: Dr Aditya Bandla
## email: adityabandla@u.nus.edu

## start of code run on NSCC (National Super Computing Center, Singapore)
## load packages
suppressPackageStartupMessages({
  library(ape)
  library(iCAMP)
  library(NST)
  library(phyloseq)
  library(tidyverse)
  library(vegan)
})

## set parallel option for vegan
options(mc.cores = 48)

## import phyloseq object
btp_fire_ps <- readRDS("2020_btp_fire_ps.rds")
btp_fire_ps

## bacterial community
comm <- btp_fire_ps %>%
  subset_taxa(., Kingdom == "Bacteria") %>%
  otu_table(.) %>%
  data.frame()

## taxonomic classifications
clas <- btp_fire_ps %>%
  subset_taxa(., Kingdom == "Bacteria") %>%
  tax_table(.) %>%
  data.frame()

## choose an archaeal ASV at random to 
## root the bacterial tree
btp_fire_ps %>%
  subset_taxa(., Kingdom == "Archaea") %>%
  tax_table() %>%
  data.frame %>%
  sample_n(1)

## ASV_248 was randomly chosen as the archaeal 
## tip for rooting the tree after rooting, this tip was pruned
tree <- read.tree("2020_btp_fire_asv_seqs_aligned.tre") %>%
  keep.tip(., c(colnames(comm), "ASV_248")) %>%
  root(., "ASV_248", resolve.root = TRUE) %>%
  drop.tip(., "ASV_248")

## import environmental data and log 
## transform all variables except pH
env <- read.csv(
  "2020_btp_depth_env_data.csv",
  header = TRUE,
  row.names = 1,
  as.is = TRUE,
  stringsAsFactors = FALSE
) %>%
  select(-one_of("psf_type", "depth"))

## perform some basic checks to see if inputs are consistent 
## with each other. We will also check if any ASVs are zero across all samples 
## and drop them.
sampid_check <- match.name(rn.list = list(comm = comm))
comm <- sampid_check$comm
comm <- comm[, colSums(comm) > 0, drop = FALSE]

spid_check <- match.name(cn.list = list(comm = comm), rn.list = list(clas = clas), tree.list = list(tree = tree))
comm <- spid_check$comm
clas <- spid_check$clas
tree <- spid_check$tree

## compute between-ASVs phylogenetic distances
nworker <- 48
memory.G <- 192
if (!file.exists(file.path("1_pd", "1_bacteria", "pd.desc"))) {
  pdist.big(
    tree = tree,
    wd = file.path(getwd(), "1_pd", "1_bacteria"),
    nworker = nworker,
    memory.G = memory.G,
    output = TRUE,
    pd.spname.file = "pd_taxon_name.csv"
  ) %>%
    write.csv(.,
      file.path(
        getwd(), "1_pd", "1_bacteria", "2020_btp_bacteria_pd.csv"
      ),
      quote = FALSE, 
      row.names = FALSE
    )
} else {
  pd_big <- list()
  pd_big$tip.label <- read.csv(
    file.path(getwd(), "1_pd", "1_bacteria", "pd_taxon_name.csv"),
    row.names = 1,
    stringsAsFactors = FALSE
  )[, 1]
  pd_big$pd.wd <- file.path(getwd(), "1_pd/1_bacteria")
  pd_big$pd.file <- "pd.desc"
  pd_big$pd.name.file <- "pd_taxon_name.csv"
}

## import phylogenetic distance matrix for 
## constructing mantel correlograms
pd_bacteria <- read.csv(file.path("1_pd", "1_bacteria", "2020_btp_bacteria_pd.csv"))
rownames(pd_bacteria) <- colnames(pd_bacteria)

## compute differences in niche optima as the absolute difference between niche optima 
## between every species-pair. The niche optima for each species is computed 
## as the abundance-weighted mean of each log-transformed (except pH) environmental 
## variable. The input needs to be a TSS normalised counts table.
niche_dif <- dniche(
  env = env, 
  comm = (comm / rowSums(comm)), 
  method = "niche.value",
  nworker = nworker, 
  out.dist = FALSE, 
  bigmemo = TRUE,
  nd.wd = file.path(getwd(), "2_niche_diff", "1_bacteria")
)

## log values do not seem to impact Mantel R or its trend
## mantel correlogram analysis for water quality parameters
## mantel correlogram uses progressive correction for p-values
## so later tests are penalised more than the earlier ones
## we are only interested mainly in the distance classes 0-0.5
## pH
mantel_correlog_ph <- mantel.correlog(
  niche_dif$nd$pH[rownames(pd_bacteria), rownames(pd_bacteria)], 
  pd_bacteria,
  cutoff = FALSE,
  n.clas = 96
)

## function to construct mantel correlograms by looping through
## variables stored in the niche optima list
construct_correlograms <- function(niche_dif, pd_bacteria, cutoff = FALSE, n.clas = 96) {
  results <- list()
  for (var in names(niche_dif$nd)) {
    cor_result <- mantel.correlog(niche_dif$nd[[var]][rownames(pd_bacteria), rownames(pd_bacteria)],
      pd_bacteria,
      cutoff = cutoff,
      n.clas = n.clas
    )
    results[[var]] <- cor_result
  }
  return(results)
}

## explore thresholds for phylogenetic signal and 
## bin size 
ds <- 0.2
bin.size.limit <- 24

## perform phylogenetic binning
phylobin <- taxa.binphy.big(
  tree = tree, 
  pd.desc = pd_big$pd.file, 
  pd.spname = pd_big$tip.label,
  pd.wd = pd_big$pd.wd, 
  ds = ds,
  bin.size.limit = bin.size.limit,
  nworker = nworker
)

## trim data
## abcut = 20 retains 77% ASVs, abcut = 10 retains 92% ASVs
sp.bin <- phylobin$sp.bin[, 3, drop = FALSE]
sp.ra <- colMeans(comm / rowSums(comm))
abcut <- 0
commc <- comm[, colSums(comm) >= abcut, drop = FALSE]
dim(commc)
spname.use <- colnames(commc)

## test within-bin phylogenetic signal
binps <- ps.bin(
  sp.bin = sp.bin,
  sp.ra = sp.ra,
  spname.use = spname.use,
  pd.desc = pd_big$pd.file,
  pd.spname = pd_big$tip.label,
  pd.wd = pd_big$pd.wd,
  nd.list = niche_dif$nd,
  nd.spname = niche_dif$names,
  ndbig.wd = niche_dif$nd.wd,
  cor.method = c("pearson"),
  r.cut = 0.05,
  p.cut = 0.05,
  min.spn = 5
)

## infer community assembly mechanism by phylogenetic-bin-based 
## null model analysis
sig.index <- "Confidence"
rand.time = 1000
prefix = "2020_btp_bacteria"
icres <- icamp.big(
  comm = comm,
  pd.desc = pd_big$pd.file,
  pd.spname = pd_big$tip.label,
  pd.wd = pd_big$pd.wd,
  rand = rand.time, 
  tree = tree,
  prefix = prefix,
  ds = ds, 
  pd.cut = NA,
  sp.check = TRUE,
  phylo.rand.scale = "within.bin",
  taxa.rand.scale = "across.all",
  phylo.metric = "bMPD",
  sig.index = sig.index,
  bin.size.limit = 24,
  nworker = nworker,
  memory.G = memory.G,
  rtree.save = FALSE,
  detail.save = TRUE,
  qp.save = FALSE,
  detail.null = TRUE,
  ignore.zero = TRUE,
  output.wd = getwd(),
  correct.special = TRUE,
  unit.sum = rowSums(comm),
  special.method = "depend",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975,
  omit.option = "no",
  meta.ab = NULL
)

## import sample grouping information
treat <- btp_fire_ps %>%
  sample_data(.) %>%
  data.frame() %>%
  mutate(psf_type_depth = paste0(psf_type, "_", depth))

## summarise iCAMP results for each bin
icbin <- icamp.bins(
  icamp.detail = icres$detail, 
  treat = select(treat, psf_type_depth),
  clas = clas, 
  silent = FALSE, 
  boot = TRUE,
  rand.time = rand.time, 
  between.group = FALSE
)

## export results
write.csv(
  icbin$Pt,
  file = paste0(prefix, "_process_importance_each_group.csv"),
  row.names = FALSE
)

write.csv(
  icbin$Ptk,
  file = paste0(prefix, "_process_importance_each_bin_each_group.csv"),
  row.names = FALSE
)

write.csv(
  icbin$Ptuv,
  file = paste0(prefix, "_process_importance_each_turnover.csv"),
  row.names = FALSE
)

write.csv(
  icbin$BPtk,
  file = paste0(prefix, "_bin_contribute_to_process_each_group.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    ID = rownames(icbin$Class.Bin),
    icbin$Class.Bin,
    stringsAsFactors = FALSE
  ),
  file = paste0(prefix, "_taxon_bin.csv"),
  row.names = FALSE
)

write.csv(
  icbin$Bin.TopClass,
  file = paste0(prefix, "_bin_top_taxon.csv"),
  row.names = FALSE
)

## extract results for bootstrap
icamp_result <- icres$CbMPDiCBraya
icamp_result_group <- icamp_result %>%
  filter(str_detect(sample1, "_B_")) %>%
  filter(str_detect(sample2, "_B_"))

## define grouping for bootstrap test
## loop through different depths manually
group <- treat %>%
  filter(depth == "95-100")

## perform bootstrap
## for every depth
icboot_deep <- icamp.boot(
  icamp.result = icamp_result_group,
  treat = select(group, psf_type),
  rand.time = rand.time,
  compare = TRUE,
  silent = FALSE,
  between.group = FALSE,
  ST.estimation = TRUE
)

## archaeal community
## 94% tips are retained for testing phylogenetic signal
comm <- btp_fire_ps %>%
  subset_taxa(., Kingdom == "Archaea") %>%
  subset_samples(., depth != "0-5") %>%
  filter_taxa(., function(x) sum(x > 0) >= 2, TRUE) %>%
  otu_table(.) %>%
  data.frame()

## taxonomic classifications
clas <- btp_fire_ps %>%
  subset_taxa(., Kingdom == "Archaea") %>%
  subset_samples(., depth != "0-5") %>%
  filter_taxa(., function(x) sum(x > 0) >= 2, TRUE) %>%
  tax_table(.) %>%
  data.frame()

## choose an bacterial ASV at random to 
## root the archaeal tree
btp_fire_ps %>%
  subset_taxa(., Kingdom == "Bacteria") %>%
  tax_table() %>%
  data.frame %>%
  sample_n(1)

## ASV_5891 was randomly chosen as the bacterial 
## tip for rooting the tree after rooting, this tip was pruned
tree <- read.tree("2020_btp_fire_asv_seqs_aligned.tre") %>%
  keep.tip(., c(colnames(comm), "ASV_5891")) %>%
  root(., "ASV_5891", resolve.root = TRUE) %>%
  drop.tip(., "ASV_5891")

## import environmental data 
env <- read.csv(
  "2020_btp_depth_env_data.csv",
  header = TRUE,
  row.names = 1,
  as.is = TRUE,
  stringsAsFactors = FALSE
) %>%
  filter(depth != "0-5") %>%
  select(-one_of("psf_type", "depth")) 

## perform some basic checks to see if inputs are consistent 
## with each other. We will also check if any ASVs are zero across all samples 
## and drop them.
sampid_check <- match.name(rn.list = list(comm = comm))
comm <- sampid_check$comm
comm <- comm[, colSums(comm) > 0, drop = FALSE]

spid_check <- match.name(cn.list = list(comm = comm), rn.list = list(clas = clas), tree.list = list(tree = tree))
comm <- spid_check$comm
clas <- spid_check$clas
tree <- spid_check$tree

## construct phylogenetic distance matrix for
## mantel correlograms
## delete the backing files so that they don't clash 
## with the files generated with the whole community data
pd_big <- pdist.big(
  tree = tree,
  wd = file.path(getwd(), "1_pd", "2_archaea"),
  nworker = nworker,
  memory.G = memory.G,
  output = TRUE,
  pd.spname.file = "pd_taxon_name.csv"
) %>%
  write.csv(.,
    file.path(
      getwd(), "1_pd", "2_archaea", "2020_btp_archaea_pd.csv"
    ),
    quote = FALSE,
    row.names = FALSE
  )

## import phylogenetic distance matrix for 
## constructing mantel correlograms
pd_archaea <- read.csv(file.path("1_pd", "2_archaea", "2020_btp_archaea_pd.csv"))
rownames(pd_archaea) <- colnames(pd_archaea)

## compute differences in niche optima as the absolute difference between niche optima 
## between every species-pair. The niche optima for each species is computed 
## as the abundance-weighted mean of each log-transformed (except pH) environmental 
## variable. The input needs to be a TSS normalised counts table.
niche_dif <- dniche(
  env = env, 
  comm = (comm / rowSums(comm)), 
  method = "niche.value",
  nworker = nworker, 
  out.dist = TRUE, 
  bigmemo = FALSE,
  nd.wd = file.path(getwd(), "2_niche_diff", "2_archaea")
)

## function to construct mantel correlograms by looping through
## variables stored in the niche optima list
construct_correlograms <- function(niche_dif, pd_archaea, cutoff = FALSE, n.clas = 75) {
  results <- list()
  for (var in names(niche_dif$nd)) {
    cor_result <- mantel.correlog(niche_dif$nd[[var]][rownames(pd_archaea), rownames(pd_archaea)],
      pd_archaea,
      cutoff = cutoff,
      n.clas = n.clas
    )
    results[[var]] <- cor_result
  }
  return(results)
}

## construct correlograms
archaea_correlog <- construct_correlograms(niche_dif, pd_archaea)

## explore thresholds for phylogenetic signal and 
## bin size 
ds <- 0.2
bin.size.limit <- 48

## perform phylogenetic binning
phylobin <- taxa.binphy.big(
  tree = tree, 
  pd.desc = pd_big$pd.file, 
  pd.spname = pd_big$tip.label,
  pd.wd = pd_big$pd.wd, 
  ds = ds,
  bin.size.limit = bin.size.limit,
  nworker = nworker
)

## trim data
sp.bin <- phylobin$sp.bin[, 3, drop = FALSE]
sp.ra <- colMeans(comm / rowSums(comm))
abcut <- 0
commc <- comm[, colSums(comm) >= abcut, drop = FALSE]
dim(commc)
spname.use <- colnames(commc)

## test within-bin phylogenetic signal
binps <- ps.bin(
  sp.bin = sp.bin,
  sp.ra = sp.ra,
  spname.use = spname.use,
  pd.desc = pd_big$pd.file,
  pd.spname = pd_big$tip.label,
  pd.wd = pd_big$pd.wd,
  nd.list = niche_dif$nd,
  nd.spname = niche_dif$names,
  ndbig.wd = niche_dif$nd.wd,
  cor.method = c("pearson"),
  r.cut = 0.05,
  p.cut = 0.05,
  min.spn = 5
)

## infer community assembly mechanism by phylogenetic-bin-based 
## null model analysis
sig.index <- "Confidence"
rand.time = 1000
ds = 0.2
prefix = "2020_btp_archaea"
icres <- icamp.big(
  comm = comm,
  pd.desc = pd_big$pd.file,
  pd.spname = pd_big$tip.label,
  pd.wd = pd_big$pd.wd,
  rand = rand.time, 
  tree = tree,
  prefix = prefix,
  ds = ds, 
  pd.cut = NA,
  sp.check = TRUE,
  phylo.rand.scale = "within.bin",
  taxa.rand.scale = "across.all",
  phylo.metric = "bMPD",
  sig.index = sig.index,
  bin.size.limit = 12,
  nworker = nworker,
  memory.G = memory.G,
  rtree.save = FALSE,
  detail.save = TRUE,
  qp.save = FALSE,
  detail.null = TRUE,
  ignore.zero = TRUE,
  output.wd = getwd(),
  correct.special = TRUE,
  unit.sum = rowSums(comm),
  special.method = "depend",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975,
  omit.option = "no",
  meta.ab = NULL
)

## import sample grouping information
treat <- btp_fire_ps %>%
  sample_data(.) %>%
  data.frame() %>%
  mutate(psf_type_depth = paste0(psf_type, "_", depth))

## summarise iCAMP results for each bin
icbin <- icamp.bins(
  icamp.detail = icres$detail, 
  treat = select(treat, psf_type_depth),
  clas = clas, 
  silent = FALSE, 
  boot = TRUE,
  rand.time = rand.time, 
  between.group = FALSE
)

## export results
write.csv(
  icbin$Pt,
  file = paste0(prefix, "_process_importance_each_group.csv"),
  row.names = FALSE
)

write.csv(
  icbin$Ptk,
  file = paste0(prefix, "_process_importance_each_bin_each_group.csv"),
  row.names = FALSE
)

write.csv(
  icbin$Ptuv,
  file = paste0(prefix, "_process_importance_each_turnover.csv"),
  row.names = FALSE
)

write.csv(
  icbin$BPtk,
  file = paste0(prefix, "_bin_contribute_to_process_each_group.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    ID = rownames(icbin$Class.Bin),
    icbin$Class.Bin,
    stringsAsFactors = FALSE
  ),
  file = paste0(prefix, "_taxon_bin.csv"),
  row.names = FALSE
)

write.csv(
  icbin$Bin.TopClass,
  file = paste0(prefix, "_bin_top_taxon.csv"),
  row.names = FALSE
)

## extract results for bootstrap
icamp_result <- icres$CbMPDiCBraya
icamp_result_group <- icamp_result %>%
  filter(str_detect(sample1, "_B_")) %>%
  filter(str_detect(sample2, "_B_"))

## define grouping for bootstrap test
## loop through different depths manually
group <- treat %>%
  filter(depth == "95-100")

## perform bootstrap
## for every depth
icboot_deep <- icamp.boot(
  icamp.result = icamp_result_group,
  treat = select(group, psf_type),
  rand.time = rand.time,
  compare = TRUE,
  silent = FALSE,
  between.group = FALSE,
  ST.estimation = TRUE
)

## load packages
suppressPackageStartupMessages({
  library(emmeans)
  library(ggh4x)
  library(ggpattern)
  library(ggsignif)
  library(Hmisc)
  library(janitor)
  library(patchwork)
  library(phyloseq)
  library(scales)
  library(tidyverse)
})

## paths to directories
repo <- file.path("/Users/abandla/Desktop/2_research/1_manuscripts/2_2020_brunei_peat_fire")
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

## import correlograms
correlog_bacteria <- readRDS(file.path(data, "6_icamp", "3_correlograms", "230714_mantel_correlog_bacteria.rds"))
correlog_archaea <- readRDS(file.path(data, "6_icamp", "3_correlograms", "230715_mantel_correlog_archaea.rds"))

## import phyloseq object
btp_fire_ps <- readRDS(file.path(data, "3_phyloseq", "2020_btp_fire_ps.rds"))
btp_fire_ps

## unnest for plotting
correlog_bacteria_df <- do.call(rbind, lapply(names(correlog_bacteria), function(var_name) {
  correlog_bacteria[[var_name]]$mantel.res %>%
    data.frame() %>%
    mutate(variable = var_name)
}))

## unnest for plotting
correlog_archaea_df <- do.call(rbind, lapply(names(correlog_archaea), function(var_name) {
  correlog_archaea[[var_name]]$mantel.res %>%
    data.frame() %>%
    mutate(variable = var_name)
}))

## plot correlograms
correlog_bacteria_figure <- correlog_bacteria_df %>%
  filter(variable %in% c("pH", "water_temperature", "EC", "DO")) %>%
  filter(!is.na(Pr.corrected.)) %>%
  mutate(
    sig = case_when(
      Pr.corrected. <= 0.05 ~ "0.05",
      Pr.corrected. > 0.05 & Pr.corrected. <= 0.1 ~ "0.1",
      TRUE ~ "ns"
    )
  ) %>%
  ggplot(., aes(x = class.index, y = Mantel.cor)) +
  geom_point(aes(fill = sig), shape = 21, size = 3) +
  geom_line() +
  facet_wrap(~variable,
    nrow = 1,
    labeller = as_labeller(c(
      "DO" = "DO",
      "EC" = "EC",
      "pH" = "pH",
      "water_temperature" = "Temperature"
    ))
  ) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    strip.text = element_text(size = 16),
    strip.background = element_blank(),
    panel.spacing.x = unit(2, "lines")
  ) +
  xlab("PD class index") +
  ylab("Mantel R") +
  scale_fill_manual(
    values = c("0.05" = "#000000", "ns" = "#ffffff", "0.1" = "#d3d3d3")
  ) +
  geom_hline(yintercept = 0, color = "red")

## plot correlograms
correlog_archaea_figure <- correlog_archaea_df %>%
  filter(variable %in% c("pH", "water_temperature", "EC", "DO")) %>%
  filter(!is.na(Pr.corrected.)) %>%
  mutate(
    sig = case_when(
      Pr.corrected. <= 0.05 ~ "0.05",
      Pr.corrected. > 0.05 & Pr.corrected. <= 0.1 ~ "0.1",
      TRUE ~ "ns"
    )
  ) %>%
  ggplot(., aes(x = class.index, y = Mantel.cor)) +
  geom_point(aes(fill = sig), shape = 21, size = 3) +
  geom_line() +
  facet_wrap(~factor(variable),
    nrow = 1,
    labeller = as_labeller(c(
      "DO" = "DO",
      "EC" = "EC",
      "pH" = "pH",
      "water_temperature" = "Temperature"
    ))
  ) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    strip.text = element_text(size = 16),
    strip.background = element_blank(),
    panel.spacing.x = unit(2, "lines")
  ) +
  xlab("PD class index") +
  ylab("Mantel R") +
  scale_fill_manual(
    values = c("0.05" = "#000000", "ns" = "#ffffff", "0.1" = "#d3d3d3")
  ) +
  geom_hline(yintercept = 0, color = "red")

## view plot
options(repr.plot.width = 12, repr.plot.height = 8)
correlog_bacteria_figure / correlog_archaea_figure

## export plot
pdf(file.path(figures, "2_sf_correlograms.pdf"), width = 12, height = 8)
correlog_bacteria_figure / correlog_archaea_figure
dev.off()

## import process bootstrap results
boot_bacteria <- readRDS(file.path(data, "6_icamp", "5_null_model", "1_bacteria", "230720_icboot_results_bacteria.rds"))
boot_archaea <- readRDS(file.path(data, "6_icamp", "5_null_model", "2_archaea", "230720_icboot_results_archaea.rds"))

## unnest
boot_bacteria_df <- do.call(rbind, lapply(names(boot_bacteria), function(var_name) {
  boot_bacteria[[var_name]]$summary %>%
    select(Group, Process, Observed, Mean, Stdev) %>%
    mutate(Mean = as.numeric(Mean), Stdev = as.numeric(Stdev)) %>%
    mutate(depth = var_name)
})) %>%
  mutate(domain = "Bacteria")

## unnest
boot_archaea_df <- do.call(rbind, lapply(names(boot_archaea), function(var_name) {
  boot_archaea[[var_name]]$summary %>%
    select(Group, Process, Observed, Mean, Stdev) %>%
    mutate(Mean = as.numeric(Mean), Stdev = as.numeric(Stdev)) %>%
    mutate(depth = var_name)
})) %>%
  mutate(domain = "Archaea")

## unnest
boot_bacteria_compare_df <- do.call(rbind, lapply(names(boot_bacteria), function(var_name) {
  boot_bacteria[[var_name]]$compare %>%
    mutate(depth = var_name)
})) %>%
  mutate(domain = "Bacteria")

## unnest
boot_archaea_compare_df <- do.call(rbind, lapply(names(boot_archaea), function(var_name) {
  boot_archaea[[var_name]]$compare %>%
    mutate(depth = var_name)
})) %>%
  mutate(domain = "Archaea")

## import bin contributions to processes
bin_process_bacteria <- read.csv(
  file.path(data, "6_icamp", "5_null_model", "1_bacteria", "2020_btp_bacteria_bin_contribute_to_process_each_group.csv")
)

## import bin contributions to processes
bin_process_archaea <- read.csv(
  file.path(data, "6_icamp", "5_null_model", "2_archaea", "2020_btp_archaea_bin_contribute_to_process_each_group.csv")
)

## combine and plot
processes_combined_figure <- bind_rows(boot_bacteria_df, boot_archaea_df) %>%
  filter(Process %nin% c("Heterogeneous.Selection", "Homogenizing.Dispersal")) %>%
  ggplot(., aes(x = depth, y = Mean, color = Group)) +
  geom_point(position = position_dodge(width = 0.8), size = 4, shape = 1, stroke = 1.1) +
  geom_errorbar(aes(ymin = Mean - Stdev, ymax = Mean + Stdev),
    width = 0.4,
    position = position_dodge(width = 0.8),
    linewidth = 0.8
  ) +
  geom_line(aes(group = Group), position = position_dodge(width = 0.8), linewidth = 0.8) +
  facet_grid2(
    domain ~ Process,
    scales = "free_y",
    independent = "y",
    labeller = as_labeller(
      default = label_wrap_gen(10), c(
        "Dispersal.Limitation" = "Dispersal Limitation",
        "Drift.and.Others" = "Drift",
        "Homogeneous.Selection" = "Homogeneous Selection",
        "Stochasticity" = "Stochasticity",
        "Bacteria" = "Bacteria",
        "Archaea" = "Archaea"
      )
    )
  ) +
  coord_cartesian(clip = "off") +
  btp_theme +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "#f0f0f0"),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_line(color = "#d3d3d3"),
    panel.grid.minor = element_line(color = "#d3d3d3"),
    strip.placement = "outside",
    axis.text.x = element_text(size = 15)
  ) +
  scale_x_discrete(labels = c("Surface", "Mid", "Deep")) +
  scale_color_manual(
    "PSF",
    values = c("#e66101", "#1a9641"),
    labels = c("Burnt", "Intact")
  ) +
  facetted_pos_scales(y = list(
    Process == "Dispersal.Limitation" ~ scale_y_continuous(limits = c(0, 0.8), n.breaks = 3),
    Process == "Drift.and.Others" ~ scale_y_continuous(limits = c(0, 0.8), n.breaks = 3),
    Process == "Heterogeneous.Selection" ~ scale_y_continuous(limits = c(0, 0.1), n.breaks = 3),
    Process == "Homogeneous.Selection" ~ scale_y_continuous(limits = c(0, 1.2), n.breaks = 3, breaks = c(0, 0.5, 1)),
    Process == "Homogenizing.Dispersal" ~ scale_y_continuous(limits = c(-0.005, 0.2), n.breaks = 3),
    Process == "Stochasticity" ~ scale_y_continuous(limits = c(0, 1), n.breaks = 3)
  )) +
  xlab("") +
  ylab("% Turnovers") +
  geom_signif(
    data = data.frame(
      Process = c("Dispersal.Limitation", "Homogeneous.Selection", "Stochasticity"),
      domain = c("Bacteria")
    ),
    aes(
      y_position = c(0.7, 0.85, 0.85),
      xmin = c(2.8, 2.8, 2.8),
      xmax = c(3.2, 3.2, 3.2),
      annotations = c("**", "*", "*")
    ),
    tip_length = 0.05,
    manual = T,
    textsize = 8,
    family = "Helvetica",
    inherit.aes = FALSE
  ) +
  geom_signif(
    data = data.frame(
      Process = c("Drift.and.Others", "Drift.and.Others", "Stochasticity", "Homogeneous.Selection"),
      domain = c("Archaea")
    ),
    aes(
      y_position = c(0.7, 0.7, 0.85, 1.05),
      xmin = c(0.8, 2.8, 2.8, 2.8),
      xmax = c(1.2, 3.2, 3.2, 3.2),
      annotations = c("*", "**", "*", "*")
    ),
    tip_length = 0.05,
    manual = T,
    textsize = 8,
    family = "Helvetica",
    inherit.aes = FALSE
  )

## view plot
options(repr.plot.width = 12, repr.plot.height = 6)
processes_combined_figure

## export figure
pdf(file.path(figures, "4_processes_community_level.pdf"), width = 12, height = 6)
processes_combined_figure
dev.off()

## community-level processes 
## relative importance for archaea
community_process_importance <- read.csv(
  file.path(
    data,
    "6_icamp",
    "5_null_model",
    "2_archaea",
    "2020_btp_archaea_process_importance_each_group.csv"
  ),
  header = TRUE
) %>%
  select(-one_of("Method", "GroupBasedOn")) %>%
  pivot_longer(-Group, names_to = "Index", values_to = "process_per") %>%
  mutate(process_per = as.numeric(process_per) * 100) %>%
  mutate(bin = "Community") %>%
  select(Group, Index, bin, process_per)

## process importance figure
process_importance_figure <- read.csv(
  file.path(
    data,
    "6_icamp",
    "5_null_model",
    "2_archaea",
    "2020_btp_archaea_process_importance_each_bin_each_group.csv"
  )
) %>%
  filter(!str_detect(Index, "Dominant")) %>%
  select(-one_of("Method", "GroupBasedOn")) %>%
  pivot_longer(-c("Group", "Index"), names_to = "bin", values_to = "process_per") %>%
  mutate(process_per = as.numeric(process_per) * 100) %>%
  bind_rows(community_process_importance) %>%
  filter(str_detect(Group, "_95-100")) %>%
  separate(Group, c("Group", "Depth"), "_") %>%
  mutate(dominant_process = case_when(
    Group == "Burnt" & bin == "Community" & Index == "HoS" ~ "Dominant",
    Group == "Burnt" & bin == "bin1" & Index == "HoS" ~ "Dominant",
    Group == "Burnt" & bin == "bin5" & Index == "HoS" ~ "Dominant",
    Group == "Intact" & bin == "bin4" & Index == "DR" ~ "Dominant",
    Group == "Intact" & bin == "bin5" & Index == "HoS" ~ "Dominant",
    TRUE ~ "NA"
  )) %>%
  ggplot(., aes(x = fct_rev(Group), y = process_per, fill = Index, pattern = dominant_process)) +
  geom_bar(stat = "identity", width = 1) +
  geom_bar_pattern(
    color = "black",
    width = 1,
    stat = "identity",
    pattern_color = "black",
    pattern_alpha = 0.2,
    pattern_spacing = 0.08,
    show.legend = TRUE
  ) +
  scale_pattern_manual(values = c("Dominant" = "stripe", "NA" = "none")) +
  facet_wrap(~bin,
    nrow = 3,
    labeller = as_labeller(c(
      "Community" = "Community",
      "bin1" = "Methanomicrobia",
      "bin2" = "Nitrososphaeria-1",
      "bin3" = "Nitrososphaeria-2",
      "bin4" = "Nitrososphaeria-3",
      "bin5" = "Bathyarchaeia",
      "bin6" = "Thermoplasmata"
    ))
  ) +
  coord_flip(clip = "off") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), labels = c("Burnt", "Intact")) +
  btp_theme +
  theme(
    legend.position = "none",
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(0, 0, -27, 0),
    strip.background = element_rect(fill = NA),
    strip.text.x = element_text(size = 16),
    panel.spacing.x = unit(2, "lines")
  ) +
  scale_fill_manual(
    "Process",
    values = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3"),
    labels = c("DL", "DR", "HD", "VS", "HS")
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(color = "black", size = 1.5)
    ),
    pattern = guide_legend("",
      nrow = 1, byrow = TRUE,
      override.aes = list(
        fill = NA,
        pattern = "stripe",
        color = "black",
        pattern_alpha = 0.4,
        pattern_spacing = 0.01,
        size = 1.5
      )
    )
  ) +
  xlab("") +
  ylab("% Turnovers")

## view plot
options(repr.plot.width = 11, repr.plot.height = 4)
process_importance_figure

## export figure
pdf(file.path(figures, "1_main", "5_process_importance_archaea.pdf"), width = 11, height = 4)
process_importance_figure
dev.off()
