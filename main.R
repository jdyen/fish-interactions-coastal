# VEFMAP fish population modelling:
#   Analysis of fish data to determine strength of interactions
#   among species or groups of species in coastal-draining or upper
#   reaches of MDB rivers (focus on River Blackfish)
#
# Author: Jian Yen (jian.yen [at] deeca.vic.gov.au)
#
# Last updated: 26 April 2024

# flags for slow steps
reload_data <- FALSE
sample_again <- FALSE

# need some packages
library(qs)
library(aae.db)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(brms)
library(ggplot2)
library(bayesplot)
library(patchwork)

# load helper scripts
source("R/data.R")
source("R/lookup.R")
source("R/utils.R")

# load data sets
cpue <- fetch_fish(target = "southern", recompile = reload_data)

# filter out comparison surveys
cpue <- cpue |> filter(category != "comparison")

# and some groups where some size categories do not occur
cpue_sum <- cpue |>
  group_by(species, category) |>
  summarise(catch = sum(catch)) |>
  ungroup() |>
  mutate(include = catch >= 10) |>
  select(-catch)
cpue <- cpue |>
  left_join(cpue_sum, by = c("species", "category")) |>
  filter(include) |>
  select(-include)

# collapse to a single survey per site and year
cpue <- cpue |>
  select(waterbody, id_site, survey_year, species, category, gear_type, reach_no, catch, effort_h) |>
  group_by(waterbody, id_site, survey_year, species, category, gear_type, reach_no) |>
  summarise(
    catch = sum(catch), 
    effort_h = sum(effort_h)
  ) |>
  ungroup()

# add previous year's catch, filter to non-missing survey pairs
cpue <- cpue |>
  mutate(survey_year_minus1 = survey_year - 1) |>
  left_join(
    cpue |> select(waterbody, reach_no, id_site, survey_year, species, category, gear_type, catch, effort_h),
    by = c("waterbody", "reach_no", "id_site", "species", "category", "gear_type", "survey_year_minus1" = "survey_year"),
    suffix = c("", "_ym1")
  ) |>
  mutate(
    log_cpue_tm1 = log(catch_ym1 + 1) - log(effort_h_ym1),
    log_effort_h = log(effort_h)
  )

# filter to paired obs only
cpue <- cpue |> filter(!is.na(log_cpue_tm1))

# add a combined species_category variable to avoid `:` in model call
cpue <- cpue |> mutate(species_category = paste(species, category, sep = "_"))

# fit model
if (sample_again) {
  
  # sample from brms model
  stan_seed <- 2024-03-13
  mod <- brm(
    bf(
      catch ~
        log_cpue_tm1 +
        species + category +
        (1 | waterbody + waterbody:reach_no + id_site) +
        (-1 + species | waterbody + waterbody:reach_no + id_site) +
        (-1 + category | waterbody + waterbody:reach_no + id_site) +
        (-1 + species_category | waterbody + waterbody:reach_no + id_site) +
        (1 | survey_year + waterbody:survey_year + 
           waterbody:reach_no:survey_year + 
           gear_type) +
        offset(log_effort_h),
      shape ~ (
        1 | species + species_category +
          waterbody + waterbody:reach_no +
          gear_type + 
          species:waterbody + species:waterbody:reach_no +
          category:waterbody +
          species_category:waterbody +
          species:gear_type
      )
    ),
    data = cpue,
    family = negbinomial(),
    chains = 4,
    cores = 4,
    seed = stan_seed,
    iter = 3000, 
    warmup = 1000,
    control = list(adapt_delta = 0.8, max_treedepth = 10),
    backend = "rstan",
    threads = threading(2),
    init_r = 1.0
  )
  
  # save to file
  qsave(mod, file = "outputs/fitted/draws.qs")
  
} else {
  
  # load saved version
  mod <- qread("outputs/fitted/draws.qs")
  
}

# only repeat these steps if it's a newly fitted model
if (sample_again) {
  
  # basic diagnostics
  pp_plot <- pp_check(mod, group = "species", type = "dens_overlay_grouped") +
    scale_x_log10() +
    theme(plot.background = element_rect(fill = "white"))
  pp_max <- bayesplot::pp_check(
    mod, group = "species", type = "stat_grouped", stat = "max"
  ) +
    theme(plot.background = element_rect(fill = "white"))
  pp_pzero_grouped <- bayesplot::pp_check(
    mod, 
    group = "species",
    type = "stat_grouped", 
    stat = \(x) mean(x == 0)
  ) +
    theme(plot.background = element_rect(fill = "white"))
  rhat <- brms::rhat(mod)
  neff <- brms::neff_ratio(mod)
  standard_diag <- tibble(
    stat = c(rep("Rhat", length(rhat)), rep("Neff ratio", length(neff))),
    par = c(names(rhat), names(neff)),
    value = c(rhat, neff)
  ) |>
    ggplot(aes(x = value)) +
    geom_histogram() +
    xlab("Value") +
    ylab("Count") +
    facet_wrap( ~ stat, scales = "free")
  
  # combine into a single plot
  hmc_diagnostics <- 
    (pp_plot + theme(legend.position = "none")) /
    (pp_max + scale_x_log10() + theme(legend.position = "none")) /
    (pp_pzero_grouped + theme(legend.position = "none")) +
    patchwork::plot_annotation(tag_levels = "a")
  
  # save these
  ggsave(
    filename = "outputs/figures/all-diagnostics.png",
    plot = hmc_diagnostics,
    device = ragg::agg_png,
    width = 8.5, 
    height = 13,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    filename = "outputs/figures/diagnostics.png",
    plot = standard_diag,
    device = ragg::agg_png,
    width = 6, 
    height = 6,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  # model fit
  fitted_vals <- posterior_epred(mod)
  fitted_vals <- t(fitted_vals)
  colnames(fitted_vals) <- paste("pred", seq_len(ncol(fitted_vals)), sep = "_")
  r2 <- mod$data |>
    select(catch, species, category, waterbody) |>
    mutate(id = seq_len(nrow(mod$data))) |>
    left_join(fitted_vals |> as_tibble() |> mutate(id = seq_len(nrow(fitted_vals))), by = "id") |>
    pivot_longer(cols = contains("pred"), names_prefix = "pred_", names_to = "iter") |>
    group_by(species, waterbody, category, iter) |>
    summarise(cor = cor(catch, value) ^ 2) |>
    group_by(species, category, waterbody) |>
    summarise(
      cor_mean = mean(cor, na.rm = TRUE),
      cor_mid = median(cor, na.rm = TRUE),
      cor_lower = quantile(cor, probs = 0.1, na.rm = TRUE),
      cor_upper = quantile(cor, probs = 0.9, na.rm = TRUE)
    )
  write.csv(r2, file = "outputs/tables/r2.csv")
  r2_all <- mod$data |>
    mutate(est = apply(fitted_vals, 1, median)) |>
    group_by(species) |>
    summarise(cor = cor(catch, est) ^ 2) |>
    ungroup()
  write.csv(r2_all, file = "outputs/tables/r2-all.csv")
  
}

# extract variance components so we can calculate correlations to plot
#    (at three levels and over a hierarchy of waterbodies, reaches, and sites)
var_comp <- VarCorr(mod, probs = c(0.1, 0.25, 0.75, 0.9))

# grab correlations and clean them up
terms_corr <- c(
  "waterbody",
  "waterbody:reach_no",
  "id_site"
)
idx_sp <- c(2:9)
idx_cat <- c(10:12)
idx_all <- c(13:31)
cor_sp <- lapply(var_comp[terms_corr], extract_corr, idx = idx_sp)
cor_cat <- lapply(var_comp[terms_corr], extract_corr, idx = idx_cat)
cor_all <- lapply(var_comp[terms_corr], extract_corr, idx = idx_all)
cor_sp <- mapply(
  \(x, y) x |> mutate(hier = y),
  x = cor_sp,
  y = names(cor_sp), 
  SIMPLIFY = FALSE
)
cor_cat <- mapply(
  \(x, y) x |> mutate(hier = y),
  x = cor_cat,
  y = names(cor_cat), 
  SIMPLIFY = FALSE
)
cor_all <- mapply(
  \(x, y) x |> mutate(hier = y),
  x = cor_all,
  y = names(cor_all), 
  SIMPLIFY = FALSE
)

# flatten these and average over each level/hierarchy to compare values
correlation_level_hier <- bind_rows(
  cor_sp |> bind_rows() |>
    filter(value != 1) |> 
    mutate(level = "species"),
  cor_cat |> bind_rows() |> 
    filter(value != 1) |>
    mutate(level = "category"),
  cor_all |> bind_rows() |> 
    filter(value != 1) |> 
    mutate(level = "species_category")
) |>
  group_by(level, hier) |>
  summarise(
    q10 = quantile(value, probs = 0.1),
    q90 = quantile(value, probs = 0.9),
    q25 = quantile(value, probs = 0.25),
    q75 = quantile(value, probs = 0.75),
    q50 = median(value)
  ) |>
  ungroup()

# and plot the values across levels and hierarchy
corr_level_hier_plot <- correlation_level_hier |>
  mutate(
    level = factor(
      level, 
      levels = rev(c("category", "species", "species_category", "observed")),
      labels = rev(c("Size class", "Species", "Species x size class", "Observed data"))
    ),
    hier = factor(
      hier,
      levels = rev(c("waterbody", "waterbody:reach_no", "id_site")),
      labels = rev(c("Waterbody", "Reach", "Site"))
    )
  ) |>
  ggplot(aes(x = q50, xmin = q10, xmax = q90, y = level, col = hier)) +
  geom_pointrange(size = 1, linewidth = 1, position = position_dodge(0.4)) +
  geom_pointrange(aes(xmin = q25, xmax = q75), linewidth = 2, position = position_dodge(0.4)) +
  ylab("") + xlab("Pairwise correlations") +
  scale_color_brewer(palette = "Set2", name = "Level") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
  filename = "outputs/figures/correlations-at-hier.png",
  plot = corr_level_hier_plot,
  device = ragg::agg_png,
  width = 6.5, 
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# grab the correlations to plot at site level (finest level,
#   most relevant to the pop models and species interactions)
corr_sp_hier_plot <- cor_sp$id_site |>
  mutate(
    name1 = gsub("species", "", name1),
    name2 = gsub("species", "", name2),
    name1 = gsub("Gadopsis", "Gadopsis ", name1),
    name1 = gsub("Galaxias", "Galaxias ", name1),
    name1 = gsub("Carp", "Carp/ ", name1),
    name1 = gsub("Nannoperca", "Nannoperca ", name1),
    name1 = gsub("^Perca", "Perca ", name1),
    name1 = gsub("Pseudaphritis", "Pseudaphritis ", name1),
    name1 = gsub("Salmo", "Salmo ", name1),
    name2 = gsub("Gadopsis", "Gadopsis ", name2),
    name2 = gsub("Galaxias", "Galaxias ", name2),
    name2 = gsub("Carp", "Carp/ ", name2),
    name2 = gsub("Nannoperca", "Nannoperca ", name2),
    name2 = gsub("^Perca", "Perca ", name2),
    name2 = gsub("Pseudaphritis", "Pseudaphritis ", name2),
    name2 = gsub("Salmo", "Salmo ", name2),
    value = ifelse(name1 == name2, NA, value),
    value_short = ifelse(pos_y > pos_x, round(value, 2), "")
  ) |>
  ggplot(aes(x = name1, y = name2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value_short)) +
  xlab("") + ylab("") +
  scale_fill_distiller(
    palette = "RdBu", 
    direction = 1, 
    limits = c(-1, 1), 
    name = "Pairwise correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "bottom"
  )
corr_cat_hier_plot <- cor_cat$id_site |>
  mutate(
    name1 = gsub("category", "", name1),
    name2 = gsub("category", "", name2),
    name1 = paste0(toupper(substr(name1, 1, 1)), tolower(substr(name1, 2, nchar(name1)))),
    name2 = paste0(toupper(substr(name2, 1, 1)), tolower(substr(name2, 2, nchar(name2)))),
    value = ifelse(name1 == name2, NA, value),
    value_short = ifelse(pos_y > pos_x, round(value, 2), "")
  ) |>
  ggplot(aes(x = name1, y = name2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value_short)) +
  xlab("") + ylab("") +
  scale_fill_distiller(
    palette = "RdBu", 
    direction = 1, 
    limits = c(-1, 1), 
    name = "Pairwise correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "bottom"
  )
corr_all_hier_plot <- cor_all$id_site |>
  mutate(
    name1 = gsub("species_category", "", name1),
    name2 = gsub("species_category", "", name2),
    name1 = gsub("Gadopsis", "Gadopsis ", name1),
    name1 = gsub("Galaxias", "Galaxias ", name1),
    name1 = gsub("Carp", "Carp/ ", name1),
    name1 = gsub("Nannoperca", "Nannoperca ", name1),
    name1 = gsub("^Perca", "Perca ", name1),
    name1 = gsub("Pseudaphritis", "Pseudaphritis ", name1),
    name1 = gsub("Salmo", "Salmo ", name1),
    name2 = gsub("Gadopsis", "Gadopsis ", name2),
    name2 = gsub("Galaxias", "Galaxias ", name2),
    name2 = gsub("Carp", "Carp/ ", name2),
    name2 = gsub("Nannoperca", "Nannoperca ", name2),
    name2 = gsub("^Perca", "Perca ", name2),
    name2 = gsub("Pseudaphritis", "Pseudaphritis ", name2),
    name2 = gsub("Salmo", "Salmo ", name2),
    name1 = gsub("_", " (", name1),
    name1 = paste0(name1, ")"),
    name2 = gsub("_", " (", name2),
    name2 = paste0(name2, ")"),
    value = ifelse(name1 == name2, NA, value),
    value_short = ifelse(pos_y > pos_x, round(value, 2), "")
  ) |>
  ggplot(aes(x = name1, y = name2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value_short)) +
  xlab("") + ylab("") +
  scale_fill_distiller(
    palette = "RdBu", 
    direction = 1, 
    limits = c(-1, 1), 
    name = "Pairwise correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "bottom"
  )

# combine the species and category plots into a single plot
corr_sp_size_plot <- (corr_sp_hier_plot + theme(legend.position = "none")) /
  corr_cat_hier_plot +
  plot_annotation(tag_levels = "a")

# save these to file
ggsave(
  filename = "outputs/figures/correlations-species-size.png",
  plot = corr_sp_size_plot,
  device = ragg::agg_png,
  width = 8, 
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  filename = "outputs/figures/correlations-all.png",
  plot = corr_all_hier_plot,
  device = ragg::agg_png,
  width = 10, 
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

# and pull out the variance intercepts for reference
terms <- c(
  "waterbody",
  "waterbody:reach_no",
  "id_site",
  "survey_year",
  "survey_year:waterbody",
  "survey_year:waterbody:reach_no",
  "gear_type"
)
var_int <- lapply(var_comp[terms], \(x) x$sd["Intercept", ])
var_int <- bind_rows(var_int) |> mutate(level = terms)
var_int_plot <- var_int |>
  mutate(
    level = factor(
      level, 
      levels = rev(
        c(
          "waterbody", "waterbody:reach_no", "id_site", "survey_year",
          "survey_year:waterbody", "survey_year:waterbody:reach_no",
          "gear_type"
        )
      ),
      labels = rev(
        c("Waterbody", "Reach", "Site", "Year", "Year x Waterbody",
          "Year x Reach", "Gear type"))
    )
  ) |>
  ggplot(aes(x = Estimate, xmin = Q10, xmax = Q90, y = level)) +
  geom_pointrange(size = 1, linewidth = 1) +
  geom_pointrange(aes(xmin = Q25, xmax = Q75), linewidth = 2) +
  ylab("") + xlab("Standard deviation") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
  filename = "outputs/figures/var-components.png",
  plot = var_int_plot,
  device = ragg::agg_png,
  width = 6.5, 
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)
