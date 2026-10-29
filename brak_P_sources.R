
# Brak / NUAR: where does the change in P (drain_effekt) between
# del v2 (4 June) and del v3 (28 July) come from?
#
# How P is built (Note "Parameterværdier for afstrømning", 03.08.2026, p. 4-8):
#   P per 100 m cell = mean P-factor over 22 years (2000-2021) of Daisy percolation for
#     (DMI 10 km weather grid) x (typejord = AC_jord, dominant in the cell) x (crop),
#   computed with one of two parameter sets: sand (formula 2) or clay (formula 3).
# So P should be ONE value per DMI_10km x AC_jord for each crop, and a change in P
# can only come from:
#   A. the soil profile given to the cell (AC_jord)           -> check: did AC_jord change?
#   B. the weather grid given to the cell (DMI_10km)          -> check: did DMI_10km change?
#   C. the sand/clay parameter set used for a profile         -> a whole profile moves, in all grids,
#                                                                the same direction
#   D. the Daisy percolation itself (weather, evaporation,    -> a change that differs from grid to grid
#      soil hydraulics), i.e. the lookup table was recomputed
#   E. cell-level filling / imputation                        -> P no longer one value per combination
#
# Sections
#   0. Read v2 and v3 (both crops), only the fields needed
#   1. A and B: did the cell keys change?
#   2. E: is P one value per DMI_10km x AC_jord in each version?
#   3. Lookup table per combination, and the change per combination
#   4. C vs D: how much of the change is explained by profile, by grid, or by both
#   5. Per profile: direction, size, and how uniform across grids (profile-level shift?)
#   6. Sand/clay: does the change line up with the cell's JB class (sand = JB 1-3, 11)?
#   7. Plots and map (per DMI grid)


library(tidyverse)
library(sf)
library(terra)
library(data.table)

brak  <- "C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak"
dlv   <- file.path(brak, "Delivered_AU")
out   <- file.path(brak, "compare_deliverables", "P_sources")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 10))

gdb <- c(v2 = file.path(dlv, "GeoDB_20260604/GeoDB_20260604/Geodatabase_20260604.gdb"),
         v3 = file.path(dlv, "GeoDB_20260728/GeoDB_20260728/Geodatabase_20260728.gdb"))
lab <- c(v2 = "del v2 (4 June)", v3 = "del v3 (28 July)")
crops <- c(SB = "Vaarbyg", WW = "VinterHvede")
keep  <- c("gridcode", "AC_jord", "DMI_10km", "Georegion", "JB2024", "drain_effekt")
sand_jb <- c(1, 2, 3, 11)        # sand parameter set (Note, p. 5)


# ---- 0. Read --------------------------------------------------------------------
rd <- function(g, crop) {
  lay <- st_layers(g)$name
  lay <- lay[str_detect(lay, crops[[crop]])]
  sql <- paste("SELECT", paste(keep, collapse = ", "), "FROM", lay)
  as_tibble(terra::vect(g, layer = lay, query = sql, what = "attributes"))
}

d <- list()
for (cr in names(crops)) {
  a <- rd(gdb[["v2"]], cr); b <- rd(gdb[["v3"]], cr)
  d[[cr]] <- a %>%
    inner_join(b, by = "gridcode", suffix = c("_v2", "_v3")) %>%
    filter(JB2024_v2 != -99) %>%
    mutate(crop = cr,
           P_ratio = drain_effekt_v3 / drain_effekt_v2,
           P_pct   = 100 * (P_ratio - 1),
           changed = abs(P_pct) > 0.1,
           sand    = JB2024_v3 %in% sand_jb)
  rm(a, b); gc()
}
d <- bind_rows(d)
count(d, crop)


# ---- 1. A and B: did the cell keys change? ----------------------------------------
keys <- d %>%
  group_by(crop) %>%
  summarise(cells = n(),
            AC_jord_changed  = sum(AC_jord_v2  != AC_jord_v3,  na.rm = TRUE),
            DMI_10km_changed = sum(DMI_10km_v2 != DMI_10km_v3, na.rm = TRUE),
            Georegion_changed = sum(Georegion_v2 != Georegion_v3, na.rm = TRUE),
            JB2024_changed   = sum(JB2024_v2   != JB2024_v3,   na.rm = TRUE),
            P_changed_pct    = 100 * mean(changed, na.rm = TRUE))
keys
write_csv(keys, file.path(out, "P1_cell_keys_changed.csv"))
# If AC_jord and DMI_10km are unchanged, sources A and B are ruled out.

d <- d %>% rename(AC_jord = AC_jord_v3, DMI_10km = DMI_10km_v3, Georegion = Georegion_v3,
                  JB2024 = JB2024_v3)


# ---- 2. E: one P per DMI_10km x AC_jord? -----------------------------------------
uniq <- d %>%
  group_by(crop, DMI_10km, AC_jord) %>%
  summarise(cells = n(),
            nP_v2 = n_distinct(round(drain_effekt_v2, 5)),
            nP_v3 = n_distinct(round(drain_effekt_v3, 5)),
            .groups = "drop")
uniq %>%
  group_by(crop) %>%
  summarise(combinations = n(),
            multi_v2 = sum(nP_v2 > 1), cells_multi_v2 = sum(cells[nP_v2 > 1]),
            multi_v3 = sum(nP_v3 > 1), cells_multi_v3 = sum(cells[nP_v3 > 1])) %>%
  print() %>%
  write_csv(file.path(out, "P2_one_value_per_combination.csv"))
# If a combination has more than one P, look at what else varies inside it
# (JB2024 sand/clay class? Georegion?):
inside <- d %>%
  semi_join(uniq %>% filter(nP_v3 > 1 | nP_v2 > 1), by = c("crop", "DMI_10km", "AC_jord")) %>%
  group_by(crop, DMI_10km, AC_jord, sand, Georegion) %>%
  summarise(cells = n(), nP_v2 = n_distinct(round(drain_effekt_v2, 5)),
            nP_v3 = n_distinct(round(drain_effekt_v3, 5)), .groups = "drop")
inside %>% group_by(crop) %>%
  summarise(subgroups = n(), still_multi_v2 = sum(nP_v2 > 1), still_multi_v3 = sum(nP_v3 > 1)) %>%
  print()
# If splitting by sand/clay class makes P unique again, the parameter set is chosen
# per cell (JB2024), not per profile.
write_csv(inside, file.path(out, "P2_inside_multi_combinations.csv"))


# ---- 3. Lookup table and change per combination ------------------------------------
combo <- d %>%
  group_by(crop, Georegion, DMI_10km, AC_jord) %>%
  summarise(cells = n(),
            share_sand = mean(sand),
            P_v2 = median(drain_effekt_v2), P_v3 = median(drain_effekt_v3),
            .groups = "drop") %>%
  mutate(P_pct = 100 * (P_v3 / P_v2 - 1),
         changed = abs(P_pct) > 0.1)
write_csv(combo, file.path(out, "P3_change_per_grid_x_profile.csv"))


# ---- 4. C vs D: profile, grid, or both? ----------------------------------------------
# log(P_v3 / P_v2) explained by profile only, grid only, both (cell-weighted)
r2 <- combo %>%
  filter(is.finite(P_pct)) %>%
  mutate(y = log(P_v3 / P_v2), AC = factor(AC_jord), DMI = factor(DMI_10km)) %>%
  group_by(crop) %>%
  group_modify(~ tibble(
    R2_profile  = summary(lm(y ~ AC,       data = .x, weights = cells))$r.squared,
    R2_grid     = summary(lm(y ~ DMI,      data = .x, weights = cells))$r.squared,
    R2_both     = summary(lm(y ~ AC + DMI, data = .x, weights = cells))$r.squared,
    R2_sandclay = summary(lm(y ~ I(share_sand > 0.5), data = .x, weights = cells))$r.squared))
r2
write_csv(r2, file.path(out, "P4_variance_explained.csv"))
# R2_profile close to 1  -> whole profiles shifted (source C: parameter set / profile hydraulics)
# R2_grid close to 1     -> the weather grids shifted (source D: weather / Daisy rerun)
# R2_both >> each        -> both; the remainder (1 - R2_both) is profile x grid specific (D)

# Cross-crop check: does the same combination change the same way in both crops?
xc <- combo %>%
  select(crop, Georegion, DMI_10km, AC_jord, cells, P_pct) %>%     # a DMI grid can span two Georegions
  pivot_wider(names_from = crop, values_from = c(P_pct, cells))
cor(xc$P_pct_SB, xc$P_pct_WW, use = "complete.obs")   # ~1 = same cause for both crops


# ---- 5. Per profile: direction, size, uniformity across grids -------------------------
wmed <- \(x, w) { o <- order(x); x <- x[o]; w <- w[o]; x[which(cumsum(w) >= sum(w) / 2)[1]] }   # weighted median
prof <- combo %>%
  group_by(crop, Georegion, AC_jord) %>%
  summarise(grids = n(),
            share_sand = weighted.mean(share_sand, cells),
            pct_cells_changed = 100 * weighted.mean(changed, cells),
            median_pct = wmed(P_pct, cells),
            min_pct = min(P_pct, na.rm = TRUE), max_pct = max(P_pct, na.rm = TRUE),
            sd_pct_across_grids = sqrt(weighted.mean((P_pct - weighted.mean(P_pct, cells, na.rm = TRUE))^2,
                                                     cells, na.rm = TRUE)),
            cells = sum(cells),            # must be last: after this line, cells is a single total
            .groups = "drop") %>%
  mutate(pattern = case_when(pct_cells_changed < 5 ~ "unchanged",
                             sd_pct_across_grids < 2 ~ "uniform shift (whole profile)",
                             TRUE ~ "varies by grid")) %>%
  arrange(crop, desc(abs(median_pct) * cells))
print(prof, n = 60)
write_csv(prof, file.path(out, "P5_change_per_profile.csv"))

prof %>% count(crop, pattern, wt = cells, name = "cells") %>%
  group_by(crop) %>% mutate(pct = 100 * cells / sum(cells)) %>% print()


# ---- 6. Sand / clay ---------------------------------------------------------------
sc <- d %>%
  group_by(crop, Georegion, sand) %>%
  summarise(cells = n(), pct_changed = 100 * mean(changed),
            median_pct = median(P_pct, na.rm = TRUE),
            q05 = quantile(P_pct, 0.05, na.rm = TRUE), q95 = quantile(P_pct, 0.95, na.rm = TRUE),
            .groups = "drop")
print(sc, n = Inf)
write_csv(sc, file.path(out, "P6_change_by_sand_clay.csv"))


# ---- 7. Plots -------------------------------------------------------------------
# 7a. Change per profile, ordered, coloured by pattern
p <- prof %>%
  filter(crop == "SB") %>%
  mutate(AC_jord = fct_reorder(factor(AC_jord), median_pct)) %>%
  ggplot(aes(median_pct, AC_jord, colour = pattern)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_errorbarh(aes(xmin = min_pct, xmax = max_pct), height = 0) +
  geom_point(aes(size = cells)) +
  facet_wrap(~ Georegion, scales = "free_y", nrow = 1, labeller = label_both) +
  labs(x = "Change in P, del v2 -> del v3 (%; point = median, line = range over DMI grids)",
       y = "Soil profile (AC_jord)", colour = NULL, size = "Cells",
       title = "Spring barley: change in P per soil profile") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6))
ggsave(file.path(out, "P7a_SB_change_per_profile.png"), p, width = 14, height = 8, bg = "white")

# 7b. Grid x profile heatmap for the largest profiles (is the change constant along a row?)
top <- prof %>% filter(crop == "SB") %>% slice_max(cells, n = 25) %>% pull(AC_jord)
p <- combo %>%
  filter(crop == "SB", AC_jord %in% top) %>%
  ggplot(aes(factor(DMI_10km), factor(AC_jord), fill = P_pct)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b", midpoint = 0,
                       limits = c(-30, 30), oob = scales::squish, name = "% change") +
  facet_grid(~ Georegion, scales = "free_x", space = "free_x") +
  labs(x = "DMI 10 km grid", y = "Soil profile (AC_jord)",
       title = "Spring barley: change in P per grid x profile (25 largest profiles)",
       subtitle = "A row of one colour = the whole profile shifted; colour varying along a row = weather/Daisy change") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(out, "P7b_SB_grid_x_profile.png"), p, width = 14, height = 7, bg = "white")

# 7c. Cross-crop: same change in SB and WW?
p <- ggplot(xc, aes(P_pct_SB, P_pct_WW)) +
  geom_abline(linetype = "dashed") +
  geom_point(alpha = 0.3, size = 0.8) +
  labs(x = "Change in P, spring barley (%)", y = "Change in P, winter wheat (%)",
       title = "Same grid x profile, both crops")
ggsave(file.path(out, "P7c_SB_vs_WW.png"), p, width = 6, height = 6, bg = "white")

cat("Done. Results in", out, "\n")
