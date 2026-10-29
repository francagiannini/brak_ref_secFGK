# =============================================================================
# Brak / N leaching (NUAR): AU deliverables, version by version
#
#   del v1   (29 May)     GeoDB_20260529   delivered 29.05.2026  first version
#   del v2   (4 June)     GeoDB_20260604   delivered 08.06.2026  second version  -> SEGES Tif_shp_20260604
#   del v3   (28 July)    GeoDB_20260728   delivered 03.08.2026  third (final)    -> SEGES Tif_shp_20260728
#   del corr (20 October) GeoDB_20261020   delivered 26.10.2026  allocation correction (spring barley only)
#
# The leaching field is called NLES5_NUAR_kgNha in the files; it is the NUAR
# N leaching, and all plots and tables call it "N leaching".
# Cells with no soil data (JB2024 = -99, 374,919 cells) are left out of every
# comparison, table and figure. In v1/v2 they hold 0 in all fields, from v3 they are empty.
#
# Sections
#   0. Setup: paths, version table, labels, small helpers
#   1. Inventory: layers, rows, fields per version (schema differences)
#   2. Pick the 10% sample (same gridcodes in every version, stratified by op_id_int)
#   3. Full dataset, one to one, step by step (v1 -> v2 -> v3 -> corr):
#        3a. % cells changed per field, size of change, co-change with N leaching
#        3b. LIST of the fields that changed at each step (+ fields added/removed)
#   4. Summaries per version (full dataset): N leaching and P, by Georegion / JB
#   5. Plots
#        5a. sample: scatter per step, change by JB x Georegion, density by JB
#        5b. sample: density of P and N leaching by Georegion, all versions
#        5c. full dataset: P and N leaching change per step, by Georegion (table + plot)
#        5d. full dataset: MAPS of the change in N leaching and in P, per step
#   6. Side check against SEGES (done 26.10.2026; set run_seges_check <- TRUE to redo)
#   7. Timeline plot
#
# Needs: tidyverse, sf, terra, data.table.
# =============================================================================

library(tidyverse)
library(sf)
library(terra)
library(data.table)   # only for fast reading/writing

brak  <- "C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak"
dlv   <- file.path(brak, "Delivered_AU")
seges <- file.path(brak, "SEGES sammenligning til AU/SEGES sammenligning til AU/Input")
out   <- file.path(brak, "compare_deliverables")
plots <- file.path(out, "plots")
dir.create(plots, recursive = TRUE, showWarnings = FALSE)

run_seges_check <- FALSE    # section 6 already ran; its csv files are reused for the timeline

theme_set(theme_bw(base_size = 10))

# ---- 0. Version table ----------------------------------------------------------
versions <- tribble(
  ~ver,   ~label,                   ~file_date,   ~delivered,   ~note,
  "v1",   "del v1 (29 May)",        "2026-05-29", "2026-05-29", "First version",
  "v2",   "del v2 (4 June)",        "2026-06-04", "2026-06-08", "Second version",
  "v3",   "del v3 (28 July)",       "2026-07-28", "2026-08-03", "Third (final) version",
  "corr", "del corr (20 October)",  "2026-10-20", "2026-10-26", "Allocation correction") %>%
  mutate(file_date = as.Date(file_date), delivered = as.Date(delivered),
         stamp = format(file_date, "%Y%m%d"),
         gdb   = file.path(dlv, paste0("GeoDB_", stamp), paste0("GeoDB_", stamp),
                           paste0("Geodatabase_", stamp, ".gdb")))
file.exists(versions$gdb)                     # all TRUE

lab      <- set_names(versions$label, versions$ver)
ver_cols <- c(v1 = "#a6761d", v2 = "#1b9e77", v3 = "#7570b3", corr = "#d95f02")
names(ver_cols) <- lab[names(ver_cols)]       # colours keyed by label, for the plots

crops <- c(SB = "Spring barley + catch crop", WW = "Winter wheat")
crop_pattern <- c(SB = "Vaarbyg", WW = "VinterHvede")

L_lab <- "N leaching (NUAR, kg N/ha)"
P_lab <- "P (drain_effekt)"

# Kept in memory for every version (full dataset): summaries, maps, SEGES match
slim_vars <- c("NLES5_NUAR_kgNha", "drain_effekt", "N_effekt", "soil_effekt",
               "N_min_fert_foraar_kgNha", "Nhusdyr_kgNHa")
drop_vars <- c("OBJECTID", "Shape_Length", "Shape_Area", "Shape_Leng")
eps    <- 1e-6     # a "change" in the GDB (doubles): anything above floating-point noise
eps_L  <- 0.01     # kg N/ha: N leaching counted as changed above this
no_soil <- -99     # JB2024 code for cells without soil data

# Helper: read the attribute table of one layer, no geometry
read_attr <- function(src, lay, fields = NULL) {
  sql <- paste("SELECT", if (is.null(fields)) "*" else paste(fields, collapse = ", "), "FROM", lay)
  d <- tryCatch(terra::vect(src, layer = lay, query = sql, what = "attributes"),
                error = function(e) st_drop_geometry(st_read(src, query = sql, quiet = TRUE)))
  as_tibble(d) %>% select(-any_of(drop_vars))
}

# Helper: the layer of one crop in one geodatabase
layer_of <- function(gdb, crop) {
  l <- st_layers(gdb)$name
  l <- l[str_detect(l, crop_pattern[[crop]])]
  if (length(l) > 1) stop("More than one ", crop, " layer in ", gdb, ": ", paste(l, collapse = ", "))
  if (length(l) == 0) NA_character_ else l
}

# Helper: compare two full versions of one crop, cell by cell (matched on gridcode),
# cells with soil data only
compare_full <- function(a, b, crop, step) {
  idx  <- match(a$gridcode, b$gridcode)
  keep <- !is.na(idx) & coalesce(a$JB2024 != no_soil, FALSE)
  vars <- setdiff(intersect(names(a), names(b)), "gridcode")
  La <- a$NLES5_NUAR_kgNha[keep]; Lb <- b$NLES5_NUAR_kgNha[idx[keep]]
  L_changed <- abs(Lb - La) > eps_L & !is.na(La) & !is.na(Lb)
  map_dfr(vars, \(v) {
    xa <- a[[v]][keep]; xb <- b[[v]][idx[keep]]
    if (is.numeric(xa) && is.numeric(xb)) {
      d   <- xb - xa
      chg <- (is.na(xa) != is.na(xb)) | (!is.na(d) & abs(d) > eps)
      rel <- !is.na(d) & abs(d) > 0.001 * pmax(abs(xa), 1e-9)
      tibble(variable = v, type = "numeric",
             pct_changed = 100 * mean(chg), pct_changed_rel0.1 = 100 * mean(rel),
             mean_before = mean(xa, na.rm = TRUE), mean_after = mean(xb, na.rm = TRUE),
             mean_diff = mean(d, na.rm = TRUE), median_diff = median(d, na.rm = TRUE),
             q05_diff = quantile(d, 0.05, na.rm = TRUE), q95_diff = quantile(d, 0.95, na.rm = TRUE),
             max_abs_diff = max(abs(d), na.rm = TRUE),
             pct_changed_where_L_changed = 100 * mean(chg[L_changed]))
    } else {
      chg <- coalesce(as.character(xa) != as.character(xb), TRUE)
      tibble(variable = v, type = "text", pct_changed = 100 * mean(chg),
             pct_changed_rel0.1 = 100 * mean(chg),
             pct_changed_where_L_changed = 100 * mean(chg[L_changed]))
    }
  }) %>%
    mutate(crop = crop, step = step, cells_compared = sum(keep),
           pct_cells_L_changed = 100 * mean(L_changed), .before = 1)
}


# ---- 1. Inventory: layers, rows, fields ------------------------------------------
inv <- versions %>%
  select(ver, label, gdb) %>%
  crossing(crop = names(crops)) %>%
  mutate(layer = map2_chr(gdb, crop, layer_of)) %>%
  filter(!is.na(layer)) %>%
  arrange(crop, factor(ver, levels = versions$ver))     # chronological (crossing() sorts alphabetically)
inv %>% select(label, crop, layer)

# Field names per version and crop (first 5 rows only)
fields <- inv %>%
  mutate(names = map2(gdb, layer, ~ names(st_read(.x, query = paste("SELECT * FROM", .y, "LIMIT 5"),
                                                  quiet = TRUE)))) %>%
  select(ver, label, crop, names) %>%
  unnest(names) %>%
  filter(!names %in% c(drop_vars, "Shape", "geom", "geometry"))

schema <- fields %>%
  select(-ver) %>%
  mutate(present = "x") %>%
  pivot_wider(names_from = label, values_from = present, values_fill = "") %>%
  select(crop, names, any_of(versions$label)) %>%
  arrange(crop, names)
print(schema, n = Inf)
write_csv(schema, file.path(out, "01_schema_fields_per_version.csv"))


# ---- 2. Sample: 10% of the cells present in every version, by op_id_int ---------
set.seed(42)
pick <- list()
for (cr in unique(inv$crop)) {
  ids <- inv %>% filter(crop == cr) %>%
    mutate(d = map2(gdb, layer, ~ read_attr(.x, .y, c("gridcode", "op_id_int"))))
  common <- Reduce(intersect, map(ids$d, "gridcode"))
  ref    <- last(ids$d) %>% filter(gridcode %in% common)
  pick[[cr]] <- ref %>%
    group_by(op_id_int) %>%
    filter(row_number() %in% sample(n(), ceiling(0.10 * n()))) %>%
    pull(gridcode)
  cat(cr, ": common cells", length(common), "- sample", length(pick[[cr]]), "\n")
}
rm(ids, ref); gc()


# ---- 3. Full dataset, step by step ----------------------------------------------
cells <- list(); steps <- list(); smp <- list(); slim <- list()
for (cr in unique(inv$crop)) {
  inv_cr <- inv %>% filter(crop == cr)
  prev <- NULL; prev_ver <- NULL
  for (i in seq_len(nrow(inv_cr))) {
    v   <- inv_cr$ver[i]
    cur <- read_attr(inv_cr$gdb[i], inv_cr$layer[i])
    cat(cr, lab[[v]], ":", nrow(cur), "rows\n")

    cells[[paste(cr, v)]] <- tibble(crop = cr, version = lab[[v]], rows = nrow(cur),
                                    gridcodes = n_distinct(cur$gridcode),
                                    duplicated = sum(duplicated(cur$gridcode)),
                                    no_soil_cells = sum(cur$JB2024 == no_soil, na.rm = TRUE))
    smp[[paste(cr, v)]]  <- cur %>% filter(gridcode %in% pick[[cr]]) %>%
      mutate(crop = cr, ver = v, version = lab[[v]])
    slim[[cr]][[v]] <- cur %>% select(gridcode, Georegion, JB2024, any_of(slim_vars))

    if (!is.null(prev)) {
      step <- paste(lab[[prev_ver]], "->", lab[[v]])
      cells[[paste(cr, v)]] <- cells[[paste(cr, v)]] %>%
        mutate(step = step,
               only_in_previous = sum(!prev$gridcode %in% cur$gridcode),
               only_in_this     = sum(!cur$gridcode %in% prev$gridcode))
      steps[[paste(cr, v)]] <- compare_full(prev, cur, cr, step)
    }
    prev <- cur; prev_ver <- v
    rm(cur); gc()
  }
  rm(prev); gc()
}

cells <- bind_rows(cells)
cells
write_csv(cells, file.path(out, "03_cells_per_version_added_removed.csv"))

step_levels <- unique(bind_rows(steps)$step)                      # chronological
steps <- bind_rows(steps)
write_csv(steps, file.path(out, "03_step_changes_full_dataset.csv"))

# 3a. Short table: fields that changed by more than 0.1% (relative) in at least 0.1% of cells
steps_short <- steps %>%
  filter(pct_changed_rel0.1 >= 0.1) %>%
  select(crop, step, variable, pct_changed_rel0.1, mean_before, mean_after, median_diff,
         q05_diff, q95_diff, max_abs_diff, pct_changed_where_L_changed) %>%
  arrange(crop, factor(step, levels = step_levels), desc(pct_changed_rel0.1))
print(steps_short, n = Inf)
write_csv(steps_short, file.path(out, "03_step_changes_short.csv"))

# 3b. One line per step: the list of fields that changed, plus fields added / removed
fields_by_ver <- split(fields$names, paste(fields$crop, fields$ver))
changed_list <- steps_short %>%
  group_by(crop, step) %>%
  summarise(n_fields_changed = n(),
            fields_changed = paste0(variable, " (", sprintf("%.1f", pct_changed_rel0.1), "%)",
                                    collapse = ", "),
            .groups = "drop") %>%
  right_join(steps %>% distinct(crop, step, cells_compared, pct_cells_L_changed), by = c("crop", "step")) %>%
  mutate(n_fields_changed = coalesce(n_fields_changed, 0L),
         from = versions$ver[match(str_remove(step, " -> .*"), versions$label)],
         to   = versions$ver[match(str_remove(step, ".* -> "), versions$label)],
         fields_removed = map2_chr(paste(crop, from), paste(crop, to),
                                   ~ paste(setdiff(fields_by_ver[[.x]], fields_by_ver[[.y]]), collapse = ", ")),
         fields_added   = map2_chr(paste(crop, from), paste(crop, to),
                                   ~ paste(setdiff(fields_by_ver[[.y]], fields_by_ver[[.x]]), collapse = ", "))) %>%
  select(crop, step, cells_compared, pct_cells_L_changed, n_fields_changed,
         fields_changed, fields_added, fields_removed) %>%
  arrange(crop, factor(step, levels = step_levels))
print(changed_list, width = Inf)
write_csv(changed_list, file.path(out, "03_fields_changed_per_step.csv"))

# Heatmap: % of cells changed (> 0.1% relative), field x step
p <- steps %>%
  mutate(crop_name = crops[crop],
         step = factor(step, levels = step_levels),
         variable = fct_reorder(variable, pct_changed_rel0.1, .fun = max)) %>%
  ggplot(aes(step, variable, fill = pct_changed_rel0.1)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(pct_changed_rel0.1 >= 0.05, sprintf("%.1f", pct_changed_rel0.1), "")), size = 2.6) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0, 100), name = "% cells\nchanged") +
  facet_wrap(~ crop_name, scales = "free_x") +
  labs(x = NULL, y = NULL,
       title = "Share of cells where each field changed (> 0.1%), version to version",
       subtitle = "Full dataset, cells with soil data (JB -99 excluded)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(plots, "03_heatmap_pct_cells_changed.png"), p, width = 10, height = 7)


# ---- 4. Summaries per version (full dataset, JB -99 excluded) ---------------------
soil <- \(d) d %>% filter(JB2024 != no_soil)

summ <- map_dfr(names(slim), \(cr) map_dfr(names(slim[[cr]]), \(v)
  soil(slim[[cr]][[v]]) %>%
    summarise(cells = n(),
              L_mean = mean(NLES5_NUAR_kgNha, na.rm = TRUE),
              L_median = median(NLES5_NUAR_kgNha, na.rm = TRUE),
              L_q05 = quantile(NLES5_NUAR_kgNha, 0.05, na.rm = TRUE),
              L_q95 = quantile(NLES5_NUAR_kgNha, 0.95, na.rm = TRUE),
              P_mean = mean(drain_effekt, na.rm = TRUE),
              P_median = median(drain_effekt, na.rm = TRUE)) %>%
    mutate(crop = cr, version = lab[[v]], .before = 1)))
summ
write_csv(summ, file.path(out, "04_summary_per_version.csv"))

summ_geo <- map_dfr(names(slim), \(cr) map_dfr(names(slim[[cr]]), \(v)
  soil(slim[[cr]][[v]]) %>%
    group_by(Georegion) %>%
    summarise(cells = n(), L_mean = mean(NLES5_NUAR_kgNha, na.rm = TRUE),
              L_median = median(NLES5_NUAR_kgNha, na.rm = TRUE),
              P_mean = mean(drain_effekt, na.rm = TRUE), .groups = "drop") %>%
    mutate(crop = cr, version = lab[[v]], .before = 1)))
write_csv(summ_geo, file.path(out, "04_summary_per_version_Georegion.csv"))

summ_jb <- map_dfr(names(slim), \(cr) map_dfr(names(slim[[cr]]), \(v)
  soil(slim[[cr]][[v]]) %>%
    group_by(JB2024) %>%
    summarise(cells = n(), L_mean = mean(NLES5_NUAR_kgNha, na.rm = TRUE),
              L_median = median(NLES5_NUAR_kgNha, na.rm = TRUE),
              P_mean = mean(drain_effekt, na.rm = TRUE), .groups = "drop") %>%
    mutate(crop = cr, version = lab[[v]], .before = 1)))
write_csv(summ_jb, file.path(out, "04_summary_per_version_JB.csv"))

p <- summ %>%
  left_join(versions %>% select(version = label, file_date), by = "version") %>%
  mutate(crop_name = crops[crop]) %>%
  ggplot(aes(file_date, L_median, colour = crop_name)) +
  geom_pointrange(aes(ymin = L_q05, ymax = L_q95), position = position_dodge(width = 4)) +
  geom_line(position = position_dodge(width = 4)) +
  scale_x_date(breaks = versions$file_date, labels = versions$label) +
  labs(x = NULL, y = paste(L_lab, "- median, 5-95%"), colour = NULL,
       title = "National N leaching per deliverable (full dataset, JB -99 excluded)") +
  theme(legend.position = "bottom")
ggsave(file.path(plots, "04_Nleaching_per_version.png"), p, width = 9, height = 5)


# ---- 5. Plots ----------------------------------------------------------------------
# 5a. On the 10% sample (JB -99 excluded from the figures)
smp <- bind_rows(smp)
fwrite(smp, file.path(out, "sample10pct_deliverables.csv"), sep = ";")
smp <- smp %>%
  filter(JB2024 != no_soil) %>%
  mutate(version = factor(version, levels = versions$label),
         JB2024 = factor(JB2024), Georegion = factor(Georegion), crop_name = crops[crop])

# Before and after for the same cell, per step
step_pairs <- smp %>%
  select(crop, gridcode, ver, version, JB2024, Georegion, NLES5_NUAR_kgNha, drain_effekt) %>%
  arrange(crop, gridcode, version) %>%
  group_by(crop, gridcode) %>%
  mutate(L_before = lag(NLES5_NUAR_kgNha), P_before = lag(drain_effekt),
         step = paste(lag(as.character(version)), "->", version)) %>%
  ungroup() %>%
  filter(!is.na(L_before)) %>%
  mutate(step = factor(step, levels = step_levels),
         dL = NLES5_NUAR_kgNha - L_before,
         dP_pct = 100 * (drain_effekt / P_before - 1))

for (cr in unique(step_pairs$crop)) {
  d  <- step_pairs %>% filter(crop == cr)
  ds <- d %>% slice_sample(n = min(nrow(d), 150000))     # lighter scatter

  p <- ggplot(ds, aes(L_before, NLES5_NUAR_kgNha)) +
    geom_abline(linetype = "dashed", colour = "grey40") +
    geom_point(alpha = 0.08, size = 0.3) +
    facet_wrap(~ step) + coord_equal() +
    labs(x = "N leaching, earlier version (kg N/ha)", y = "N leaching, later version (kg N/ha)",
         title = paste(crops[[cr]], "- N leaching, version to version (sample)"))
  ggsave(file.path(plots, paste0("05a_", cr, "_scatter_Nleaching_per_step.png")), p, width = 11, height = 4.5)

  p <- ggplot(ds, aes(P_before, drain_effekt)) +
    geom_abline(linetype = "dashed", colour = "grey40") +
    geom_point(alpha = 0.08, size = 0.3) +
    facet_wrap(~ step) + coord_equal() +
    labs(x = "P (drain_effekt), earlier version", y = "P (drain_effekt), later version",
         title = paste(crops[[cr]], "- percolation term P, version to version (sample)"))
  ggsave(file.path(plots, paste0("05a_", cr, "_scatter_P_per_step.png")), p, width = 11, height = 4.5)

  p <- ggplot(d, aes(JB2024, dL)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
    facet_grid(step ~ Georegion, labeller = labeller(Georegion = label_both)) +
    labs(x = "JB2024", y = "Change in N leaching (kg N/ha)",
         title = paste(crops[[cr]], "- change in N leaching per step, by JB and Georegion (sample)"))
  ggsave(file.path(plots, paste0("05a_", cr, "_change_Nleaching_JB_Georegion.png")), p, width = 13, height = 8)
}

p <- ggplot(smp, aes(NLES5_NUAR_kgNha, colour = version)) +
  geom_density() +
  scale_colour_manual(values = ver_cols) +
  facet_grid(crop_name ~ JB2024, scales = "free_y") +
  coord_cartesian(xlim = c(0, 150)) +
  labs(x = L_lab, colour = NULL,
       title = "Distribution of N leaching per deliverable, by JB (sample)") +
  theme(legend.position = "bottom")
ggsave(file.path(plots, "05a_density_Nleaching_by_JB.png"), p, width = 14, height = 6)

# 5b. Density of P and N leaching by Georegion, all versions (as in the earlier analysis)
for (cr in unique(smp$crop)) {
  d <- smp %>% filter(crop == cr)

  p <- ggplot(d, aes(drain_effekt, colour = version)) +
    geom_density(linewidth = 0.7, adjust = 0.7) +
    facet_wrap(~ Georegion, nrow = 1, labeller = label_both) +
    scale_colour_manual(values = ver_cols) +
    labs(x = P_lab, y = "Density", colour = NULL,
         title = paste(crops[[cr]], "- distribution of P by Georegion (sample)")) +
    theme(legend.position = "bottom")
  ggsave(file.path(plots, paste0("05b_", cr, "_density_P_by_Georegion.png")), p, width = 13, height = 3.8)

  p <- ggplot(d, aes(NLES5_NUAR_kgNha, colour = version)) +
    geom_density(linewidth = 0.7, adjust = 0.7) +
    facet_wrap(~ Georegion, nrow = 1, labeller = label_both) +
    scale_colour_manual(values = ver_cols) +
    coord_cartesian(xlim = c(0, quantile(d$NLES5_NUAR_kgNha, 0.995, na.rm = TRUE))) +
    labs(x = L_lab, y = "Density", colour = NULL,
         title = paste(crops[[cr]], "- distribution of N leaching by Georegion (sample)")) +
    theme(legend.position = "bottom")
  ggsave(file.path(plots, paste0("05b_", cr, "_density_Nleaching_by_Georegion.png")), p, width = 13, height = 3.8)
}

# 5c. P and N leaching change per step, by Georegion (full dataset)
pair_full <- \(cr, from, to) {
  a <- soil(slim[[cr]][[from]]); b <- slim[[cr]][[to]]
  idx <- match(a$gridcode, b$gridcode)
  tibble(gridcode = a$gridcode, Georegion = a$Georegion, JB2024 = a$JB2024,
         P_before = a$drain_effekt, P_after = b$drain_effekt[idx],
         L_before = a$NLES5_NUAR_kgNha, L_after = b$NLES5_NUAR_kgNha[idx]) %>%
    filter(!is.na(idx)) %>%
    mutate(dP = P_after - P_before,
           dP_pct = ifelse(P_before > 0, 100 * (P_after / P_before - 1), NA_real_),
           dL = L_after - L_before,
           dL_pct = ifelse(L_before > 0, 100 * (L_after / L_before - 1), NA_real_))
}
summ_change <- \(d) d %>%
  summarise(cells = n(),
            P_pct_cells_changed = 100 * mean(abs(dP_pct) > 0.1, na.rm = TRUE),
            P_median_pct = median(dP_pct, na.rm = TRUE),
            P_q05_pct = quantile(dP_pct, 0.05, na.rm = TRUE),
            P_q95_pct = quantile(dP_pct, 0.95, na.rm = TRUE),
            L_pct_cells_changed = 100 * mean(abs(dL) > eps_L, na.rm = TRUE),
            L_mean_before = mean(L_before, na.rm = TRUE),
            L_mean_after = mean(L_after, na.rm = TRUE),
            L_median_diff = median(dL, na.rm = TRUE),
            L_q05_diff = quantile(dL, 0.05, na.rm = TRUE),
            L_q95_diff = quantile(dL, 0.95, na.rm = TRUE),
            L_median_pct = median(dL_pct, na.rm = TRUE),
            .groups = "drop")

pl_geo <- list()
for (cr in names(slim)) {
  vv <- names(slim[[cr]])                       # chronological (filled in that order in section 3)
  for (i in seq_along(vv)[-1]) {
    step <- paste(lab[[vv[i - 1]]], "->", lab[[vv[i]]])
    d <- pair_full(cr, vv[i - 1], vv[i])
    pl_geo[[paste(cr, step)]] <- bind_rows(
      d %>% mutate(Georegion = as.character(Georegion)) %>% group_by(Georegion) %>% summ_change(),
      d %>% summ_change() %>% mutate(Georegion = "All")) %>%
      mutate(crop = cr, step = step, .before = 1)
    rm(d); gc()
  }
}
pl_geo <- bind_rows(pl_geo) %>%
  mutate(step = factor(step, levels = step_levels)) %>%
  arrange(crop, step, Georegion)
print(pl_geo, n = Inf, width = Inf)
write_csv(pl_geo, file.path(out, "05c_P_and_Nleaching_change_by_Georegion.csv"))

# Median change and 5-95% range per Georegion and step
p <- pl_geo %>%
  filter(Georegion != "All") %>%
  select(crop, step, Georegion, P_median_pct, P_q05_pct, P_q95_pct, L_median_pct = L_median_diff,
         L_q05_pct = L_q05_diff, L_q95_pct = L_q95_diff) %>%
  pivot_longer(-c(crop, step, Georegion), names_to = c("term", ".value"), names_pattern = "(P|L)_(.*)_pct") %>%
  mutate(term = recode(term, P = "P (drain_effekt), % change", L = "N leaching, change in kg N/ha"),
         crop_name = crops[crop]) %>%
  ggplot(aes(Georegion, median, ymin = q05, ymax = q95, colour = step)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.3) +
  facet_grid(term ~ crop_name, scales = "free_y") +
  labs(x = "Georegion", y = "Median and 5-95% range", colour = NULL,
       title = "Change in P and in N leaching per step, by Georegion (full dataset)") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2))
ggsave(file.path(plots, "05c_change_P_Nleaching_by_Georegion.png"), p, width = 11, height = 7)

# 5d. Maps of the change per step (full dataset, 500 m cells = mean of 25 x 100 m cells)
# Cell centres come from the geodatabase geometry; read once in chunks and cached.
xy_file <- file.path(out, "gridcode_xy.csv")
if (file.exists(xy_file)) {
  xy <- fread(xy_file) %>% as_tibble()
} else {
  gi <- inv %>% filter(crop == "SB", ver == "v3")
  n_all <- nrow(slim$SB$v3); chunk <- 400000; start <- 0; got <- 0; xy <- list()
  while (got < n_all) {
    x <- st_read(gi$gdb, quiet = TRUE,
                 query = sprintf("SELECT gridcode FROM %s WHERE OBJECTID > %d AND OBJECTID <= %d",
                                 gi$layer, start, start + chunk))
    if (nrow(x) > 0) {
      cc <- st_coordinates(suppressWarnings(st_centroid(st_geometry(x))))
      xy[[length(xy) + 1]] <- tibble(gridcode = x$gridcode, x = cc[, 1], y = cc[, 2])
      got <- got + nrow(x)
    }
    start <- start + chunk
    cat("coordinates:", got, "of", n_all, "\n")
    if (start > 10 * n_all) stop("OBJECTID range not as expected")
    rm(x); gc()
  }
  xy <- bind_rows(xy)
  fwrite(xy, xy_file)
}
crs_dk <- st_crs(st_read(inv$gdb[1], query = paste("SELECT * FROM", inv$layer[1], "LIMIT 1"), quiet = TRUE))$wkt
tmpl <- rast(ext(min(xy$x) - 250, max(xy$x) + 250, min(xy$y) - 250, max(xy$y) + 250),
             resolution = 500, crs = crs_dk)

to_grid <- \(d, col) {                     # mean of a column on the 500 m grid, as a data frame
  i <- match(d$gridcode, xy$gridcode)
  ok <- !is.na(i) & !is.na(d[[col]]) & is.finite(d[[col]])
  r <- rasterize(cbind(xy$x[i[ok]], xy$y[i[ok]]), tmpl, values = d[[col]][ok], fun = "mean")
  as.data.frame(r, xy = TRUE, na.rm = TRUE) %>% setNames(c("x", "y", "value"))
}

for (cr in names(slim)) {
  vv <- names(slim[[cr]])
  mL <- list(); mP <- list()
  for (i in seq_along(vv)[-1]) {
    step <- paste(lab[[vv[i - 1]]], "->", lab[[vv[i]]])
    d <- pair_full(cr, vv[i - 1], vv[i])
    mL[[step]] <- to_grid(d, "dL") %>% mutate(step = step)
    mP[[step]] <- to_grid(d, "dP_pct") %>% mutate(step = step)
    rm(d); gc()
  }
  mL <- bind_rows(mL) %>% mutate(step = factor(step, levels = step_levels))
  mP <- bind_rows(mP) %>% mutate(step = factor(step, levels = step_levels))
  limL <- quantile(abs(mL$value), 0.99, na.rm = TRUE)
  limP <- quantile(abs(mP$value), 0.99, na.rm = TRUE)

  p <- ggplot(mL, aes(x, y, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(low = "#2166ac", mid = "grey92", high = "#b2182b", midpoint = 0,
                         limits = c(-limL, limL), oob = scales::squish, name = "kg N/ha") +
    facet_wrap(~ step, nrow = 1) + coord_equal() + theme_void() +
    labs(title = paste(crops[[cr]], "- change in N leaching per step"),
         subtitle = "Mean change per 500 m cell; red = more leaching, blue = less. JB -99 excluded")
  ggsave(file.path(plots, paste0("05d_", cr, "_map_change_Nleaching.png")), p,
         width = 5 * length(levels(droplevels(mL$step))), height = 5.5, dpi = 200)

  p <- ggplot(mP, aes(x, y, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(low = "#2166ac", mid = "grey92", high = "#b2182b", midpoint = 0,
                         limits = c(-limP, limP), oob = scales::squish, name = "% change") +
    facet_wrap(~ step, nrow = 1) + coord_equal() + theme_void() +
    labs(title = paste(crops[[cr]], "- change in P (drain_effekt) per step"),
         subtitle = "Mean % change per 500 m cell; red = higher P, blue = lower. JB -99 excluded")
  ggsave(file.path(plots, paste0("05d_", cr, "_map_change_P.png")), p,
         width = 5 * length(levels(droplevels(mP$step))), height = 5.5, dpi = 200)
}

# ---- 6. Side check against SEGES ---------------------------------------------------
# Ran on 26.10.2026: both SEGES packages are identical to del v2 / del v3.
# Set run_seges_check <- TRUE at the top to run it again.
if (run_seges_check) {
seges_files <- tribble(
  ~seges,          ~crop, ~dir,                                  ~base,
  "SEGES 4 June",  "SB",  "Tif_shp_20260604/Tif_shp_20260604",   "DK_Grid_VaarbygEfterafg",
  "SEGES 4 June",  "WW",  "Tif_shp_20260604/Tif_shp_20260604",   "DK_Grid_VinterHvede",
  "SEGES 28 July", "SB",  "Tif_shp_20260728/Tif_shp_20260728",   "DK_Grid_VaarbygEfterafg_20260728",
  "SEGES 28 July", "WW",  "Tif_shp_20260728/Tif_shp_20260728",   "DK_Grid_VinterHvede_20260728") %>%
  mutate(tif  = file.path(.env$seges, dir, paste0(base, ".tif")),
         shp1 = file.path(.env$seges, dir, paste0(base, "_1_af_2.shp")),
         shp2 = file.path(.env$seges, dir, paste0(base, "_2_af_2.shp")),
         expected = c("v2", "v2", "v3", "v3"))
file.exists(c(seges_files$tif, seges_files$shp1, seges_files$shp2))

# Shapefile field names are cut to 10 characters: match them to the GDB names
key10 <- \(x) toupper(substr(x, 1, 10))

shp <- list()
for (i in seq_len(nrow(seges_files))) {
  f <- seges_files[i, ]
  shp[[paste(f$seges, f$crop)]] <- bind_rows(
    read_attr(f$shp1, tools::file_path_sans_ext(basename(f$shp1))),
    read_attr(f$shp2, tools::file_path_sans_ext(basename(f$shp2))))
  cat(f$seges, f$crop, ":", nrow(shp[[paste(f$seges, f$crop)]]), "rows\n")
}

# 6a. Match matrix: every SEGES map against every deliverable, on the key fields
match_tbl <- list()
for (i in seq_len(nrow(seges_files))) {
  f <- seges_files[i, ]
  s <- shp[[paste(f$seges, f$crop)]]
  for (v in names(slim[[f$crop]])) {
    g    <- slim[[f$crop]][[v]]
    idx  <- match(g$gridcode, s$gridcode)
    keep <- !is.na(idx)
    for (gv in setdiff(names(g), "gridcode")) {
      sv <- names(s)[match(key10(gv), key10(names(s)))]
      if (is.na(sv) || !is.numeric(s[[sv]])) next
      d <- s[[sv]][idx[keep]] - g[[gv]][keep]
      match_tbl[[length(match_tbl) + 1]] <- tibble(
        seges = f$seges, crop = f$crop, deliverable = lab[[v]], expected = v == f$expected,
        field_gdb = gv, field_shp = sv,
        cells_gdb = nrow(g), cells_shp = nrow(s), cells_matched = sum(keep),
        mean_abs_diff = mean(abs(d), na.rm = TRUE), max_abs_diff = max(abs(d), na.rm = TRUE),
        pct_equal_1e3 = 100 * mean(abs(d) <= 0.0006, na.rm = TRUE))   # shapefile keeps ~3 decimals
    }
  }
}
match_tbl <- bind_rows(match_tbl)
match_tbl %>%
  filter(field_gdb %in% c("NLES5_NUAR_kgNha", "drain_effekt")) %>%
  select(seges, crop, deliverable, field_gdb, pct_equal_1e3, mean_abs_diff) %>%
  arrange(seges, crop, field_gdb, desc(pct_equal_1e3)) %>%
  print(n = Inf)
write_csv(match_tbl, file.path(out, "06a_match_SEGES_vs_each_deliverable.csv"))

p <- match_tbl %>%
  mutate(deliverable = factor(deliverable, levels = versions$label), crop_name = crops[crop]) %>%
  ggplot(aes(deliverable, field_gdb, fill = pct_equal_1e3)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f", pct_equal_1e3)), size = 2.6) +
  scale_fill_distiller(palette = "Greens", direction = 1, limits = c(0, 100), name = "% cells\nequal") +
  facet_grid(crop_name ~ seges) +
  labs(x = NULL, y = NULL, title = "Which deliverable is behind each SEGES map? (% of cells equal)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(plots, "06a_match_SEGES_vs_deliverables.png"), p, width = 11, height = 7)

# 6b. Expected pairs, every field one to one (full dataset)
pair_tbl <- list()
for (i in seq_len(nrow(seges_files))) {
  f   <- seges_files[i, ]
  gi  <- inv %>% filter(ver == f$expected, crop == f$crop)
  g   <- read_attr(gi$gdb, gi$layer)
  s   <- shp[[paste(f$seges, f$crop)]]
  idx <- match(g$gridcode, s$gridcode); keep <- !is.na(idx)
  map_shp <- tibble(field_gdb = setdiff(names(g), "gridcode")) %>%
    mutate(field_shp = names(s)[match(key10(field_gdb), key10(names(s)))])
  pair_tbl[[i]] <- map_shp %>%
    mutate(seges = f$seges, crop = f$crop, deliverable = lab[[f$expected]],
           cells_gdb = nrow(g), cells_shp = nrow(s), cells_matched = sum(keep),
           cells_only_gdb = sum(!keep), cells_only_shp = sum(!s$gridcode %in% g$gridcode),
           .before = 1) %>%
    mutate(res = map2(field_gdb, field_shp, \(gv, sv) {
      if (is.na(sv)) return(tibble(status = "not in shapefile"))
      a <- g[[gv]][keep]; b <- s[[sv]][idx[keep]]
      if (is.numeric(a) && is.numeric(b)) {
        d <- b - a
        tibble(status = "compared", mean_abs_diff = mean(abs(d), na.rm = TRUE),
               max_abs_diff = max(abs(d), na.rm = TRUE),
               pct_equal_1e3 = 100 * mean(abs(d) <= 0.0006, na.rm = TRUE),
               pct_equal_rel = 100 * mean(abs(d) <= 1e-4 * pmax(abs(a), 1), na.rm = TRUE))
      } else {
        tibble(status = "compared", pct_equal_1e3 = 100 * mean(as.character(a) == as.character(b), na.rm = TRUE))
      }
    })) %>%
    unnest(res)
  rm(g); gc()
}
pair_tbl <- bind_rows(pair_tbl)
print(pair_tbl %>% select(seges, crop, deliverable, field_gdb, field_shp, status, pct_equal_1e3, max_abs_diff), n = Inf)
write_csv(pair_tbl, file.path(out, "06b_expected_pairs_every_field.csv"))
# A match: cells_only_gdb = cells_only_shp = 0 and pct_equal ~ 100 for every field
# (shapefile values are rounded, so max_abs_diff up to ~0.0005 is still equal)

# 6c. GeoTIFF vs GDB on sampled cell centres (1 in 200 cells, with geometry)
tif_tbl <- list()
for (i in seq_len(nrow(seges_files))) {
  f  <- seges_files[i, ]
  r  <- rast(f$tif)
  print(r)
  for (v in names(slim[[f$crop]])) {
    gi  <- inv %>% filter(ver == v, crop == f$crop)
    pts <- st_read(gi$gdb, query = sprintf("SELECT * FROM %s WHERE OBJECTID %% 200 = 0", gi$layer),
                   quiet = TRUE)
    pts <- st_transform(pts, crs(r))
    val <- terra::extract(r, vect(st_centroid(st_geometry(pts))))[, 2]
    a   <- st_drop_geometry(pts) %>% select(where(is.numeric), -any_of(c(drop_vars, "gridcode")))
    tif_tbl[[length(tif_tbl) + 1]] <- tibble(field_gdb = names(a)) %>%
      mutate(seges = f$seges, crop = f$crop, deliverable = lab[[v]], expected = v == f$expected,
             cells = sum(!is.na(val)),
             mean_abs_diff = map_dbl(field_gdb, ~ mean(abs(val - a[[.x]]), na.rm = TRUE)),
             pct_equal_0.01 = map_dbl(field_gdb, ~ 100 * mean(abs(val - a[[.x]]) <= 0.01, na.rm = TRUE)),
             .before = 1)
  }
}
tif_tbl <- bind_rows(tif_tbl)
# Which GDB field does each TIF hold, and from which deliverable?
tif_best <- tif_tbl %>%
  group_by(seges, crop, deliverable) %>%
  slice_min(mean_abs_diff, n = 1) %>%
  ungroup() %>%
  arrange(seges, crop, mean_abs_diff)
print(tif_best, n = Inf)
write_csv(tif_tbl,  file.path(out, "06c_TIF_vs_deliverables_all_fields.csv"))
write_csv(tif_best, file.path(out, "06c_TIF_vs_deliverables_best_field.csv"))


}   # end of section 6

# ---- 7. Timeline plot ------------------------------------------------------------
# Result of the side check, to print on the arrows (N leaching, % cells equal)
chk <- read_csv(file.path(out, "06a_match_SEGES_vs_each_deliverable.csv"), show_col_types = FALSE) %>%
  filter(expected, field_gdb == "NLES5_NUAR_kgNha") %>%
  group_by(seges) %>%
  summarise(txt = paste0(crop, " ", sprintf("%.1f", pct_equal_1e3), "% equal", collapse = "\n"))

tl_au <- versions %>% mutate(row = "AU deliverables")
en_day <- \(d) paste(as.integer(format(d, "%d")), month.abb[as.integer(format(d, "%m"))])   # English, whatever the locale
tl_se <- tibble(label = c("SEGES Tif_shp 4 June", "SEGES Tif_shp 28 July"),
                file_date = as.Date(c("2026-06-04", "2026-07-28")),
                from = as.Date(c("2026-06-04", "2026-07-28")),
                seges = c("SEGES 4 June", "SEGES 28 July"),
                row = "SEGES maps") %>%
  left_join(chk, by = "seges")

p <- ggplot() +
  geom_line(data = tl_au, aes(file_date, row, group = row), colour = "grey60") +
  geom_segment(data = tl_se, aes(x = from, xend = file_date, y = "AU deliverables", yend = "SEGES maps"),
               arrow = arrow(length = unit(2, "mm")), colour = "grey40") +
  geom_point(data = tl_au, aes(file_date, row, colour = label), size = 4) +
  geom_text(data = tl_au %>% filter(ver != "v2"),
            aes(file_date, row, label = paste0(label, "\ndelivered ", en_day(delivered))),
            vjust = -1, size = 3) +
  geom_text(data = tl_au %>% filter(ver == "v2"),            # v1 and v2 are 6 days apart: v2 below the line
            aes(file_date, row, label = paste0(label, "\ndelivered ", en_day(delivered))),
            vjust = 1.8, hjust = -0.05, size = 3) +
  geom_point(data = tl_se, aes(file_date, row), shape = 15, size = 4) +
  geom_text(data = tl_se, aes(file_date, row, label = paste0(label, "\n", coalesce(txt, ""))),
            vjust = 1.6, size = 3) +
  scale_colour_manual(values = ver_cols, guide = "none") +
  scale_x_date(breaks = seq(as.Date("2026-06-01"), as.Date("2026-11-01"), by = "month"),
               labels = \(d) month.abb[as.integer(format(d, "%m"))],
               expand = expansion(mult = 0.08)) +
  scale_y_discrete(limits = c("SEGES maps", "AU deliverables"), expand = expansion(add = 0.8)) +
  labs(x = NULL, y = NULL, title = "AU deliverables and the SEGES maps built on them") +
  theme(panel.grid.major.y = element_blank())
ggsave(file.path(plots, "07_timeline.png"), p, width = 11, height = 4)

cat("\nDone. Results in", out, "\n")
