# =============================================================================
# Brak / NUAR: which P-value delivery (Note A, "P-værdier på kort") is behind
# the P term (drain_effekt) in each N leaching deliverable?
#
# SGAV's understanding (request of 28.10.2026):
#   P delivery 29 May 2026  -> leaching data delivered 8 June  (del v2)
#   P delivery 3 Aug 2026   -> leaching data delivered 3 Aug   (del v3)
# This script tests that, for both crops, on the full dataset.
#
# Sections
#   0. Setup: put the P-value geodatabases in Brak/Delivered_AU_P (any subfolders)
#   1. Read the P deliveries (gridcode + all P_* fields)
#   2. Read drain_effekt from the leaching deliverables (v1, v2, v3, corr)
#   3. Match matrix: every leaching version x every P delivery x every P field
#   4. Best match per leaching version, and where the remaining cells differ
#   5. P deliveries compared with each other (what changed between them)
#   6. Plot
# Output: Brak/compare_deliverables/P_deliveries
# =============================================================================

library(tidyverse)
library(sf)
library(terra)

brak  <- "C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak"
P_dir <- file.path(brak, "Delivered_AU_P")      # <- folder with the P-value geodatabases
dlv   <- file.path(brak, "Delivered_AU")
out   <- file.path(brak, "compare_deliverables", "P_deliveries")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 10))

tol <- 0.0006          # P is delivered with 3 decimals: equal within rounding
no_soil <- -99

months <- c("January", "February", "March", "April", "May", "June", "July",
            "August", "September", "October", "November", "December")
nice_date <- \(d) paste(as.integer(format(d, "%d")), months[as.integer(format(d, "%m"))])

# Read the attribute table of one layer (no geometry)
read_attr <- function(src, lay, fields = NULL) {
  sql <- paste("SELECT", if (is.null(fields)) "*" else paste(fields, collapse = ", "), "FROM", lay)
  d <- tryCatch(terra::vect(src, layer = lay, query = sql, what = "attributes"),
                error = function(e) st_drop_geometry(st_read(src, query = sql, quiet = TRUE)))
  as_tibble(d)
}


# ---- 1. P deliveries ----------------------------------------------------------
P_gdbs <- list.dirs(P_dir, recursive = TRUE, full.names = TRUE)
P_gdbs <- P_gdbs[str_detect(P_gdbs, "\\.gdb$")]
P_gdbs                                            # check: one per P delivery

P_tbl <- tibble(gdb = P_gdbs) %>%
  mutate(date  = as.Date(str_extract(basename(gdb), "\\d{8}"), "%Y%m%d"),
         label = paste("P del", nice_date(date))) %>%
  arrange(date)
P_tbl                                             # if a date is NA, set it by hand

P <- list()
for (i in seq_len(nrow(P_tbl))) {
  lays <- st_layers(P_tbl$gdb[i])
  print(lays)
  lay  <- lays$name[which.max(lays$features)]     # the big table
  hdr  <- names(st_read(P_tbl$gdb[i], query = paste("SELECT * FROM", lay, "LIMIT 1"), quiet = TRUE))
  id   <- hdr[tolower(hdr) == "gridcode"]
  pf   <- hdr[str_detect(hdr, regex("^p_", ignore_case = TRUE))]
  cat(P_tbl$label[i], ": layer", lay, "- id", id, "- P fields:", paste(pf, collapse = ", "), "\n")
  P[[P_tbl$label[i]]] <- read_attr(P_tbl$gdb[i], lay, c(id, pf)) %>%
    rename(gridcode = all_of(id))
}
sapply(P, nrow)
sapply(P, \(d) sum(duplicated(d$gridcode)))       # should be 0


# ---- 2. drain_effekt from the leaching deliverables --------------------------
versions <- tribble(
  ~ver,   ~label,                   ~stamp,
  "v1",   "del v1 (29 May)",        "20260529",
  "v2",   "del v2 (4 June)",        "20260604",
  "v3",   "del v3 (28 July)",       "20260728",
  "corr", "del corr (20 October)",  "20261020") %>%
  mutate(gdb = file.path(dlv, paste0("GeoDB_", stamp), paste0("GeoDB_", stamp),
                         paste0("Geodatabase_", stamp, ".gdb")))
crops <- c(SB = "Vaarbyg", WW = "VinterHvede")

L <- list()
for (i in seq_len(nrow(versions))) for (cr in names(crops)) {
  l <- st_layers(versions$gdb[i])$name
  l <- l[str_detect(l, crops[[cr]])]
  if (length(l) != 1) next
  L[[paste(cr, versions$ver[i])]] <- read_attr(versions$gdb[i], l,
                                               c("gridcode", "JB2024", "Georegion", "AC_jord", "drain_effekt")) %>%
    filter(JB2024 != no_soil) %>%
    mutate(crop = cr, version = versions$label[i])
}
names(L)


# ---- 3. Match matrix ---------------------------------------------------------------
match_tbl <- list()
for (k in names(L)) {
  g <- L[[k]]
  for (pl in names(P)) {
    pd  <- P[[pl]]
    idx <- match(g$gridcode, pd$gridcode)
    for (pf in setdiff(names(pd), "gridcode")) {
      d <- pd[[pf]][idx] - g$drain_effekt
      match_tbl[[length(match_tbl) + 1]] <- tibble(
        crop = g$crop[1], version = g$version[1], P_delivery = pl, P_field = pf,
        cells = nrow(g), cells_in_P = sum(!is.na(idx)),
        pct_equal = 100 * mean(abs(d) <= tol, na.rm = TRUE),
        mean_abs_diff = mean(abs(d), na.rm = TRUE),
        max_abs_diff = max(abs(d), na.rm = TRUE))
    }
  }
}
match_tbl <- bind_rows(match_tbl)
write_csv(match_tbl, file.path(out, "P1_match_leaching_vs_P_deliveries.csv"))


# ---- 4. Best match per leaching version ---------------------------------------
best <- match_tbl %>%
  group_by(crop, version) %>%
  slice_max(pct_equal, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(crop, factor(version, levels = versions$label), desc(pct_equal))
print(best, n = Inf, width = Inf)
write_csv(best, file.path(out, "P2_best_match_per_version.csv"))
# A match = pct_equal ~ 100 and cells_in_P = cells.

# Where do the cells that do NOT match the best P delivery sit?
top <- best %>% group_by(crop, version) %>% slice(1) %>% ungroup()
diff_where <- list()
for (i in seq_len(nrow(top))) {
  g  <- L[[paste(top$crop[i], versions$ver[match(top$version[i], versions$label)])]]
  pd <- P[[top$P_delivery[i]]]
  d  <- pd[[top$P_field[i]]][match(g$gridcode, pd$gridcode)] - g$drain_effekt
  off <- g %>% mutate(diff = d) %>% filter(is.na(diff) | abs(diff) > tol)
  diff_where[[i]] <- off %>%
    mutate(profile_group = case_when(AC_jord %in% c(998, 999) ~ "998/999",
                                     TRUE ~ "other profiles")) %>%
    group_by(Georegion, profile_group) %>%
    summarise(cells = n(), median_diff = median(diff, na.rm = TRUE),
              missing_in_P = sum(is.na(diff)), .groups = "drop") %>%
    mutate(crop = top$crop[i], version = top$version[i],
           P_delivery = top$P_delivery[i], P_field = top$P_field[i], .before = 1)
}
diff_where <- bind_rows(diff_where)
print(diff_where, n = Inf)
write_csv(diff_where, file.path(out, "P3_cells_not_matching_best_P.csv"))


# ---- 5. P deliveries compared with each other ------------------------------------
pp <- list()
if (length(P) > 1) for (i in 1:(length(P) - 1)) for (j in (i + 1):length(P)) {
  a <- P[[i]]; b <- P[[j]]
  idx <- match(a$gridcode, b$gridcode)
  for (pf in intersect(setdiff(names(a), "gridcode"), names(b))) {
    d <- b[[pf]][idx] - a[[pf]]
    pp[[length(pp) + 1]] <- tibble(from = names(P)[i], to = names(P)[j], P_field = pf,
                                   cells_common = sum(!is.na(idx)),
                                   only_in_from = sum(is.na(idx)),
                                   only_in_to = sum(!b$gridcode %in% a$gridcode),
                                   pct_changed = 100 * mean(abs(d) > tol, na.rm = TRUE),
                                   median_pct_change = median(100 * d / a[[pf]], na.rm = TRUE))
  }
}
pp <- bind_rows(pp)
print(pp, n = Inf)
write_csv(pp, file.path(out, "P4_P_deliveries_compared.csv"))


# ---- 6. Plot: % of cells equal, leaching version x P delivery (expected P field) --
exp_field <- c(SB = "P_Vaarbygudl", WW = "p_Vhvede")    # check the names printed in section 1
p <- match_tbl %>%
  filter(tolower(P_field) == tolower(exp_field[crop])) %>%
  mutate(version = factor(version, levels = versions$label),
         P_delivery = factor(P_delivery, levels = P_tbl$label),
         crop = recode(crop, SB = "Spring barley + catch crop (P_Vaarbygudl)", WW = "Winter wheat (p_Vhvede)")) %>%
  ggplot(aes(P_delivery, version, fill = pct_equal)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f", pct_equal)), size = 3) +
  scale_fill_distiller(palette = "Greens", direction = 1, limits = c(0, 100), name = "% cells\nequal") +
  facet_wrap(~ crop) +
  labs(x = "P-value delivery", y = "N leaching deliverable",
       title = "Which P delivery is behind the P term of each leaching deliverable?",
       subtitle = "% of cells where drain_effekt equals the delivered P (within 0.0006); JB -99 excluded")
ggsave(file.path(out, "P5_match_leaching_vs_P_deliveries.png"), p, width = 10, height = 4.5, bg = "white")

cat("Done. Results in", out, "\n")
