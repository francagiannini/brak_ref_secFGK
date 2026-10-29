# =============================================================================
# Brak / NUAR leaching: spring barley with catch crop
# Comparison of three versions (may = 060726, aug = 270726, oct = 121026)
# on the 10% sample (same gridcodes in all three, stratified by opland108)
#
# Sections
#   0. Setup and read the sample
#   1. Long format (one row per cell x version)
#   2. General summary: what changed, outputs and inputs per version
#   3. OUTPUT plots (NLES5_NUAR_kgNha), facets by JB2024 and Georegion
#   4. INPUT plots (N terms, drain/P, soil/S, manure), facets by JB2024 and Georegion
#   5. Summary by opland108 (coastal catchment) + plots (+ optional map)
#   6. Christen's consistency checks
#   7. Recompute L from its terms (does the equation close?)
# =============================================================================

library(data.table)   # only for fast reading
library(tidyverse)

out   <- "C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak"
plots <- file.path(out, "plots2")
dir.create(plots, showWarnings = FALSE)

theme_set(theme_bw(base_size = 10))
ver_cols <- c(may = "#1b9e77", aug = "#7570b3", oct = "#d95f02")


# ---- 0. Read the sample ------------------------------------------------------
smp <- fread(file.path(out, "sample10pct_may_aug_oct20.csv"), sep = ";") %>%
  as_tibble()

dim(smp)   # expected: 385252 x 86


# ---- 1. Long format: one row per cell x version ------------------------------
# Columns are named <variable>_<version>; tidyr splits them back.
smp_long <- smp %>%
  pivot_longer(-gridcode,
               names_to      = c(".value", "version"),
               names_pattern = "(.*)_(may|aug|oct)$") %>%
  mutate(version   = factor(version, levels = c("may", "aug", "oct")),
         JB2024    = factor(JB2024),
         Georegion = factor(Georegion)) |> filter(!is.na(version))
# ammonia_orgNfert only exists in oct; it is NA for may and aug

# Cell-level table of changes between versions (wide)
chg <- smp %>%
  transmute(gridcode,
            opland108 = opland108_oct,
            JB2024    = factor(JB2024_oct),
            Georegion = factor(Georegion_oct),
            AC_jord   = AC_jord_oct,
            ammonia   = ammonia_orgNfert_oct,
            d_drain_may_aug = drain_effekt_aug - drain_effekt_may,
            d_drain_aug_oct = drain_effekt_oct - drain_effekt_aug,
            d_Nmin_aug_oct  = N_min_fert_foraar_kgNha_oct - N_min_fert_foraar_kgNha_aug,
            d_Neff_aug_oct  = N_effekt_oct - N_effekt_aug,
            L_may = NLES5_NUAR_kgNha_may,
            L_aug = NLES5_NUAR_kgNha_aug,
            L_oct = NLES5_NUAR_kgNha_oct,
            d_L_may_aug   = L_aug - L_may,
            d_L_aug_oct   = L_oct - L_aug,
            pct_L_may_aug = 100 * (L_aug / L_may - 1),
            pct_L_aug_oct = 100 * (L_oct / L_aug - 1),
            P_changed     = abs(d_drain_may_aug) > 1e-6)


# ---- 2. General summary ------------------------------------------------------

# 2a. Which variables changed between versions, and in how many cells
#     (variable x version table, then compare the version columns)
# (this builds a long table of ~32 M rows; needs ~2 GB RAM)
changed <- smp %>%
  select(gridcode, ends_with(c("_may", "_aug", "_oct"))) %>%
  select(-starts_with("ammonia_orgNfert")) %>%
  pivot_longer(-gridcode,
               names_to      = c("variable", "version"),
               names_pattern = "(.*)_(may|aug|oct)$") %>%
  pivot_wider(names_from = version, values_from = value) %>%
  group_by(variable) %>%
  summarise(cells          = n(),
            n_may_aug      = sum(abs(aug - may) > 1e-6, na.rm = TRUE),
            n_aug_oct      = sum(abs(oct - aug) > 1e-6, na.rm = TRUE),
            mean_d_may_aug = mean(aug - may, na.rm = TRUE),
            mean_d_aug_oct = mean(oct - aug, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(pct_may_aug = round(100 * n_may_aug / cells, 2),
         pct_aug_oct = round(100 * n_aug_oct / cells, 2)) %>%
  arrange(desc(n_aug_oct), desc(n_may_aug))

print(changed, n = Inf)
write_csv2(changed, file.path(out, "summary_changed_variables.csv"))

# 2b. Output summary per version (whole sample, by JB, by Georegion)
out_summary <- smp_long %>%
  group_by(version) %>%
  summarise(cells = n(),
            mean  = mean(NLES5_NUAR_kgNha, na.rm = T), sd = sd(NLES5_NUAR_kgNha, na.rm = T),
            p05   = quantile(NLES5_NUAR_kgNha, 0.05, na.rm = T),
            median = median(NLES5_NUAR_kgNha, na.rm = T),
            p95   = quantile(NLES5_NUAR_kgNha, 0.95, na.rm = T),
            .groups = "drop")
out_summary

out_by_jb <- smp_long %>%
  group_by(JB2024, version) %>%
  summarise(cells = n(), mean_L = mean(NLES5_NUAR_kgNha, na.rm = T), .groups = "drop") %>%
  pivot_wider(names_from = version, values_from = mean_L, names_prefix = "L_") %>%
  mutate(pct_may_aug = 100 * (L_aug / L_may - 1),
         pct_aug_oct = 100 * (L_oct / L_aug - 1))
out_by_jb

out_by_geo <- smp_long %>%
  group_by(Georegion, version) %>%
  summarise(cells = n(), mean_L = mean(NLES5_NUAR_kgNha), .groups = "drop") %>%
  pivot_wider(names_from = version, values_from = mean_L, names_prefix = "L_") %>%
  mutate(pct_may_aug = 100 * (L_aug / L_may - 1),
         pct_aug_oct = 100 * (L_oct / L_aug - 1))
out_by_geo

write_csv2(out_by_jb,  file.path(out, "summary_L_by_JB.csv"))
write_csv2(out_by_geo, file.path(out, "summary_L_by_Georegion.csv"))

# 2c. Input summary per version (mean of each model term)
in_summary <- smp_long %>%
  group_by(version) %>%
  summarise(across(c(N_min_fert_foraar_kgNha, Nhusdyr_kgNHa, orgN_husdyr_kgNha,
                     org_N_2024_tNha, N_effekt, neffekt_power, crop_effekt,
                     drain_effekt, soil_effekt, yeareffekt),
                   mean),
            .groups = "drop")
in_summary


# ---- 3. OUTPUT plots: NLES5_NUAR_kgNha ---------------------------------------

# 3a. Distribution per version, facet JB2024
p_out_jb <- ggplot(smp_long, aes(NLES5_NUAR_kgNha, colour = version)) +
  geom_density() +
  facet_wrap(~ JB2024, scales = "free_y", labeller = label_both) +
  scale_colour_manual(values = ver_cols) +
  labs(x = "Leaching (kg N/ha)", y = "Density",
       title = "Leaching per cell by version, by JB2024")
ggsave(file.path(plots, "out_01_density_L_by_JB.png"), p_out_jb, width = 10, height = 7)

# 3b. Distribution per version, facet Georegion
p_out_geo <- ggplot(smp_long, aes(NLES5_NUAR_kgNha, colour = version)) +
  geom_density() +
  facet_wrap(~ Georegion, scales = "free_y", labeller = label_both) +
  scale_colour_manual(values = ver_cols) +
  labs(x = "Leaching (kg N/ha)", y = "Density",
       title = "Leaching per cell by version, by Georegion")
ggsave(file.path(plots, "out_02_density_L_by_Georegion.png"), p_out_geo, width = 10, height = 5)

# 3c. Boxplots, JB2024 x version
p_out_box <- ggplot(smp_long, aes(JB2024, NLES5_NUAR_kgNha, fill = version)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  scale_fill_manual(values = ver_cols) +
  labs(x = "JB2024", y = "Leaching (kg N/ha)", title = "Leaching by JB2024 and version")
ggsave(file.path(plots, "out_03_box_L_JB_version.png"), p_out_box, width = 10, height = 5)

# 3d. Aug vs Oct per cell, facet JB2024
p_out_scatter <- ggplot(chg, aes(L_aug, L_oct)) +
  geom_bin_2d(bins = 80) +
  geom_abline(linetype = "dashed", colour = "red") +
  facet_wrap(~ JB2024, labeller = label_both) +
  scale_fill_viridis_c(trans = "log10") +
  labs(x = "Aug leaching (kg N/ha)", y = "Oct leaching (kg N/ha)",
       title = "Aug vs Oct leaching per cell, by JB2024")
ggsave(file.path(plots, "out_04_scatter_aug_oct_by_JB.png"), p_out_scatter, width = 10, height = 8)

# 3e. % change Aug -> Oct, by JB2024 and by Georegion
p_out_pct_jb <- ggplot(chg, aes(JB2024, pct_L_aug_oct)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  labs(x = "JB2024", y = "% change in leaching, Aug -> Oct",
       title = "Effect of adding manure ammonium, by JB2024")
ggsave(file.path(plots, "out_05_pct_aug_oct_by_JB.png"), p_out_pct_jb, width = 8, height = 4)

p_out_pct_geo <- ggplot(chg, aes(Georegion, pct_L_aug_oct)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  labs(x = "Georegion", y = "% change in leaching, Aug -> Oct",
       title = "Effect of adding manure ammonium, by Georegion")
ggsave(file.path(plots, "out_06_pct_aug_oct_by_Georegion.png"), p_out_pct_geo, width = 6, height = 4)

# 3f. May -> Aug: only cells where P changed
p_out_mayaug <- chg %>%
  filter(P_changed) %>%
  ggplot(aes(pct_L_may_aug)) +
  geom_histogram(bins = 60) +
  facet_wrap(~ Georegion, scales = "free_y", labeller = label_both) +
  labs(x = "% change in leaching, May -> Aug", y = "Cells",
       title = "Cells where drain_effekt changed May -> Aug (P imputation)")
ggsave(file.path(plots, "out_07_pct_may_aug_Pchanged_by_Georegion.png"), p_out_mayaug, width = 10, height = 5)


# ---- 4. INPUT plots ----------------------------------------------------------

# 4a. N chain: N_min -> N_effekt -> neffekt_power, by version and JB2024
p_in_N <- smp_long %>%
  select(version, JB2024, N_min_fert_foraar_kgNha, N_effekt, neffekt_power) %>%
  pivot_longer(-c(version, JB2024), names_to = "term", values_to = "value") %>%
  ggplot(aes(JB2024, value, fill = version)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  facet_wrap(~ term, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = ver_cols) +
  labs(x = "JB2024", y = NULL, title = "N terms by JB2024 and version")
ggsave(file.path(plots, "in_01_N_terms_by_JB.png"), p_in_N, width = 10, height = 9)

# 4b. Change in N_min Aug -> Oct vs ammonia added (should lie on the 1:1 line)
p_in_Nmin_ammonia <- ggplot(chg, aes(ammonia, d_Nmin_aug_oct)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_abline(linetype = "dashed", colour = "red") +
  facet_wrap(~ JB2024, labeller = label_both) +
  labs(x = "ammonia_orgNfert (kg N/ha)", y = "N_min_fert_foraar Oct - Aug (kg N/ha)",
       title = "Is the October change in mineral N exactly the manure ammonium?")
ggsave(file.path(plots, "in_02_dNmin_vs_ammonia_by_JB.png"), p_in_Nmin_ammonia, width = 10, height = 8)

# 4c. drain_effekt (P) by version and Georegion
p_in_drain <- ggplot(smp_long, aes(Georegion, drain_effekt, fill = version)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  scale_fill_manual(values = ver_cols) +
  labs(x = "Georegion", y = "drain_effekt (P)", title = "P factor by Georegion and version")
ggsave(file.path(plots, "in_03_drain_by_Georegion.png"), p_in_drain, width = 8, height = 4)

# 4d. Where P changed May -> Aug: Georegion x JB2024, and soil profile 999
p_in_Pchg <- chg %>%
  count(Georegion, JB2024, P_changed) %>%
  group_by(Georegion, JB2024) %>%
  mutate(share = n / sum(n)) %>%
  filter(P_changed) %>%
  ggplot(aes(JB2024, Georegion, fill = 100 * share)) +
  geom_tile() +
  geom_text(aes(label = n), size = 2.5, colour = "white") +
  scale_fill_viridis_c(name = "% of cells") +
  labs(title = "Share of cells where P changed May -> Aug (label = cells)")
ggsave(file.path(plots, "in_04_Pchanged_tile_Georegion_JB.png"), p_in_Pchg, width = 9, height = 4)

p_in_Pchg_999 <- chg %>%
  filter(P_changed) %>%
  mutate(AC_jord_999 = AC_jord == 999) %>%
  ggplot(aes(d_drain_may_aug, fill = AC_jord_999)) +
  geom_histogram(bins = 60) +
  facet_wrap(~ Georegion, scales = "free_y", labeller = label_both) +
  labs(x = "drain_effekt Aug - May", y = "Cells", fill = "AC_jord = 999",
       title = "Size of the P change May -> Aug")
ggsave(file.path(plots, "in_05_dDrain_may_aug.png"), p_in_Pchg_999, width = 10, height = 5)

# 4e. Soil: S factor vs clay, coloured by JB2024 (versions should be identical)
p_in_soil <- smp_long %>%
  filter(version == "oct") %>%
  ggplot(aes(ler_tex2024, soil_effekt, colour = JB2024)) +
  geom_point(alpha = 0.2, size = 0.4) +
  labs(x = "Clay, ler_tex2024 (%)", y = "soil_effekt (S)",
       title = "S = exp(-zeta * clay); colour shows JB/clay mismatches")
ggsave(file.path(plots, "in_06_soil_vs_clay.png"), p_in_soil, width = 8, height = 5)

# 4f. Topsoil organic N vs ctotal / CN, by JB2024 (floor at 3, cap at 10)
p_in_NT <- smp_long %>%
  filter(version == "oct") %>%
  ggplot(aes(ctotal2024 / CN_ration, org_N_2024_tNha)) +
  geom_point(alpha = 0.1, size = 0.4) +
  geom_abline(linetype = "dashed", colour = "red") +
  facet_wrap(~ JB2024, labeller = label_both) +
  labs(x = "ctotal2024 / CN_ration", y = "org_N_2024_tNha (NT)",
       title = "NT = C / (C:N), limited to 3-10 t N/ha")
ggsave(file.path(plots, "in_07_NT_check_by_JB.png"), p_in_NT, width = 10, height = 8)


# ---- 5. Summary by opland108 (coastal catchment) -----------------------------
by_opl <- chg %>%
  group_by(opland108) %>%
  summarise(cells       = n(),
            ammonia     = first(ammonia),
            L_may       = mean(L_may),
            L_aug       = mean(L_aug),
            L_oct       = mean(L_oct),
            pct_may_aug = 100 * (L_aug / L_may - 1),
            pct_aug_oct = 100 * (L_oct / L_aug - 1),
            n_P_changed = sum(P_changed),
            n_Nmin_same = sum(abs(d_Nmin_aug_oct) < 1e-6),
            .groups = "drop") %>%
  left_join(smp %>%
              group_by(opland108 = opland108_oct) %>%
              summarise(Nhusdyr   = first(Nhusdyr_kgNHa_oct),
                        orgN      = first(orgN_husdyr_kgNha_oct),
                        Nmin_aug  = mean(N_min_fert_foraar_kgNha_aug),
                        Nmin_oct  = mean(N_min_fert_foraar_kgNha_oct),
                        main_geo  = names(which.max(table(Georegion_oct))),
                        .groups = "drop"),
            by = "opland108") %>%
  arrange(desc(pct_aug_oct))

print(by_opl, n = 20)
write_csv2(by_opl, file.path(out, "summary_by_opland108.csv"))

# 5a. Mean leaching per catchment, three versions
p_opl_dumbbell <- by_opl %>%
  select(opland108, L_may, L_aug, L_oct) %>%
  pivot_longer(-opland108, names_to = "version", values_to = "L",
               names_prefix = "L_") %>%
  mutate(version   = factor(version, levels = c("may", "aug", "oct")),
         opland108 = fct_reorder(factor(opland108), L)) %>%
  ggplot(aes(L, opland108, colour = version)) +
  geom_line(aes(group = opland108), colour = "grey70") +
  geom_point(size = 1.5) +
  scale_colour_manual(values = ver_cols) +
  labs(x = "Mean leaching (kg N/ha)", y = "opland108",
       title = "Mean leaching per coastal catchment") +
  theme(axis.text.y = element_text(size = 5))
ggsave(file.path(plots, "opl_01_dumbbell.png"), p_opl_dumbbell, width = 6, height = 12)

# 5b. % change vs manure ammonium, coloured by main Georegion
p_opl_ammonia <- ggplot(by_opl, aes(ammonia, pct_aug_oct, colour = main_geo)) +
  geom_point(size = 2) +
  geom_text(aes(label = opland108), size = 2.3, vjust = -0.8, check_overlap = TRUE,
            show.legend = FALSE) +
  labs(x = "ammonia_orgNfert (kg N/ha)", y = "% change in mean leaching, Aug -> Oct",
       colour = "Main Georegion",
       title = "Catchments with more manure ammonium change most")
ggsave(file.path(plots, "opl_02_ammonia_vs_pct.png"), p_opl_ammonia, width = 8, height = 5)

# 5c. Manure composition per catchment (organic vs ammonium)
p_opl_manure <- by_opl %>%
  mutate(opland108 = fct_reorder(factor(opland108), Nhusdyr)) %>%
  select(opland108, orgN, ammonia) %>%
  pivot_longer(-opland108, names_to = "part", values_to = "kgN") %>%
  ggplot(aes(kgN, opland108, fill = part)) +
  geom_col() +
  labs(x = "Manure N (kg N/ha)", y = "opland108", fill = NULL,
       title = "Manure N per catchment: organic + ammonium = Nhusdyr") +
  theme(axis.text.y = element_text(size = 5))
ggsave(file.path(plots, "opl_03_manure_parts.png"), p_opl_manure, width = 6, height = 12)

# 5d. OPTIONAL map by catchment (needs the 108 coastal catchment polygons)
# library(sf)
# opl_sf <- st_read("PATH/TO/kystvandoplande.shp")
# names(opl_sf)                                             # find the ID matching opland108
# map_dat <- opl_sf %>% left_join(by_opl, by = c("ID_FIELD" = "opland108"))
# p_map <- ggplot(map_dat) +
#   geom_sf(aes(fill = pct_aug_oct), colour = "white", linewidth = 0.1) +
#   scale_fill_viridis_c(name = "% change\nAug -> Oct") +
#   theme_void()
# ggsave(file.path(plots, "map_pct_aug_oct.png"), p_map, width = 7, height = 7)


# ---- 6. Christen's checks (on the sample) ------------------------------------
# 6a. One value per opland108 for Nhusdyr, orgN (all versions) and ammonia (oct)
check_manure <- smp_long %>%
  group_by(version, opland108) %>%
  summarise(n_Nhusdyr = n_distinct(round(Nhusdyr_kgNHa, 5)),
            n_orgN    = n_distinct(round(orgN_husdyr_kgNha, 5)),
            n_ammonia = n_distinct(round(ammonia_orgNfert, 5)),
            .groups = "drop")

check_manure %>%
  group_by(version) %>%
  summarise(catchments     = n(),
            bad_Nhusdyr    = sum(n_Nhusdyr > 1),
            bad_orgN       = sum(n_orgN > 1),
            bad_ammonia    = sum(n_ammonia > 1))

# 6b. ammonia = Nhusdyr - orgN (oct)
smp %>%
  summarise(max_abs_diff = max(abs(ammonia_orgNfert_oct -
                                     (Nhusdyr_kgNHa_oct - orgN_husdyr_kgNha_oct))))

# 6c. One drain_effekt per DMI_10km x AC_jord
check_drain <- smp_long %>%
  group_by(version, DMI_10km, AC_jord) %>%
  summarise(n_drain = n_distinct(round(drain_effekt, 5)), cells = n(), .groups = "drop")

check_drain %>%
  group_by(version) %>%
  summarise(combinations = n(), bad = sum(n_drain > 1), cells_in_bad = sum(cells[n_drain > 1]))

check_drain %>% filter(n_drain > 1) %>% arrange(desc(n_drain))

# 6d. The cells that break the "October only changes N" expectation
chg %>% filter(abs(d_drain_aug_oct) > 1e-6)                       # P changed Aug -> Oct
chg %>% filter(abs(d_Nmin_aug_oct - ammonia) > 0.1)               # N_min change != ammonia


# ---- 7. Recompute L from its terms --------------------------------------------
# L = (yeareffekt + neffekt_power * drain_effekt * soil_effekt * rho) * (1 - EA)
# EA = 0.45 for spring barley with catch crop. rho is estimated from the data.
EA <- 0.45

rho_est <- smp_long %>%
  filter(drain_effekt > 0, soil_effekt > 0) %>%
  mutate(rho = (NLES5_NUAR_kgNha / (1 - EA) - yeareffekt) /
           (neffekt_power * drain_effekt * soil_effekt)) %>%
  group_by(version) %>%
  summarise(rho_median = median(rho, na.rm = TRUE),
            rho_p01 = quantile(rho, 0.01, na.rm = TRUE),
            rho_p99 = quantile(rho, 0.99, na.rm = TRUE))
rho_est   # one constant across cells and versions = the equation closes

rho <- rho_est$rho_median[rho_est$version == "oct"]

recalc <- smp_long %>%
  mutate(L_hat = (yeareffekt + neffekt_power * drain_effekt * soil_effekt * rho) * (1 - EA),
         resid = NLES5_NUAR_kgNha - L_hat)

recalc %>%
  group_by(version) %>%
  summarise(mean_resid = mean(resid), max_abs_resid = max(abs(resid)),
            n_resid_gt_0.5 = sum(abs(resid) > 0.5))

p_recalc <- ggplot(recalc, aes(resid)) +
  geom_histogram(bins = 100) +
  facet_grid(version ~ Georegion, labeller = label_both, scales = "free_y") +
  labs(x = "Reported - recomputed leaching (kg N/ha)", y = "Cells",
       title = "Does L close from its delivered terms? (neffekt_power has 1 decimal)")
ggsave(file.path(plots, "chk_01_recompute_residuals.png"), p_recalc, width = 11, height = 6)

