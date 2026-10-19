
# Brak / NUAR leaching: continuous WINTER WHEAT (VISA quota map)
# Same checks as brak_sample_analysis.R, adapted to winter wheat
#
# Winter wheat: mineral fertiliser at the N norm, no manure, no catch crop
#   -> expected: Nhusdyr = orgN = 0, N_min_fert_foraar = N_norm_jb,
#                crop_effekt = 0, EA = 0, and NO change between versions
#                except P (the imputation), because the October error was manure.
#
# Sections
#   0. Read the wheat sample (made with brak_wheat_sample.R)
#   1. Long format
#   2. What changed between versions
#   3. Checks specific to wheat (norm, no manure, rules)
#   4. Output and input plots, facets by JB2024 and Georegion
#   5. Summary by opland108
#   6. Recompute L (EA = 0)
#   7. Wheat vs spring barley in the same cells (needs the barley sample)


library(data.table)   # fast reading only
library(tidyverse)

base  <- "O:/Tech_AGRO/Klima/cdb/_12784_oplandsmodel/data/Til_eva"
out   <- "C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak"
plots <- file.path(out, "plots_wheat")
dir.create(plots, showWarnings = FALSE)

theme_set(theme_bw(base_size = 10))
ver_cols <- c(may = "#1b9e77", aug = "#7570b3", oct = "#d95f02")


# ---- 0. Read the wheat sample (made with brak_wheat_sample.R) ----------------
wsmp <- fread(file.path(out, "sample10pct_wheat_may_aug.csv"), sep = ";") %>%
  as_tibble()
dim(wsmp)

# Versions present in the sample (from the column suffixes)
versions <- intersect(c("may", "aug", "oct"),
                      unique(str_extract(names(wsmp), "(may|aug|oct)$")))
versions


# ---- 1. Long format -----------------------------------------------------------
w_long <- wsmp %>%
  pivot_longer(-gridcode,
               names_to      = c(".value", "version"),
               names_pattern = "(.*)_(may|aug|oct)$") %>%
  mutate(version   = factor(version, levels = intersect(c("may", "aug", "oct"), versions)),
         JB2024    = factor(JB2024),
         Georegion = factor(Georegion))


# ---- 2. What changed between versions ----------------------------------------
w_changed <- wsmp %>%
  pivot_longer(-gridcode,
               names_to      = c("variable", "version"),
               names_pattern = "(.*)_(may|aug|oct)$") %>%
  pivot_wider(names_from = version, values_from = value) %>%
  group_by(variable) %>%
  summarise(cells     = n(),
            n_may_aug = if (all(c("may", "aug") %in% versions)) sum(abs(aug - may) > 1e-6, na.rm = TRUE) else NA_integer_,
            n_aug_oct = if (all(c("aug", "oct") %in% versions)) sum(abs(oct - aug) > 1e-6, na.rm = TRUE) else NA_integer_,
            .groups = "drop") %>%
  arrange(desc(n_aug_oct), desc(n_may_aug))
print(w_changed, n = Inf)
write_csv2(w_changed, file.path(out, "wheat_summary_changed_variables.csv"))

# Output summary by version, JB2024 and Georegion
w_out_by_jb <- w_long %>%
  group_by(JB2024, version) %>%
  summarise(cells = n(), mean_L = mean(NLES5_NUAR_kgNha), .groups = "drop") %>%
  pivot_wider(names_from = version, values_from = mean_L, names_prefix = "L_")
w_out_by_geo <- w_long %>%
  group_by(Georegion, version) %>%
  summarise(cells = n(), mean_L = mean(NLES5_NUAR_kgNha), .groups = "drop") %>%
  pivot_wider(names_from = version, values_from = mean_L, names_prefix = "L_")
w_out_by_jb; w_out_by_geo
write_csv2(w_out_by_jb,  file.path(out, "wheat_summary_L_by_JB.csv"))
write_csv2(w_out_by_geo, file.path(out, "wheat_summary_L_by_Georegion.csv"))


# ---- 3. Checks specific to winter wheat --------------------------------------
# 3a. Fertilisation = norm, no manure, crop effect 0 (all should be 0 / TRUE)
w_long %>%
  group_by(version) %>%
  summarise(cells          = n(),
            n_Nmin_not_norm = sum(abs(N_min_fert_foraar_kgNha - N_norm_jb) > 0.05),
            max_Nhusdyr    = max(Nhusdyr_kgNHa),
            max_orgN       = max(orgN_husdyr_kgNha),
            crop_effekt    = paste(unique(crop_effekt), collapse = ", "),
            yeareffekt     = paste(unique(yeareffekt), collapse = ", "))

# 3b. Norm by JB (should match Note B Table 1: 171 / 178 / 206 / 219 / 128)
w_long %>%
  filter(version == last(versions)) %>%
  count(JB2024, jb_Nnorm, N_norm_jb)

# 3c. The drain_effekt rule is checked on the full files in brak_wheat_sample.R


# ---- 4. Plots -----------------------------------------------------------------
# 4a. Leaching by version, facet JB2024
p_w_jb <- ggplot(w_long, aes(NLES5_NUAR_kgNha, colour = version)) +
  geom_density() +
  facet_wrap(~ JB2024, scales = "free_y", labeller = label_both) +
  scale_colour_manual(values = ver_cols) +
  labs(x = "Leaching (kg N/ha)", y = "Density", title = "Winter wheat leaching by version, by JB2024")
ggsave(file.path(plots, "w_out_01_density_L_by_JB.png"), p_w_jb, width = 10, height = 7)

# 4b. Leaching by version, facet Georegion
p_w_geo <- ggplot(w_long, aes(NLES5_NUAR_kgNha, colour = version)) +
  geom_density() +
  facet_wrap(~ Georegion, scales = "free_y", labeller = label_both) +
  scale_colour_manual(values = ver_cols) +
  labs(x = "Leaching (kg N/ha)", y = "Density", title = "Winter wheat leaching by version, by Georegion")
ggsave(file.path(plots, "w_out_02_density_L_by_Georegion.png"), p_w_geo, width = 10, height = 5)

# 4c. Boxplots JB2024 x version
p_w_box <- ggplot(w_long, aes(JB2024, NLES5_NUAR_kgNha, fill = version)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  scale_fill_manual(values = ver_cols) +
  labs(x = "JB2024", y = "Leaching (kg N/ha)", title = "Winter wheat leaching by JB2024 and version")
ggsave(file.path(plots, "w_out_03_box_L_JB_version.png"), p_w_box, width = 10, height = 5)

# 4d. Inputs: N norm and N_effekt by JB2024 (should be one value per JB for N_min)
p_w_N <- w_long %>%
  filter(version == last(versions)) %>%
  select(JB2024, N_min_fert_foraar_kgNha, N_effekt, neffekt_power) %>%
  pivot_longer(-JB2024, names_to = "term", values_to = "value") %>%
  ggplot(aes(JB2024, value)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  facet_wrap(~ term, ncol = 1, scales = "free_y") +
  labs(x = "JB2024", y = NULL, title = "Winter wheat N terms by JB2024 (latest version)")
ggsave(file.path(plots, "w_in_01_N_terms_by_JB.png"), p_w_N, width = 9, height = 8)

# 4e. P by Georegion and version
p_w_drain <- ggplot(w_long, aes(Georegion, drain_effekt, fill = version)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  scale_fill_manual(values = ver_cols) +
  labs(x = "Georegion", y = "drain_effekt (P)", title = "Winter wheat P factor by Georegion and version")
ggsave(file.path(plots, "w_in_02_drain_by_Georegion.png"), p_w_drain, width = 8, height = 4)

# 4f. Where P changed May -> Aug (if both versions exist)
if (all(c("may", "aug") %in% versions)) {
  w_pchg <- wsmp %>%
    transmute(gridcode, Georegion = factor(Georegion_aug), JB2024 = factor(JB2024_aug),
              d_drain = drain_effekt_aug - drain_effekt_may,
              pct_L   = 100 * (NLES5_NUAR_kgNha_aug / NLES5_NUAR_kgNha_may - 1)) %>%
    filter(abs(d_drain) > 1e-6)
  cat("Wheat cells where P changed May -> Aug:", nrow(w_pchg), "of", nrow(wsmp), "\n")
  
  p_w_pchg <- ggplot(w_pchg, aes(pct_L)) +
    geom_histogram(bins = 60) +
    facet_wrap(~ Georegion, scales = "free_y", labeller = label_both) +
    labs(x = "% change in leaching, May -> Aug", y = "Cells",
         title = "Winter wheat: cells where drain_effekt changed May -> Aug")
  ggsave(file.path(plots, "w_out_04_pct_may_aug_Pchanged.png"), p_w_pchg, width = 10, height = 5)
}


# ---- 5. Summary by opland108 ---------------------------------------------------
w_by_opl <- w_long %>%
  group_by(opland108, version) %>%
  summarise(cells = n(), mean_L = mean(NLES5_NUAR_kgNha),
            mean_P = mean(drain_effekt), .groups = "drop") %>%
  pivot_wider(names_from = version, values_from = c(mean_L, mean_P)) %>%
  arrange(opland108)
print(w_by_opl, n = 20)
write_csv2(w_by_opl, file.path(out, "wheat_summary_by_opland108.csv"))


# ---- 6. Recompute L (EA = 0 for winter wheat) ----------------------------------
EA_w <- 0
w_rho <- w_long %>%
  filter(drain_effekt > 0, soil_effekt > 0) %>%
  mutate(rho = (NLES5_NUAR_kgNha / (1 - EA_w) - yeareffekt) /
           (neffekt_power * drain_effekt * soil_effekt)) %>%
  group_by(version) %>%
  summarise(rho_median = median(rho), rho_p01 = quantile(rho, 0.01), rho_p99 = quantile(rho, 0.99))
w_rho     # should be ~1.0876, the same as for spring barley

rho_w <- w_rho$rho_median[1]
w_long %>%
  mutate(resid = NLES5_NUAR_kgNha - (yeareffekt + neffekt_power * drain_effekt * soil_effekt * rho_w) * (1 - EA_w)) %>%
  group_by(version) %>%
  summarise(max_abs_resid = max(abs(resid)), n_gt_0.5 = sum(abs(resid) > 0.5))


# ---- 7. Winter wheat vs spring barley in the same cells -------------------------
# Uses the barley sample from brak_sample_analysis.R; matches on gridcode.
sb <- fread(file.path(out, "sample10pct_may_aug_oct.csv"), sep = ";") %>% as_tibble()
lv <- last(versions)     # latest wheat version

ww_sb <- wsmp %>%
  select(gridcode,
         L_ww    = all_of(paste0("NLES5_NUAR_kgNha_", lv)),
         P_ww    = all_of(paste0("drain_effekt_", lv)),
         Nmin_ww = all_of(paste0("N_min_fert_foraar_kgNha_", lv))) %>%
  inner_join(sb %>% select(gridcode, opland108 = opland108_oct,
                           JB2024 = JB2024_oct, Georegion = Georegion_oct,
                           L_sb_aug = NLES5_NUAR_kgNha_aug, L_sb_oct = NLES5_NUAR_kgNha_oct,
                           P_sb = drain_effekt_oct),
             by = "gridcode") %>%
  mutate(JB2024 = factor(JB2024), Georegion = factor(Georegion))
nrow(ww_sb)   # cells in both samples

# Is P different between the two crops? (it should be: different Daisy crop)
ww_sb %>% summarise(share_P_equal = mean(abs(P_ww - P_sb) < 1e-6),
                    mean_P_ww = mean(P_ww), mean_P_sb = mean(P_sb))

p_ww_sb_P <- ggplot(ww_sb, aes(P_sb, P_ww)) +
  geom_bin_2d(bins = 80) +
  geom_abline(linetype = "dashed", colour = "red") +
  facet_wrap(~ Georegion, labeller = label_both) +
  scale_fill_viridis_c(trans = "log10") +
  labs(x = "drain_effekt, spring barley + catch crop", y = "drain_effekt, winter wheat",
       title = "P for the two crops in the same cells")
ggsave(file.path(plots, "w_cmp_01_P_ww_vs_sb.png"), p_ww_sb_P, width = 10, height = 6)

# Leaching ratio SB / WW per catchment, before and after the October fix
ratio_opl <- ww_sb %>%
  group_by(opland108) %>%
  summarise(L_ww = mean(L_ww), L_sb_aug = mean(L_sb_aug), L_sb_oct = mean(L_sb_oct),
            .groups = "drop") %>%
  mutate(ratio_aug = L_sb_aug / L_ww, ratio_oct = L_sb_oct / L_ww)
ratio_opl %>% arrange(desc(ratio_oct)) %>% print(n = 15)
write_csv2(ratio_opl, file.path(out, "wheat_vs_barley_by_opland108.csv"))

p_ratio <- ratio_opl %>%
  select(opland108, ratio_aug, ratio_oct) %>%
  pivot_longer(-opland108, names_to = "barley_version", values_to = "ratio",
               names_prefix = "ratio_") %>%
  ggplot(aes(ratio, fill = barley_version)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  labs(x = "Mean leaching ratio, spring barley + catch crop / winter wheat",
       y = "Catchments", fill = "Barley version",
       title = "How the October fix moves barley leaching relative to wheat")
ggsave(file.path(plots, "w_cmp_02_ratio_sb_ww_by_catchment.png"), p_ratio, width = 8, height = 4)

