# SB Brak 
# NUAR / NLES5 leaching outputs (spring barley with catch crop)----
# Read, sample, check consistency and compare May / August / October versions
# Requires: data.table  (install.packages("data.table"))

library(data.table)

# Paths (if the å in the folder names causes trouble, see the note below)
f_may1 <- "O:/Tech_AGRO/Klima/cdb/_12784_oplandsmodel/data/Til_eva/Vaarbyg_vinterhvede060726/Nleaching_output_ data_2021vByg_csemicol.txt"
f_may2 <- "O:/Tech_AGRO/Klima/cdb/_12784_oplandsmodel/data/Til_eva/Nleaching_output_ data_2021SB_ccsemicol.txt"
f_aug <- "O:/Tech_AGRO/Klima/cdb/_12784_oplandsmodel/data/Til_eva/Vårbyg_vinterhvede_270726/Nleaching_output_ data_2021vByg_csemicol.txt"
f_oct <- "O:/Tech_AGRO/Klima/cdb/_12784_oplandsmodel/data/Til_eva/Vårbyg_m_udlæg121026/Nleaching_output_ data_2021vByg_csemicol191026_extracol.txt"

file.exists(c(f_may, f_aug, f_oct))
file.size(c(f_may, f_aug, f_oct)) / 1e6   # MB

# First lines as raw text, to see the separator and decimal mark
readLines(f_may, n = 5)
readLines(f_aug, n = 5)
readLines(f_oct, n = 5)

# First 1000 rows as a table
may <- fread(f_may, sep = ";", fill = TRUE)
aug <- fread(f_aug, sep = ";", fill = TRUE)
oct <- fread(f_oct, sep = ";")

str(may)
str(aug)
str(oct)

may[, V30 := NULL]          # empty column from the trailing ";" in the header
aug[, V30 := NULL]

c(may = nrow(may), aug = nrow(aug), oct = nrow(oct))


# Columns that differ between versions
setdiff(names(oct), names(aug))
setdiff(names(aug), names(may))
setdiff(names(oct), names(may))


# ---- 1. Gridcodes present in all three versions ----
common <- intersect(intersect(may$gridcode, aug$gridcode), oct$gridcode)
length(common)
c(only_may = sum(!may$gridcode %in% common),
  only_aug = sum(!aug$gridcode %in% common),
  only_oct = sum(!oct$gridcode %in% common))

# ---- 2. 10% sample, stratified by catchment (taken from oct) ----
set.seed(42)
pick <- oct[gridcode %in% common, .(gridcode = sample(gridcode, ceiling(0.10 * .N))), by = opland108]$gridcode
length(pick)

# ---- 3. One wide table: same gridcode, three versions side by side ----
s_may <- may[gridcode %in% pick]
s_aug <- aug[gridcode %in% pick]
s_oct <- oct[gridcode %in% pick]

setnames(s_may, setdiff(names(s_may), "gridcode"), paste0(setdiff(names(s_may), "gridcode"), "_may"))
setnames(s_aug, setdiff(names(s_aug), "gridcode"), paste0(setdiff(names(s_aug), "gridcode"), "_aug"))
setnames(s_oct, setdiff(names(s_oct), "gridcode"), paste0(setdiff(names(s_oct), "gridcode"), "_oct"))

smp <- merge(merge(s_may, s_aug, by = "gridcode"), s_oct, by = "gridcode")
dim(smp)
rm(s_may, s_aug, s_oct); gc()

# Save it, so next time you can start from here instead of the 700 MB files
fwrite(smp, "C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak/sample10pct_may_aug_oct20.csv", sep = ";")
# smp <- fread("C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak/sample10pct_may_aug_oct.csv")

# Which inputs/outputs changed, May -> Aug and Aug -> Oct (number of cells)
smp[, .(drain_may_aug = sum(drain_effekt_may != drain_effekt_aug),
        drain_aug_oct = sum(drain_effekt_aug != drain_effekt_oct),
        Nmin_aug_oct  = sum(N_min_fert_foraar_kgNha_aug != N_min_fert_foraar_kgNha_oct),
        L_may_aug     = sum(NLES5_NUAR_kgNha_may != NLES5_NUAR_kgNha_aug),
        L_aug_oct     = sum(NLES5_NUAR_kgNha_aug != NLES5_NUAR_kgNha_oct))]

head(smp)


# ---- 1. Every column: how many cells changed, May->Aug and Aug->Oct ----
cols <- setdiff(names(may), "gridcode")
res <- data.table(column = cols, may_aug = 0L, aug_oct = 0L)
for (cn in cols) {
  a <- smp[[paste0(cn, "_may")]]; b <- smp[[paste0(cn, "_aug")]]; c <- smp[[paste0(cn, "_oct")]]
  res[column == cn, may_aug := sum(abs(b - a) > 1e-6, na.rm = TRUE)]
  res[column == cn, aug_oct := sum(abs(c - b) > 1e-6, na.rm = TRUE)]
}
res

# ---- 2. The 8,604 cells where P changed May -> Aug ----
p_chg <- smp[abs(drain_effekt_aug - drain_effekt_may) > 1e-6]
p_chg[, .N, by = AC_jord_aug][order(-N)]                    # is it mostly AC_jord 999?
p_chg[, .N, by = Georegion_aug][order(-N)]
p_chg[, summary(drain_effekt_aug - drain_effekt_may)]
p_chg[, summary(100 * (NLES5_NUAR_kgNha_aug / NLES5_NUAR_kgNha_may - 1))]   # % change in L

# ---- 3. The 43 cells where P changed Aug -> Oct (should be 0) ----
smp[abs(drain_effekt_oct - drain_effekt_aug) > 1e-6,
    .(gridcode, AC_jord_aug, AC_jord_oct, DMI_10km_oct, drain_effekt_may, drain_effekt_aug, drain_effekt_oct)]

# ---- 4. The 209 cells where N_min did NOT change Aug -> Oct ----
smp[N_min_fert_foraar_kgNha_aug == N_min_fert_foraar_kgNha_oct,
    .N, by = .(opland108_oct, ammonia_orgNfert_oct, N_min_fert_foraar_kgNha_oct)]

# ---- 5. Is the change in N_min exactly the ammonia? ----
smp[, d_min := N_min_fert_foraar_kgNha_oct - N_min_fert_foraar_kgNha_aug]
smp[, summary(d_min - ammonia_orgNfert_oct)]
smp[abs(d_min - ammonia_orgNfert_oct) > 0.1, .N]
smp[abs(d_min - ammonia_orgNfert_oct) > 0.1,
    .(gridcode, opland108_oct, JB2024_oct, N_norm_jb_oct, N_min_fert_foraar_kgNha_aug,
      N_min_fert_foraar_kgNha_oct, ammonia_orgNfert_oct)][1:20]

# ---- 6. Christen's checks on the sample ----
smp[, .(n_Nhusdyr = uniqueN(Nhusdyr_kgNHa_oct), n_orgN = uniqueN(orgN_husdyr_kgNha_oct),
        n_ammonia = uniqueN(ammonia_orgNfert_oct)), by = opland108_oct][
          n_Nhusdyr > 1 | n_orgN > 1 | n_ammonia > 1]
smp[, .(n_drain = uniqueN(drain_effekt_oct)), by = .(DMI_10km_oct, AC_jord_oct)][n_drain > 1, .N]
smp[, .(n_drain = uniqueN(drain_effekt_may)), by = .(DMI_10km_may, AC_jord_may)][n_drain > 1, .N]

# ---- 7. The 54,807 cells that are new in Aug (full file, not the sample) ----
new_aug <- aug[!gridcode %in% may$gridcode]
new_aug[, .N, by = AC_jord][order(-N)][1:15]
new_aug[, .N, by = Georegion][order(-N)]
new_aug[, summary(drain_effekt)]


# Brak
# WW NUAR / NLES5 leaching outputs (continuous WINTER WHEAT)----
# Read, sample and save the May / August / October versions
# Same steps as the spring barley sample; same gridcodes as the barley sample
# Requires: data.table  (install.packages("data.table"))

# Brak
# NUAR / NLES5 leaching outputs (continuous WINTER WHEAT)----
# Read, sample and save the May / August versions
# (there is no October winter wheat file: the October rerun was spring barley only)
# Same steps as the spring barley sample; same gridcodes as the barley sample
# Requires: data.table  (install.packages("data.table"))

library(data.table)

# Paths: same folders as spring barley, wheat file name
f_may <- "O:/Tech_AGRO/Klima/cdb/_12784_oplandsmodel/data/Til_eva/Vaarbyg_vinterhvede060726/Nleaching_output_ data_2021Wheatsemicol.txt"
f_aug <- "O:/Tech_AGRO/Klima/cdb/_12784_oplandsmodel/data/Til_eva/Vårbyg_vinterhvede_270726/Nleaching_output_ data_2021Wheatsemicol.txt"
# f_oct <- "O:/Tech_AGRO/Klima/cdb/_12784_oplandsmodel/data/Til_eva/Vårbyg_m_udlæg121026/Nleaching_output_ data_2021Wheatsemicol.txt"   # does not exist

file.exists(c(f_may, f_aug))
file.size(c(f_may, f_aug)) / 1e6   # MB

# If a path is FALSE, find the real names here:
# list.files("O:/Tech_AGRO/Klima/cdb/_12784_oplandsmodel/data/Til_eva",
#            pattern = "Wheat", recursive = TRUE, full.names = TRUE)

# First lines as raw text, to see the separator, decimal mark and header
readLines(f_may, n = 5)
readLines(f_aug, n = 5)
# readLines(f_oct, n = 5)

# Full files
w_may <- fread(f_may, sep = ";", fill = TRUE)
w_aug <- fread(f_aug, sep = ";", fill = TRUE)
# w_oct <- fread(f_oct, sep = ";", fill = TRUE)

str(w_may)
str(w_aug)
# str(w_oct)

# Remove the empty column from a trailing ";" in the header, where it exists
if ("V30" %in% names(w_may)) w_may[, V30 := NULL]
if ("V30" %in% names(w_aug)) w_aug[, V30 := NULL]
# if ("V30" %in% names(w_oct)) w_oct[, V30 := NULL]

c(may = nrow(w_may), aug = nrow(w_aug))

# Columns that differ between versions
setdiff(names(w_aug), names(w_may))
setdiff(names(w_may), names(w_aug))
# setdiff(names(w_oct), names(w_aug))


# ---- 1. Gridcodes present in both versions ----
common_w <- intersect(w_may$gridcode, w_aug$gridcode)
# common_w <- intersect(common_w, w_oct$gridcode)
length(common_w)
c(only_may = sum(!w_may$gridcode %in% common_w),
  only_aug = sum(!w_aug$gridcode %in% common_w))
# only_oct = sum(!w_oct$gridcode %in% common_w)

# ---- 2. Sample: the SAME gridcodes as the spring barley sample ----
# so wheat and barley can be compared cell by cell later
smp_sb <- fread("C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak/sample10pct_may_aug_oct.csv",
                sep = ";", select = "gridcode")
pick_w <- intersect(smp_sb$gridcode, common_w)
length(pick_w)                                  # expected close to 385,252
length(setdiff(smp_sb$gridcode, common_w))      # barley sample cells missing in wheat

# (Alternative: a new 10% sample stratified by catchment, as for barley)
# set.seed(42)
# pick_w <- w_aug[gridcode %in% common_w, .(gridcode = sample(gridcode, ceiling(0.10 * .N))), by = opland108]$gridcode

# ---- 3. One wide table: same gridcode, two versions side by side ----
s_may <- w_may[gridcode %in% pick_w]
s_aug <- w_aug[gridcode %in% pick_w]
# s_oct <- w_oct[gridcode %in% pick_w]

setnames(s_may, setdiff(names(s_may), "gridcode"), paste0(setdiff(names(s_may), "gridcode"), "_may"))
setnames(s_aug, setdiff(names(s_aug), "gridcode"), paste0(setdiff(names(s_aug), "gridcode"), "_aug"))
# setnames(s_oct, setdiff(names(s_oct), "gridcode"), paste0(setdiff(names(s_oct), "gridcode"), "_oct"))

wsmp <- merge(s_may, s_aug, by = "gridcode")
# wsmp <- merge(wsmp, s_oct, by = "gridcode")
dim(wsmp)
rm(s_may, s_aug); gc()

# ---- 4. drain_effekt rule on the FULL files (one value per DMI_10km x AC_jord) ----
w_may[, .(n_drain = uniqueN(round(drain_effekt, 5))), by = .(DMI_10km, AC_jord)][n_drain > 1, .N]
w_aug[, .(n_drain = uniqueN(round(drain_effekt, 5))), by = .(DMI_10km, AC_jord)][n_drain > 1, .N]
# w_oct[, .(n_drain = uniqueN(round(drain_effekt, 5))), by = .(DMI_10km, AC_jord)][n_drain > 1, .N]

# Save the sample, so the analysis can start from here instead of the 700 MB files
fwrite(wsmp, "C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak/sample10pct_wheat_may_aug.csv", sep = ";")
# wsmp <- fread("C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak/sample10pct_wheat_may_aug.csv")