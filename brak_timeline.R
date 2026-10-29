# =============================================================================
# Brak: timeline of all deliveries and their use
#   lane 1  P-value deliveries (Note A)        29 May, 4 June, 7 July, 3 August
#   lane 2  N leaching deliveries (geodatabase) v1, v2, v3, corr
#   lane 3  Use: SEGES maps and SGAV discharge limits
# Links P -> leaching are drawn from the result of brak_compare_Pdeliveries.R
# (P2_best_match_per_version.csv) when it exists; otherwise from SGAV's stated
# understanding (dashed).
# =============================================================================

library(tidyverse)

brak <- "C:/Users/au710823/OneDrive - Aarhus universitet/PolicySupport/Brak"
out  <- file.path(brak, "compare_deliverables", "plots")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
en_day <- \(d) paste(as.integer(format(d, "%d")), months[as.integer(format(d, "%m"))])
lanes  <- c("Use (SEGES / SGAV)", "N leaching deliveries", "P-value deliveries")

# ---- Events ---------------------------------------------------------------------
P_del <- tibble(
  lane  = "P-value deliveries",
  label = c("P 29 May", "P 4 June", "P 7 July", "P 3 Aug"),
  date  = as.Date(c("2026-05-29", "2026-06-04", "2026-07-07", "2026-08-03")),
  note  = c("data + draft note", "delivery", "revised data", "revised data + note"))

L_del <- tibble(
  lane  = "N leaching deliveries",
  ver   = c("v1", "v2", "v3", "corr"),
  label = c("del v1 (29 May)", "del v2 (4 June)", "del v3 (28 July)", "del corr (20 October)"),
  date  = as.Date(c("2026-05-29", "2026-06-04", "2026-07-28", "2026-10-20")),        # file date
  delivered = as.Date(c("2026-05-29", "2026-06-08", "2026-08-03", "2026-10-26")),
  col   = c("#a6761d", "#1b9e77", "#7570b3", "#d95f02"))

Use <- tibble(
  lane  = "Use (SEGES / SGAV)",
  label = c("SEGES map\nTif_shp 4 June", "Preliminary limits\n26 June / 7 July",
            "SEGES map\nTif_shp 28 July", "Final limits\n29 Sept"),
  date  = as.Date(c("2026-06-04", "2026-06-26", "2026-07-28", "2026-09-29")),
  from  = c("v2", "v2", "v3", "v3"))       # which leaching delivery each one uses

# ---- Links P -> leaching ----------------------------------------------------------
best_file <- file.path(brak, "compare_deliverables", "P_deliveries", "P2_best_match_per_version.csv")
if (file.exists(best_file)) {
  best <- read_csv(best_file, show_col_types = FALSE) %>%
    group_by(crop, version) %>% slice_max(pct_equal, n = 1, with_ties = FALSE) %>% ungroup() %>%
    filter(crop == "SB")                              # spring barley drives the braklægningspunkt
  links <- best %>%
    transmute(ver = L_del$ver[match(version, L_del$label)],
              P_date = as.Date(sprintf("2026-%02d-%02d",
                                       match(str_match(P_delivery, "(\\d{1,2}) (\\w+)$")[, 3], month.name),
                                       as.integer(str_match(P_delivery, "(\\d{1,2}) (\\w+)$")[, 2]))),
              txt = sprintf("%.1f%% equal", pct_equal), confirmed = TRUE)
} else {
  links <- tibble(ver = c("v2", "v3"), P_date = as.Date(c("2026-05-29", "2026-08-03")),
                  txt = "SGAV's understanding (to check)", confirmed = FALSE)
}
links <- links %>% left_join(L_del %>% select(ver, L_date = date), by = "ver")

use_links <- Use %>% left_join(L_del %>% select(from = ver, L_date = date), by = "from")

# ---- Plot ---------------------------------------------------------------------------
y <- set_names(seq_along(lanes), lanes)

p <- ggplot() +
  # lane lines
  geom_segment(data = tibble(yy = y), aes(x = as.Date("2026-05-20"), xend = as.Date("2026-10-31"), y = yy, yend = yy),
               colour = "grey85") +
  # P -> leaching links
  geom_segment(data = links, aes(x = P_date, xend = L_date, y = y[["P-value deliveries"]] - 0.08,
                                 yend = y[["N leaching deliveries"]] + 0.08, linetype = confirmed),
               arrow = arrow(length = unit(2, "mm")), colour = "grey35") +
  geom_text(data = links, aes(x = L_date, y = mean(y[c("P-value deliveries", "N leaching deliveries")]),
                              label = txt), hjust = -0.05, size = 2.8, colour = "grey30") +
  # leaching -> use links
  geom_segment(data = use_links, aes(x = L_date, xend = date, y = y[["N leaching deliveries"]] - 0.08,
                                     yend = y[["Use (SEGES / SGAV)"]] + 0.08),
               arrow = arrow(length = unit(2, "mm")), colour = "grey55") +
  # points
  geom_point(data = P_del, aes(date, y[lane]), shape = 23, size = 3.5, fill = "grey40") +
  geom_point(data = L_del, aes(date, y[lane], fill = label), shape = 21, size = 4.5) +
  geom_point(data = Use, aes(date, y[lane]), shape = 22, size = 3.5, fill = "black") +
  # labels
  geom_text(data = P_del, aes(date, y[lane], label = paste0(label, "\n", note)),
            vjust = -0.8, size = 2.8, lineheight = 0.9) +
  geom_text(data = L_del, aes(date, y[lane], label = paste0(label, "\ndelivered ", en_day(delivered))),
            vjust = 1.8, size = 2.8, lineheight = 0.9) +
  geom_text(data = Use, aes(date, y[lane], label = label), vjust = 1.7, size = 2.8, lineheight = 0.9) +
  scale_fill_manual(values = set_names(L_del$col, L_del$label), guide = "none") +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
  scale_y_continuous(breaks = y, labels = names(y), limits = c(0.5, 3.6)) +
  scale_x_date(breaks = seq(as.Date("2026-06-01"), as.Date("2026-10-01"), by = "month"),
               labels = \(d) months[as.integer(format(d, "%m"))],
               limits = as.Date(c("2026-05-15", "2026-11-05"))) +
  labs(x = NULL, y = NULL,
       title = "P-value deliveries, N leaching deliveries and their use by SEGES and SGAV, 2026") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(out, "08_timeline_P_leaching_use.png"), p, width = 12, height = 5, bg = "white")
p
