# brak_ref_secFGK

R scripts behind the DCA advisory paper **"Braklægningspunkt: version control of the NUAR N leaching deliverables"** (Giannini-Kurina & Børgesen, 2026; file no. 2026-1047352).

The scripts compare, cell by cell, the four N leaching datasets AU delivered to SGAV in 2026. They also check them against the files SEGES used and against the P-term (runoff) deliveries. The data are not part of this repository.

## What the analysis answers

1. **Timeline.** Which datasets were delivered when, and under which file names and dates.
2. **Changes from version to version.** What changed between each delivery (fields, N inputs, P term, N leaching), and where.
3. **Why.** Whether each change came from a revised order from SGAV, an AU method decision, or an AU error correction.
4. **SEGES check.** Whether the two SEGES files (Tif_shp_20260604, Tif_shp_20260728) are exactly AU's deliverables.
5. **P-term check.** Which P-term delivery lies behind the P term (drain_effekt) of each leaching deliverable.

## Deliveries compared

| Label in scripts and figures | File | File date | Delivered | Note |
|---|---|---|---|---|
| del v1 (29 May) | GeoDB_20260529 | 29.05.2026 | 29.05.2026 | First version |
| del v2 (4 June) | GeoDB_20260604 | 04.06.2026 | 08.06.2026 | Second version; SEGES Tif_shp_20260604; preliminary discharge limits |
| del v3 (28 July) | GeoDB_20260728 | 28.07.2026 | 03.08.2026 | Third version; SEGES Tif_shp_20260728; final discharge limits |
| del corr (20 October) | GeoDB_20261020 | 20.10.2026 | 26.10.2026 | Correction of the N allocation; spring barley only |

P-term deliveries: P_GeoDB_20260529, P_GeoDB_20260706 (delivered 07.07) and P_GeoDB_20260727 (delivered 03.08), as geodatabase and/or shapefile.

## Scripts

Run them in this order. Each script starts with a header that lists its sections.

| Script | What it does | Main outputs |
|---|---|---|
| `brak_compare_deliverables.R` | Inventory of layers and fields. Full-dataset comparison step by step (v1 → v2 → v3 → corr): % of cells changed per field and the list of changed fields. Summaries by Georegion and JB. Plots on a 10% sample, maps of the change per step, the SEGES check (shapefiles and GeoTIFFs) and the timeline | `01_` – `06_*.csv`, `plots/` |
| `brak_P_sources.R` | Where the v2 → v3 change in P comes from: cell keys, whether P is one value per DMI grid × soil profile, variance explained by profile / grid / sand-clay, change per profile | `P_sources/P1` – `P7*` |
| `brak_compare_Pdeliveries.R` | Matches drain_effekt in every leaching deliverable against every P field of every P-term delivery (GDB and shapefile), and compares the P deliveries with each other | `P_deliveries/P1` – `P5*` |
| `brak_timeline.R` | Timeline of P-term deliveries, leaching deliveries and their use (SEGES files, SGAV discharge limits) | `plots/08_timeline_P_leaching_use.png` |
| `brak_compare_SEGES_EOB.R` | Byte check (md5) that the SEGES and EOB copies of the SEGES files are identical | `compare_SEGES_EOB/` |

## Data (not included)

Set `brak` at the top of each script to your own folder. The scripts expect this layout:

```
Brak/
├── Delivered_AU/            # N leaching deliveries, unzipped
│   ├── GeoDB_20260529/GeoDB_20260529/Geodatabase_20260529.gdb
│   ├── GeoDB_20260604/...   (same pattern)
│   ├── GeoDB_20260728/...
│   └── GeoDB_20261020/...
├── Delivered_AU_P/          # P-term deliveries, unzipped (GDB and/or shapefiles, any subfolders)
├── SEGES sammenligning til AU/SEGES sammenligning til AU/Input/
│   ├── Tif_shp_20260604/Tif_shp_20260604/   # .tif + *_1_af_2.shp / *_2_af_2.shp
│   └── Tif_shp_20260728/Tif_shp_20260728/
└── compare_deliverables/    # created by the scripts: all results
```

Each geodatabase holds one polygon layer per crop (spring barley with catch crop: `DK_Grid_VaarbygEfterafg*`; winter wheat: `DK_Grid_VinterHvede*`), with one 100 m cell per row and `gridcode` as the key (ETRS89 / UTM 32N).

## Conventions

- **N leaching.** The field `NLES5_NUAR_kgNha` is the NUAR N leaching, and all figures call it *N leaching*.
- **Cells without soil data are excluded.** JB2024 = −99 marks 374,919 cells (0 in every field in v1–v2, empty from v3), so results cover 3,906,813 cells.
- **What counts as a change.** A field counts as changed in a cell when it moves by more than 0.1%. P is compared within 0.0006, because it is delivered with 3 decimals.
- **Full data or sample.** Comparisons, summaries and maps use the full dataset. Scatter and density plots use a 10% sample stratified by coastal catchment (`op_id_int`, seed 42).
- **Shapefile field names** are cut to 10 characters, so they are matched to the GDB names on their first 10 characters, ignoring case.

## Requirements

R ≥ 4.1 with `tidyverse`, `sf`, `terra` and `data.table`. The full comparison holds two versions of one crop in memory at a time (several GB of RAM). Reading the cell centres for the maps takes a while the first time; the result is cached in `gridcode_xy.csv`.

## Citation

Giannini-Kurina F, Børgesen CD. 2026. Braklægningspunkt: version control of the NUAR N leaching deliverables. Advisory paper from DCA – Danish Centre for Food and Agriculture, Aarhus University. Submitted 30.10.2026.

## License

GPL-3.0, see `LICENSE`.
