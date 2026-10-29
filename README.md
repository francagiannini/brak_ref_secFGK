# brak_ref_secFGK

R scripts behind the DCA advisory paper **"Braklægningspunkt: version control of the NUAR N leaching deliverables"** (Giannini-Kurina & Børgesen, 2026).

Between May and October 2026, AU delivered N leaching data (NUAR, 100 m cells) to SGAV four times, as file geodatabases, for spring barley with catch crop and winter wheat:

- del v1 (29 May)
- del v2 (4 June)
- del v3 (28 July)
- del corr (20 October, spring barley only)

These scripts compare the deliveries cell by cell to document what changed from one version to the next, and why. They also check that the files SEGES used are AU's deliverables, and which P-term (runoff) delivery each one is based on. The data themselves are not included.

## Scripts

| Script | What it does |
|---|---|
| `brak_compare_deliverables.R` | Main comparison, on the full dataset. It lists layers and fields per version and compares v1 → v2 → v3 → corr step by step: share of cells changed per field, size of the change, and a list of changed fields per step. It summarises N leaching and P by Georegion and JB, and produces plots (10% sample), maps of the change per step and the timeline. It also checks the SEGES files (Tif_shp_20260604 and Tif_shp_20260728) against every deliverable, field by field and in the GeoTIFFs. |
| `brak_P_sources.R` | Explains the change in the P term (drain_effekt) from v2 to v3. It checks whether cells changed soil profile or weather grid, whether P is one value per DMI 10 km grid × soil profile, and how much of the change the profile, the grid and the sand/clay class explain. It also gives the change per soil profile. |
| `brak_compare_Pdeliveries.R` | Finds which P-term delivery (29 May, 6 July, 27 July; geodatabase or shapefile) is behind the P term of each N leaching deliverable. It also compares the P deliveries with each other. |

Run `brak_compare_deliverables.R` first; the other two can run in any order. Set the `brak` path at the top of each script to the folder holding the deliveries.

## Notes

- **N leaching:** the field `NLES5_NUAR_kgNha` is the NUAR N leaching and is called *N leaching* in all outputs.
- **Excluded cells:** cells without soil data (JB2024 = −99) are left out of all comparisons.
- **Change thresholds:** a field counts as changed when it moves by more than 0.1%; P is compared within 0.0006 (3 decimals).
- **Requirements:** R ≥ 4.1 with `tidyverse`, `sf`, `terra` and `data.table`.

## License

GPL-3.0, see `LICENSE`.

# AI Disclaimer

The scripts were written with help from Claude (Anthropic), an AI assistant, which drafted and debugged most of the code from the author’s instructions. 
The author designed the analysis, ran every script on the delivered data, checked the outputs, and is responsible for the code and the results.