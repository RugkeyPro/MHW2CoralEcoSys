# Demo input dictionary and selection

Runtime input: `analysis/modules/conservation_demo/data/reef_cells.csv.gz`, a gzip-compressed UTF-8 CSV with 29,746 real cells and 25 fields. This is a derived scientific input dataset, not a synthetic example or a final output table. The full population of positive-area cells from the accepted archived input is retained; there is no filtering by the result of the prioritization.

| Fields | Meaning / units |
| --- | --- |
| `cell_id` | Unique source reef-grid identifier; unchanged |
| `row`, `col` | Original grid indices, retained for provenance |
| `longitude`, `latitude` | Original cell coordinates in EPSG:4326 degrees |
| `reef_area_km2` | Positive reef area per grid cell, km² |
| `mpa_fraction_with_point_buffers` | Existing MPA reef fraction, including parent point-buffer treatment; [0,1] |
| `P_valid_ssp*_2050` | True when all three original climate-member support masks are finite |
| `P_ssp*_2050` | Mean of the three OR10/MESS member masks; not a selection probability |
| `mpei_total` | Mean 2045–2055 36-tracer 0–100 m MPEI mapped to the reef cell, kg m⁻³ |
| `mpei_supported` | Whether the parent field supports a valid mapped MPEI value |
| `mpei_e_star_global_reef_minmax` | Parent global-reef min-max normalization of MPEI; secondary ranking only |
| `ssp*_2050__member` | Nine exact archived member support masks: 0, 1 or NaN |

Scenarios are SSP126/245/585; exact member labels are GFDL-ESM4 `r1i1p1f1_gr`, CNRM-ESM2-1 `r1i1p1f2_gn` and IPSL-CM6A-LR `r1i1p1f1_gn`. The same future MPEI field is paired with all scenarios; plastic trajectories are not independently SSP-specific.

Missing support is not zero support. Unsupported MPEI cells can retain spatial/habitat selection, but do not enter a supported-MPEI numerator or denominator, according to the original calculation. Reproduction requires retaining all rows, original row order and the original global normalization.

The source Parquet files include geometry and additional components. Their WKB boundary geometry and unused components are excluded from this numeric demo. Original third-party reef/MPA raw datasets retain their original reuse requirements; the software MIT license does not grant rights to those raw datasets.

Frozen results under `reference/` are for comparison, not analysis input. `provenance.json` records their exact bytes or lossless gzip compression, source origins and SHA-256. `Dataset_index.csv` gives the package roles. For research reuse, cite the repository commit and accepted September 6 lineage and consult the study team for complete upstream sources.
