# Selected-data dictionary

All CSVs are UTF-8 and use decimal numeric values, not locale-formatted strings. Percent-change fields are **percentage points as numbers**: `47.87` represents `47.87%`, not a fraction of `0.4787`. Rows are complete published result summaries, not individual organisms or in situ concentration observations.

## Shared identifiers

| Field | Meaning |
| --- | --- |
| `scenario` in Figure 4 | `ssp126_2050`, `ssp245_2050`, `ssp585_2050`: climate habitat scenarios for the 2045–2055 mean |
| `scenario` in Table S5 | `SSP1-2.6`, `SSP2-4.5`, `SSP5-8.5`, mapped explicitly to the Figure 4 identifiers |
| `region` | Global, Caribbean, GBR, Southeast Asia / Southeast_Asia. Global includes the other reporting regions; never sum these rows |
| `member` | The exact archived CMIP6 member identifier, retained without relabelling source data |
| `_mean`, `_min`, `_max` | Arithmetic member mean and minimum/maximum across three climate-member results |

## Future MPEI exposure tables

`figure4_period_matched_mpei_change_by_member_20260902.csv` has 36 member × scenario × region rows. The ensemble table has 12 rows.

| Field family | Definition / unit |
| --- | --- |
| `historical_mpei_burden_kg_m3_km2` | Historical integral of HSI × MPEI × cell area; kg m⁻³ km² |
| `future_mpei_burden_kg_m3_km2` | Corresponding future integral in future OR10 co-suitable habitat |
| `historical_mpei_pressure_kg_m3` | Historical burden / HSI-weighted habitat area; kg m⁻³ |
| `future_mpei_pressure_kg_m3` | Future burden / HSI-weighted habitat area; kg m⁻³ |
| `dynamic_mpei_burden_change_pct` | Member-wise percentage change relative to historical burden |
| `dynamic_mpei_pressure_change_pct` | Member-wise percentage change relative to historical pressure |
| `future_habitat_area_km2` | Future OR10 co-suitable habitat area; km² |
| `future_hsi_weighted_area_km2` | Future area integral of HSI; km² × HSI |

Ensemble suffixes extend these same field families. The burden unit is an exposure-index area integral, not a measured plastic mass flux.

## Continuous-quality tables

`figure4_quality_reduction_by_member_20260902.csv` has 36 rows; its ensemble counterpart has 12.

- `historical_continuous_quality_km2_hsi`: historical fixed-domain integral of continuous four-taxon geometric-mean suitability.
- `future_continuous_quality_km2_hsi`: future quality integrated on that same historical domain.
- `continuous_quality_reduction_pct`: `100 × (1 − future / historical)`.
- `continuous_quality_retention_pct`: `100 × future / historical`.
- Ensemble suffixes retain the mean and member extrema. This is not coral cover or a risk probability.

## Repaired conservation table

`Supplementary_Table_S5_20260906.csv` has 12 region × scenario rows.

| Field | Definition |
| --- | --- |
| `future_mpei_vs_climate_area_changed_pct` | Selected priority area reallocated when MPEI is used as a secondary criterion, relative to habitat-only selection |
| `area_changed_lower_95`, `area_changed_upper_95` | Published 95% percentile bootstrap limits, percentage points |
| `future_mpei_vs_climate_mean_mpei_change_pct` | Percentage change in selected reef-area and model-support-weighted mean MPEI |
| `mpei_change_lower_95`, `mpei_change_upper_95` | Published 95% percentile bootstrap limits for that change |

Equal selected reef area and model-support tiers are retained between the two configurations. Model support is the fraction of three models satisfying four-taxon OR10 and MESS criteria; it is not a selection probability. Small negative limits close to 1e-14 are numerical roundoff around zero, preserved in the source.

## Historical heatwave table

`regional_mhw_tercile_plot_data_20260807.csv` has 60 region × response × group rows.

| Field | Definition |
| --- | --- |
| `response` | Acropora, Lobophora, Scarus, Cephalopholis or community |
| `mhw_group` / `group_order` | Low/reference (0), middle (1), high (2) annual cumulative MHW exposure groups |
| `mean_unit_hsi_weighted_mpei_index` | Group mean unit HSI-weighted static 54-tracer MPEI; kg m⁻³ |
| `contrast_vs_low_pct` | Relative difference versus low-MHW group |
| `ci95_low_pct`, `ci95_high_pct` | Published 95% moving-block-bootstrap interval |
| `ci95_excludes_zero` | Parent analysis interval classification, retained in raw CSV |

This table uses static 1993–2019 MPEI and 1993–2022 HSI/MHW grouping. It is not the future 36-tracer field.

## Spatial data and browser schema

The three `ssp*_2050_ensemble_mean_quality_reduction_pct_20260902.tif` files are original accepted derived EPSG:4326 rasters. NoData cells are masked; they are not treated as zero. Browser `map_*.json` files contain one-degree display-bin means and counts. Browser `demo.json` holds normalized future, historical and conservation records plus their source manifest. Missing values are rejected by processing; no fictitious observations are inserted.

`data/provenance.json` records source-relative paths, SHA-256, bytes and CSV row/column coverage. Original machine drive letters, credentials and unrelated project files are excluded from the public subset.
