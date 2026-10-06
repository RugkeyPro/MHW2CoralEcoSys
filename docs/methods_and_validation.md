# Methods and verification boundary

## Accepted analyses

The authoritative research decision uses one global 0.083° SDM for four functional taxa: *Acropora*, *Lobophora*, *Scarus* and *Cephalopholis*. Caribbean, GBR and Southeast Asia reporting products are crops of global outputs, not independently trained regional models.

Future habitat projections use GFDL-ESM4, CNRM-ESM2-1 and IPSL-CM6A-LR with SSP126/245/585. The original Figure 4 methods define continuous co-suitability as the geometric mean of constrained taxon HSI values. Continuous quality is integrated on the fixed historical OR10 domain. Exposure instead pairs historical and future OR10 habitat with period-matched MPEI. That difference in domains is preserved; the demo does not equate the quality and exposure denominators.

## Numerical example

For each member, region and scenario:

```text
quality_reduction_pct = 100 × (1 − future_quality / historical_quality)
burden_change_pct = 100 × (future_burden / historical_burden − 1)
pressure_change_pct = 100 × (future_pressure / historical_pressure − 1)
future_pressure = future_burden / future_HSI_weighted_habitat_area
```

Member changes are averaged *after* these percentage calculations. The reported ensemble means are not ratios of ensemble totals. Means and extrema are reconciled to the accepted ensemble tables with absolute tolerance 1e-8 and relative tolerance 1e-10. No region is aggregated into another, and no historical record is mixed with the future exposure field.

The nine selected file hashes are verified before numerical processing. No source values are corrected, imputed, winsorized or tuned to obtain a desired result. Unexpected missing/duplicate scenario–region keys, missing climate members, nonfinite numbers and nonpositive historical denominators raise errors.

## Uncertainty

- Future Figure 4 whiskers: minimum and maximum of the three climate-member results, not a confidence interval.
- Historical MHW contrasts: original 95% intervals from 2,000 circular four-year moving-block bootstrap replicates; three groups of 10 years. The published summary is included, not the complete original bootstrap draws. The demo displays these intervals and does not regenerate those draws.
- Conservation: full-data estimates with published 95% percentile intervals from the September 6 repaired result. Joint climate-model and 5° spatial-block resampling used 1,000 nominal replicates; 991 effective in GBR. The demo validates summaries and displays intervals; it does not rerun the full raster bootstrap.

The global conservation gain interval includes zero. An absolute tolerance of 1e-10 percentage points is used only to recognize signed near-zero floating-point quantiles in the UI. Source values remain unchanged.

## Map display

All three included GeoTIFFs have EPSG:4326 and shape 1934 × 4338. Each has 76,635 finite, valid analysis cells. Negative quality-reduction values indicate increased quality and are preserved.

`scripts/build_map_subset.py` bins valid cell-center coordinates into one-degree longitude/latitude squares and averages their unweighted source values. Every original valid cell belongs to a displayed bin. Each displayed record is `[longitude_center, latitude_center, arithmetic_mean_pct, source_cell_count]`. The map uses a Natural Earth projection and public-domain Natural Earth land geometry. Region labels are navigation anchors; no fabricated reef polygons are shown.

These simple display-bin means are not area-weighted ecological summaries. Headlines are always obtained from the accepted area-integrated result tables. Changing the climate member affects numeric member cards; the map remains explicitly labelled an ensemble-mean raster.

## What was checked

- Source-copy SHA-256 equality and current bundled source hashes.
- 252 comparisons of member-wise transformations, ensemble means/extrema and burden/pressure identity.
- Frozen browser dataset equals the newly assembled numerical output.
- Distinct historical keys and complete 60-record source coverage.
- Map display coverage and valid coordinate ranges; negative values retained.
- Frontend value formatting, real regional/scenario switches, export escaping and near-zero interval handling.
- Browser interactions, download events, English/Chinese switching and mobile overflow.
- Production build and visual inspection of desktop, English, conservation and mobile screenshots.

The local validation report records actual commands and outcomes. GitHub CI independently reruns the numerical and browser checks after the proposed commit is uploaded.

## Limits inherited from the research archive

The future MPEI is an arithmetic mean of annual instantaneous endpoint fields, not within-year concentration averages. It is a model source trajectory shared by the SSPs; physical forcing beyond its native period was cycled in the source archive. It is not an independent future ocean reanalysis or SSP-specific plastic emissions forecast.

Historical MHW comparisons retain SST/year confounding. They cannot independently attribute all exposure differences to marine heatwaves. External MPEI does not represent internal dose or toxicological outcome. Selected-area MPEI increase under exposure-aware ranking is not pollution removal.

The full project's original manuscript, SEM provenance, environment reconstruction and three missing migration archives are not closed by this demo. No complete MaxEnt refitting, NetCDF transformation or independent scientific validation of all parent inputs was performed.
